using IndexFunArrays # for gaussian blob generation
using Random # to seed the random generator for reproducibility
using AstroStacker
using Statistics
using Interpolations # for the synthetic local warp used to build test frames

Random.seed!(7)

# a simple, dependency-free box blur used only to simulate spatially-varying (seeing) blur below
function box_blur(img, radius)
    out = zeros(size(img))
    n = 0
    for dx in -radius:radius, dy in -radius:radius
        out .+= circshift(img, (dx, dy))
        n += 1
    end
    return out ./ n
end

# a smooth, non-rigid per-frame displacement field (like a mild atmospheric-seeing distortion),
# applied via sub-pixel bilinear resampling
function local_warp(img, A, L, phase)
    sz = size(img)
    itp = extrapolate(interpolate(img, BSpline(Linear())), 0.0)
    out = similar(img)
    for x in 1:sz[1], y in 1:sz[2]
        dx = A * sin(2pi * x / L + phase[1]) * cos(2pi * y / L + phase[2])
        dy = A * cos(2pi * x / L + phase[3]) * sin(2pi * y / L + phase[4])
        out[x, y] = itp(x + dx, y + dy)
    end
    return out
end

function local_blur_mask(sz, center, radius)
    xs = 1:sz[1]
    ys = 1:sz[2]
    d2 = (xs .- center[1]) .^ 2 .+ (ys .- center[2]) .^ 2
    return exp.(-d2 ./ (2 * radius^2))
end

@testset "patch_quality_grid and quality_weight_map discriminate sharp from blurred regions" begin
    # a direct, registration-independent check of the core new mechanism: a patch quality metric
    # that responds to local sharpness, and a weight map that is correspondingly higher where the
    # image is locally sharp.
    sz = (96, 96)
    N_blobs = 90
    positions = sz .* rand(2, N_blobs)
    amps = 0.5 .+ rand(N_blobs)
    sigmas = 2.0 .+ 2.0 .* rand(2, N_blobs)
    frame_sharp = gaussian(sz, offset=positions, weight=amps, sigma=sigmas)

    half = sz[2] ÷ 2
    frame_half_blurred = copy(frame_sharp)
    frame_half_blurred[:, 1:half] .= box_blur(frame_sharp[:, 1:half], 14)

    q = patch_quality_grid(frame_half_blurred, (4, 4))
    @test size(q) == (4, 4)
    @test mean(q[:, 1:2]) < mean(q[:, 3:4]) # left (blurred) half has lower quality than the sharp right half

    nodes = AstroStacker.FindShift.get_default_markers(frame_half_blurred, (4, 4))
    w = quality_weight_map(q, nodes, sz; power=2.0)
    @test size(w) == sz
    @test mean(w[:, half+1:end]) > 5 * mean(w[:, 1:half]) # sharp half gets substantially more weight
end

@testset "test_stack_many_lucky" begin
    sz = (96, 96)
    N_blobs = 90
    positions = sz .* rand(2, N_blobs)
    amps = 0.5 .+ rand(N_blobs)
    sigmas = 2.0 .+ 2.0 .* rand(2, N_blobs)
    ground_truth = gaussian(sz, offset=positions, weight=amps, sigma=sigmas)

    # build a series of frames: each has a genuine smooth local (non-rigid) distortion -- the kind of
    # thing patch-wise TPS registration is meant to correct -- plus one randomly-located, badly blurred
    # ("bad seeing") patch, simulating spatially- and temporally-varying atmospheric seeing.
    Nframes = 8
    frames = Vector{Matrix{Float64}}(undef, Nframes)
    for n in 1:Nframes
        phase = 2pi .* rand(4)
        warped = local_warp(ground_truth, 1.5, 40.0, phase)
        blurred = box_blur(warped, 14)
        center = sz .* (0.2 .+ 0.6 .* rand(2))
        mask = local_blur_mask(sz, center, 14.0)
        frame = mask .* blurred .+ (1 .- mask) .* warped
        frame .+= 0.01 .* randn(sz...)
        frames[n] = frame
    end
    input_stack = cat(frames...; dims=3)

    out = stack_many_lucky(input_stack; grid_size=(6, 6), verbose=false)
    @test size(out.result) == sz
    @test length(out.aligned) == Nframes
    @test length(out.quality) == Nframes

    # isolate the effect of the quality-weighted blending itself (as opposed to the registration, which
    # both share) by comparing against a plain unweighted average of the very same aligned frames.
    naive_mean = dropdims(mean(cat(out.aligned...; dims=3), dims=3), dims=3)

    rms(a, b) = sqrt(mean(abs2.(a .- b)))
    lucky_err = rms(out.result, ground_truth)
    naive_err = rms(naive_mean, ground_truth)

    @test lucky_err < naive_err
end
