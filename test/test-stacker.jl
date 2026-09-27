using IndexFunArrays # for gaussian blob generation
using Random # to seed the random generator for reproducibility
using AstroStacker
using Astroalign

Random.seed!(42)
# using View5D  # recommended for visualization
@testset "test_stacker" begin
    # create inital star coordinates
    sz = (100, 100)
    mid_pos = [sz...] ./ 2
    N = 60
    star_pos = sz .* rand(2, N)
    star_amp = rand(N)
    star_shape = 0.3 .+ rand(2, N)
    # create star image
    y1 = gaussian(sz, offset=star_pos, weight=star_amp, sigma=star_shape)
    rot_mat(alpha) = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]
    # rotate the stars by random angles
    frames = []
    coords = []
    N = 10
    for alpha in 33 * rand(N)
        shift_vec = 10.0 .* (rand(2) .- 0.5)
        new_pos = rot_mat(alpha*pi/180) * (star_pos .- mid_pos) .+ mid_pos .+ shift_vec
        frame = gaussian(sz, offset=new_pos, weight=star_amp, sigma=star_shape)
        push!(frames, frame)
        push!(coords, new_pos)
    end
    frames = cat(frames..., dims=3)
    # @vt frames # visualize all frames as a 3d stack
    # align the first frame in the stack to the reference image
    y2_aligned = Astroalign.align_frames([frames[:,:,1]], y1)[1];
    # @vt y1 y2_aligned # display alignment (toggle between frames using the keys `,` and `.`)
    # check alignment by summing at the origninal source positions and the destinatoin positions
    y2_aligned[isnan.(y2_aligned)] .= 0
    y1_coords = clamp.(round.(Int, star_pos), 1, sz[1]-1)
    y2_coords = clamp.(round.(Int, coords[1]), 1, sz[1]-1)
    sum_at_true = sum(sum.([y2_aligned[c...] for c in eachslice(y1_coords, dims=2)]))
    sum_at_source =sum(sum.([y2_aligned[c...] for c in eachslice(y2_coords, dims=2)]))
    @test sum_at_true / sum_at_source > 8
    # stack all frames in the stack to the coordinates of the first frame
    f = com_psf;
    stacked, all_params = stack_many(frames; ref_slice=1, f=f, use_drizzle=false, verbose=false, box_size= (15,15));
    @test size(stacked) == (100, 100, 1, 1)
    # @vt res
    sum_at_true = sum(sum.([stacked[c...] for c in eachslice(y2_coords, dims=2)]))
    sum_at_source =sum(sum.([stacked[c...] for c in eachslice(y1_coords, dims=2)]))
    @test sum_at_true / sum_at_source > 8
end

@testset "stack_many fast/slow-mem tiering hooks (to_fast_mem/to_slow_mem)" begin
    # `copy` is a non-identity but still-CPU "fake tiering" pair: it exercises the staging/destaging
    # plumbing added to do_drizzle_warp! without needing real GPU hardware. Since it doesn't change any
    # values, stacking with it must reproduce exactly the same result as the default (identity) hooks.
    sz = (100, 100)
    N = 60
    star_pos = sz .* rand(2, N)
    star_amp = rand(N)
    star_shape = 0.3 .+ rand(2, N)
    Nframes = 6
    frames = cat((gaussian(sz, offset=star_pos .+ 3 .* (rand(2, N) .- 0.5), weight=star_amp, sigma=star_shape) for _ in 1:Nframes)..., dims=3)

    Random.seed!(123)
    stacked_default, _ = stack_many(frames; ref_slice=1, f=com_psf, use_drizzle=true, bayer_pattern="RGGB", verbose=false, box_size=(15,15))

    Random.seed!(123)
    stacked_tiered, _ = stack_many(frames; ref_slice=1, f=com_psf, use_drizzle=true, bayer_pattern="RGGB", verbose=false, box_size=(15,15),
                                    to_fast_mem=copy, to_slow_mem=copy)

    @test stacked_default == stacked_tiered
end

@testset "stack_many bayer_pattern + use_drizzle=false pre-debayers then warps per channel" begin
    # This combination debayers the whole stack up front via bin_rgb, then falls through to the ordinary
    # already-color, non-drizzle code path (registration from the mono reference channel, applied
    # independently to each of the 3 resulting color channels). Verify that's exactly what happens by
    # cross-checking against manually calling bin_rgb + stack_many(...; use_drizzle=false) directly.
    sz = (100, 100)
    N = 60
    star_pos = sz .* rand(2, N)
    star_amp = rand(N)
    star_shape = 0.3 .+ rand(2, N)
    Nframes = 5
    frames = cat((gaussian(sz, offset=star_pos .+ 3 .* (rand(2, N) .- 0.5), weight=star_amp, sigma=star_shape) for _ in 1:Nframes)..., dims=3)

    Random.seed!(77)
    result_direct, _ = stack_many(frames; ref_slice=1, f=com_psf, use_drizzle=false, bayer_pattern="RGGB", verbose=false, box_size=(15,15))
    @test size(result_direct) == (sz[1] ÷ 2, sz[2] ÷ 2, 1, 3)

    debayered = bin_rgb(frames; bayer_pattern="RGGB")
    Random.seed!(77)
    result_manual, _ = stack_many(debayered; ref_slice=1, f=com_psf, use_drizzle=false, verbose=false, box_size=(15,15))

    @test result_direct == result_manual
end

@testset "stack_many debayer_fun=debayer_interp keeps full resolution" begin
    sz = (100, 100)
    N = 60
    star_pos = sz .* rand(2, N)
    star_amp = rand(N)
    star_shape = 0.3 .+ rand(2, N)
    Nframes = 5
    frames = cat((gaussian(sz, offset=star_pos .+ 3 .* (rand(2, N) .- 0.5), weight=star_amp, sigma=star_shape) for _ in 1:Nframes)..., dims=3)

    result, _ = stack_many(frames; ref_slice=1, f=com_psf, use_drizzle=false, bayer_pattern="RGGB",
                            debayer_fun=debayer_interp, verbose=false, box_size=(15,15))
    # unlike the default debayer_fun=bin_rgb (which halves resolution), debayer_interp keeps the full
    # original (H, W) resolution.
    @test size(result) == (sz[1], sz[2], 1, 3)

    # cross-check equivalence with calling debayer_interp + stack_many(...; use_drizzle=false) manually.
    debayered = debayer_interp(frames; bayer_pattern="RGGB")
    Random.seed!(88)
    result_direct, _ = stack_many(frames; ref_slice=1, f=com_psf, use_drizzle=false, bayer_pattern="RGGB",
                                   debayer_fun=debayer_interp, verbose=false, box_size=(15,15))
    Random.seed!(88)
    result_manual, _ = stack_many(debayered; ref_slice=1, f=com_psf, use_drizzle=false, verbose=false, box_size=(15,15))
    @test result_direct == result_manual
end

@testset "stack_many errors on drizzle_supersampling != 1 with use_drizzle=false" begin
    sz = (60, 60)
    N = 30
    star_pos = sz .* rand(2, N)
    star_amp = rand(N)
    star_shape = 0.3 .+ rand(2, N)
    Nframes = 3
    frames = cat((gaussian(sz, offset=star_pos, weight=star_amp, sigma=star_shape) for _ in 1:Nframes)..., dims=3)

    # explicitly requesting supersampling without the drizzle algorithm is a likely mistake (it would
    # otherwise be silently ignored) and must error, regardless of bayer_pattern. This must error before
    # any registration/photometry is attempted, so plain content (not necessarily star-like) is fine here.
    @test_throws ErrorException stack_many(frames; use_drizzle=false, drizzle_supersampling=2.0, verbose=false)
    @test_throws ErrorException stack_many(frames; use_drizzle=false, bayer_pattern="RGGB", drizzle_supersampling=2.0, verbose=false)

    # leaving drizzle_supersampling untouched must NOT error -- its default depends on use_drizzle.
    result, _ = stack_many(frames; ref_slice=1, f=com_psf, use_drizzle=false, verbose=false, box_size=(15,15))
    @test size(result)[1:2] == sz
end