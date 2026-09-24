using AstroStacker
using CoordinateTransformations

@testset "forward_warp! kernel (nearest)" begin
    src = Float64.(reshape(1:12, 3, 4)) # 3x4
    tfm = AffineMap([1.0 0.0; 0.0 1.0], [1.0, 0.0]) # shift by (+1, 0)
    result = zeros(5, 4)
    weights = zeros(5, 4)
    AstroStacker.forward_warp!(result, weights, src, tfm; use_interp=false, supersample=1)

    expected = zeros(5, 4)
    expected[2:4, :] .= src
    @test result == expected
    @test weights[2:4, :] == ones(3, 4)
    @test all(weights[[1, 5], :] .== 0)
end

@testset "forward_warp! kernel (bilinear interp)" begin
    src = ones(4, 4)
    tfm = AffineMap([1.0 0.0; 0.0 1.0], [0.5, 0.5]) # half-pixel shift in both dims
    result = zeros(6, 6)
    weights = zeros(6, 6)
    AstroStacker.forward_warp!(result, weights, src, tfm; use_interp=true, supersample=1)

    # every source pixel splats a total weight of 1 across its 4 surrounding destination pixels,
    # so the total accumulated weight must equal the number of source pixels (mass-conserving).
    @test sum(weights) ≈ length(src)
    @test all(0 .<= weights .<= 1 + 1e-8)
    # a constant source, warped, stays constant wherever fully covered
    @test all(isapprox.(result[2:4, 2:4] ./ weights[2:4, 2:4], 1.0))
end

@testset "drizzle transform stays isbits (required for GPU kernel launch)" begin
    # GPU kernel compilation requires every argument (including the transform) to be isbits; a plain
    # Matrix/Vector-backed AffineMap anywhere in the composition chain silently breaks this. Running only
    # on the CPU KernelAbstractions backend (as the other tests in this file do) can NOT catch that class
    # of bug, since the CPU backend doesn't enforce isbits arguments -- this checks the actual type
    # produced by the production code path directly. Deliberately construct `inv_tfm` from plain
    # (non-static) Matrix/Vector, since nothing guarantees Astroalign.find_transform won't do the same.
    inv_tfm_plain = AffineMap([1.0 0.0; 0.0 1.0], [0.3, -0.7])
    tfm = AstroStacker.bayer_pixel_transform(inv_tfm_plain, 2.0, 1, 1)
    @test isbitstype(typeof(tfm))
end

@testset "forward_warp! kernel matches a naive reference scatter, including collisions" begin
    # a downscaling transform makes several distinct source pixels land on the same destination
    # pixel; the atomic-based kernel must accumulate all of them, like a plain sequential loop would.
    using Random
    Random.seed!(7)
    src = rand(6, 6)
    tfm = AffineMap([0.4 0.0; 0.0 0.4], [1.0, 1.0])
    dst_size = (4, 4)
    result = zeros(dst_size)
    weights = zeros(dst_size)
    AstroStacker.forward_warp!(result, weights, src, tfm; use_interp=false, supersample=1)

    ref_result = zeros(dst_size)
    ref_weights = zeros(dst_size)
    for x in axes(src, 1), y in axes(src, 2)
        px = round(Int, 0.4x + 1.0)
        py = round(Int, 0.4y + 1.0)
        if checkbounds(Bool, ref_result, px, py)
            ref_result[px, py] += src[x, y]
            ref_weights[px, py] += 1
        end
    end
    @test result == ref_result
    @test weights == ref_weights
    @test any(ref_weights .> 1) # sanity: this setup actually produces collisions
end
