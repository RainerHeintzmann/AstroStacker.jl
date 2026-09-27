using IndexFunArrays # for gaussian blob generation
using Random # to seed the random generator for reproducibility
using AstroStacker
using Statistics

Random.seed!(3)

@testset "stack_many_fft (mono, no drizzle) recovers a known translation" begin
    sz = (80, 80)
    N_blobs = 60
    positions = sz .* rand(2, N_blobs)
    amps = 0.5 .+ rand(N_blobs)
    sigmas = 1.5 .+ 1.5 .* rand(2, N_blobs)
    ground_truth = gaussian(sz, offset=positions, weight=amps, sigma=sigmas)

    # simple global (integer-pixel) per-frame jitter, plus a little noise -- exactly what
    # stack_many_fft's translation-only FFT registration is meant to correct. Frame 1 (ref_slice below)
    # is pinned at zero shift, since the stacked result is aligned onto *that* frame's own coordinate
    # system, not onto ground_truth directly -- comparing against ground_truth is only meaningful if the
    # reference frame itself has no offset from it.
    Nframes = 6
    frames = Vector{Matrix{Float64}}(undef, Nframes)
    frames[1] = ground_truth .+ 0.01 .* randn(sz...)
    for n in 2:Nframes
        frames[n] = circshift(ground_truth, (rand(-3:3), rand(-3:3))) .+ 0.01 .* randn(sz...)
    end
    input_stack = cat(frames...; dims=3)

    result, all_params = stack_many_fft(input_stack; ref_slice=1, use_drizzle=false, verbose=false)
    @test length(all_params) == Nframes
    result_2d = dropdims(result, dims=(3, 4))
    @test size(result_2d) == sz

    rms(a, b) = sqrt(mean(abs2.(a .- b)))
    raw_mean = dropdims(mean(input_stack, dims=3), dims=3)

    # registration should clearly beat a naive average of the un-registered, jittered frames
    @test rms(result_2d, ground_truth) < rms(raw_mean, ground_truth)
end

@testset "stack_many_fft (Bayer/drizzle path) runs via fft_find_transform" begin
    # bayer_pattern="RGGB" (explicit -- bayer_pattern defaults to `nothing`, i.e. mono/no debayering)
    # treats input_stack as a Bayer-pattern mosaic and runs the same drizzle pipeline as stack_many, just
    # with translation-only FFT registration instead of star detection + RANSAC. There's no
    # numerically-meaningful "ground truth RGB scene" comparison without a true Bayer-consistent fixture,
    # so this checks that the pipeline (registration -> drizzle warp -> outlier rejection) runs correctly
    # end-to-end and produces the expected shape, for both supported shift estimators.
    sz = (64, 64)
    N_blobs = 40
    positions = sz .* rand(2, N_blobs)
    amps = 0.5 .+ rand(N_blobs)
    sigmas = 1.5 .+ 1.5 .* rand(2, N_blobs)
    ground_truth = gaussian(sz, offset=positions, weight=amps, sigma=sigmas)

    Nframes = 5
    frames = [circshift(ground_truth, (rand(-2:2), rand(-2:2))) .+ 0.01 .* randn(sz...) for _ in 1:Nframes]
    input_stack = cat(frames...; dims=3)

    result, all_params = stack_many_fft(input_stack; ref_slice=1, verbose=false, bayer_pattern="RGGB")
    @test length(all_params) == Nframes
    @test size(result) == (sz..., 1, 3)

    result_lk, all_params_lk = stack_many_fft(input_stack; ref_slice=1, verbose=false, bayer_pattern="RGGB", shift_fun=AstroStacker.FindShift.find_shift_lk)
    @test length(all_params_lk) == Nframes
    @test size(result_lk) == (sz..., 1, 3)
end

@testset "stack_many_fft (Bayer/drizzle path, drizzle_supersampling=1) doesn't crash in remove_outliers" begin
    # Regression test: with use_drizzle=true and drizzle_supersampling exactly 1, all_masks used to stay
    # a tiny (1,1,Nimgs,1) dummy placeholder (only drizzle_supersampling != 1 triggered a real per-pixel
    # mask), causing a DimensionMismatch crash in remove_outliers regardless of how many outliers were
    # found. The fix gates real mask allocation on use_drizzle alone.
    sz = (64, 64)
    N_blobs = 40
    positions = sz .* rand(2, N_blobs)
    amps = 0.5 .+ rand(N_blobs)
    sigmas = 1.5 .+ 1.5 .* rand(2, N_blobs)
    ground_truth = gaussian(sz, offset=positions, weight=amps, sigma=sigmas)

    Nframes = 5
    frames = [circshift(ground_truth, (rand(-2:2), rand(-2:2))) .+ 0.01 .* randn(sz...) for _ in 1:Nframes]
    input_stack = cat(frames...; dims=3)

    result, all_params = stack_many_fft(input_stack; ref_slice=1, verbose=false, bayer_pattern="RGGB", drizzle_supersampling=1.0)
    @test length(all_params) == Nframes
    # drizzle_supersampling=1 keeps the output at the Bayer-subsampled (half) resolution, unlike the
    # default 2.0 (which restores full resolution).
    @test size(result) == (sz .÷ 2..., 1, 3)
end

@testset "stack_many_fft (mono, drizzle) accumulates via scatter-warp without debayering" begin
    # bayer_pattern=nothing (the default) + use_drizzle=true (also the default): mono/already-color data
    # gets the drizzle scatter-warp treatment (sub-pixel accumulation, optional supersampling) without any
    # Bayer sub-image splitting or 3-channel reconstruction.
    sz = (64, 64)
    N_blobs = 40
    positions = sz .* rand(2, N_blobs)
    amps = 0.5 .+ rand(N_blobs)
    sigmas = 1.5 .+ 1.5 .* rand(2, N_blobs)
    ground_truth = gaussian(sz, offset=positions, weight=amps, sigma=sigmas)

    Nframes = 6
    frames = Vector{Matrix{Float64}}(undef, Nframes)
    frames[1] = ground_truth .+ 0.01 .* randn(sz...)
    for n in 2:Nframes
        frames[n] = circshift(ground_truth, (rand(-2:2), rand(-2:2))) .+ 0.01 .* randn(sz...)
    end
    input_stack = cat(frames...; dims=3)

    # drizzle_supersampling=1 here: unlike a Bayer mosaic (whose 4 sub-pixel positions give built-in
    # coverage diversity), mono drizzle relies entirely on registration-derived sub-pixel shifts to fill a
    # supersampled grid without holes -- these synthetic frames use purely integer-pixel shifts, so
    # supersampling>1 would (correctly) trigger the max_uncovered_frac warning rather than reflect a bug.
    # Correctness of the scatter-warp accumulation itself is fully exercised at supersampling=1.
    result, all_params = stack_many_fft(input_stack; ref_slice=1, verbose=false, drizzle_supersampling=1.0)
    @test length(all_params) == Nframes
    @test size(result) == (sz..., 1, 1)

    result_2d = dropdims(result, dims=(3, 4))
    rms(a, b) = sqrt(mean(abs2.(a .- b)))
    raw_mean = dropdims(mean(input_stack, dims=3), dims=3)
    @test rms(result_2d, ground_truth) < rms(raw_mean, ground_truth)

    # supersampling=2 should at least run and produce the expected (upsampled) shape.
    result_super, _ = stack_many_fft(input_stack; ref_slice=1, verbose=false, drizzle_supersampling=2.0)
    @test size(result_super) == (2 .* sz..., 1, 1)
end
