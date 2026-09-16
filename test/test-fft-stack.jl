using IndexFunArrays # for gaussian blob generation
using Random # to seed the random generator for reproducibility
using AstroStacker
using Statistics

Random.seed!(3)

@testset "test_stack_many_fft" begin
    sz = (80, 80)
    N_blobs = 60
    positions = sz .* rand(2, N_blobs)
    amps = 0.5 .+ rand(N_blobs)
    sigmas = 1.5 .+ 1.5 .* rand(2, N_blobs)
    ground_truth = gaussian(sz, offset=positions, weight=amps, sigma=sigmas)

    # simple global (integer-pixel) per-frame jitter, plus a little noise -- exactly what
    # stack_many_fft's translation-only FFT registration is meant to correct.
    Nframes = 8
    frames = Vector{Matrix{Float64}}(undef, Nframes)
    for n in 1:Nframes
        frames[n] = circshift(ground_truth, (rand(-3:3), rand(-3:3))) .+ 0.01 .* randn(sz...)
    end
    input_stack = cat(frames...; dims=3)

    out = stack_many_fft(input_stack; verbose=false)
    @test size(out.result) == sz
    @test size(out.aligned) == size(input_stack)
    @test length(out.shifts) == Nframes

    rms(a, b) = sqrt(mean(abs2.(a .- b)))
    raw_mean = dropdims(mean(input_stack, dims=3), dims=3)

    # registration should clearly beat a naive average of the un-registered, jittered frames
    @test rms(out.result, ground_truth) < rms(raw_mean, ground_truth)
end
