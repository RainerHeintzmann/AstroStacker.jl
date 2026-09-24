        using AstroStacker
using Test

# args = parse_args(Base.ARGS)
# testsuite = find_tests(@__DIR__)

# runtests(Astroalign, args; testsuite, init_code)

include("test-warp-kernel.jl")
include("test-stacker.jl")
include("test-lucky-stack.jl")
include("test-fft-stack.jl")
