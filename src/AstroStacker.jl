module AstroStacker

using Astroalign
using CoordinateTransformations # Translation
using StaticArrays
using Statistics # for median, var, mean
using NDTools
using FindShift
using Interpolations
using KernelAbstractions # for backend-portable (CPU/GPU) drizzle warp kernels

export correct_dark_flat, stack_many, stack_many_lucky, stack_many_fft
export com_psf, collect_info
export bin_mono, bin_rgb
export laplacian_variance, patch_quality_grid, quality_weight_map

include("utils.jl")
include("preprocess_helpers.jl")
include("findpeaks.jl")
include("warp.jl")
include("stacker.jl")
include("quality.jl")
include("lucky_stack.jl")
include("fft_stack.jl")

end # module AstroStacker
