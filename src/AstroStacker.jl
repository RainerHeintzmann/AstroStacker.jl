module AstroStacker

using Astroalign
using CoordinateTransformations # Translation
using StaticArrays
using Statistics # for median

export correct_dark_flat, stack_many
export com_psf, collect_info
export bin_mono, bin_rgb

include("utils.jl")
include("preprocess_helpers.jl")
include("findpeaks.jl")
include("warp.jl")
include("stacker.jl")

end # module AstroStacker
