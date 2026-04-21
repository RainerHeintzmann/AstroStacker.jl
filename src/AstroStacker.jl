module AstroStacker

using Astroalign
using CoordinateTransformations # Translation
using StaticArrays

export correct_dark_flat, stack_many
export com_psf

include("utils.jl")
include("preprocess_helpers.jl")
include("findpeaks.jl")
include("warp.jl")
include("stacker.jl")

end # module AstroStacker
