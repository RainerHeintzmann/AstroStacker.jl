"""
    laplacian_variance(patch)

local sharpness ("focus") metric: the variance of the discrete Laplacian of `patch`.
Higher values indicate locally sharper (more high-spatial-frequency content) data; out-of-focus or
badly-seeing-blurred patches have a low Laplacian variance. This is the standard sharpness metric used
by lucky-imaging tools (e.g. AutoStakkert!, Registax) to rank frames/regions.
"""
function laplacian_variance(patch::AbstractMatrix)
    lap = 4 .* patch .- circshift(patch, (1, 0)) .- circshift(patch, (-1, 0)) .- circshift(patch, (0, 1)) .- circshift(patch, (0, -1))
    return var(lap)
end

"""
    patch_quality_grid(img, grid_size; patch_size=max.(size(img) .÷ grid_size, 5))

computes the [`laplacian_variance`](@ref) quality metric on a regular `grid_size` grid of patches over `img`,
using the same patch layout (`FindShift.extract_patches`) as `FindShift.find_deformations`, so that the
resulting quality grid lines up with the deformation nodes.

Returns a `grid_size` array of quality values.
"""
function patch_quality_grid(img, grid_size; patch_size=max.(size(img) .÷ grid_size, 5))
    # avoid_border=false: we need one quality value per full grid node (matching FindShift.get_default_markers),
    # not the border-omitting reduced grid that FindShift.extract_patches uses by default.
    patches, _ = FindShift.extract_patches(img, patch_size=patch_size; grid_size=grid_size, avoid_border=false)
    q = zeros(Float64, grid_size)
    for (n, ci) in enumerate(CartesianIndices(grid_size))
        q[ci] = laplacian_variance(patches[n])
    end
    return q
end

"""
    quality_weight_map(q_grid, nodes, img_size; power=1.0)

turns the sparse per-patch quality values `q_grid` (as returned by [`patch_quality_grid`](@ref), on the
`nodes` grid returned by `FindShift.get_default_markers`) into a smooth, dense per-pixel weight map of
size `img_size`, via bilinear (`Gridded(Linear())`) interpolation between patch centers. Since `nodes`
spans the full image extent (first to last pixel along each dimension), no extrapolation is needed.

`power` is applied to the (clamped non-negative) interpolated weight to control how aggressively the
blend favors the locally sharpest frame (`power=1` is a soft/linear blend, larger values approach hard
best-frame selection).
"""
function quality_weight_map(q_grid::AbstractMatrix, nodes, img_size; power=1.0)
    itp = interpolate(nodes, q_grid, Gridded(Linear()))
    w = [itp(x, y) for x in 1:img_size[1], y in 1:img_size[2]]
    return max.(w, 0) .^ power
end
