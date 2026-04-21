# Helper function. Transforms source coordinates and copies values.
function warp_assign!(result, weights, src, x,y, tfm, supersample)
    src_pos = SVector{2}(x,y) # [x, y]
    # Forward transform source coord to dest coord
    dest_coord = tfm(src_pos)

    # Map to supersampled grid
    px = round(Int, dest_coord[1] * supersample)
    py = round(Int, dest_coord[2] * supersample)

    if checkbounds(Bool, result, px, py)
        @inbounds result[px, py] += src[x, y]
        @inbounds weights[px, py] += one(eltype(weights))
    end
end

# the same but with linear interpolation
function warp_assign_interp!(result, weights, src, x,y, tfm, supersample)
    src_pos = SVector{2}(x,y) # [x, y]
    # Forward transform source coord to dest coord
    dest_coord = tfm(src_pos)

    fx = dest_coord[1] * supersample
    fy = dest_coord[2] * supersample
    # Map to supersampled grid
    px = floor(Int, fx)
    wx = fx .- px
    py = floor(Int, fy)
    wy = fy .- py

    s = src[x, y]
    if checkbounds(Bool, result, px, py)
        w = (1-wx)*(1-wy)
        @inbounds result[px, py] += w*s
        @inbounds weights[px, py] += w
    end
    if checkbounds(Bool, result, px+1, py)
        w = wx*(1-wy)
        @inbounds result[px+1, py] += w*s
        @inbounds weights[px+1, py] += w
    end
    if checkbounds(Bool, result, px, py+1)
        w = (1-wx)*wy
        @inbounds result[px, py+1] += w*s
        @inbounds weights[px, py+1] += w
    end
    if checkbounds(Bool, result, px+1, py+1)
        w = wx*wy
        @inbounds result[px+1, py+1] += w*s
        @inbounds weights[px+1, py+1] += w
    end
end

"""
    forward_warp!(result, weights, src, tfm, dest_size; supersample = 1)

Forward-mode warp: iterate over source, scatter to destination with accumulation.
"""
function forward_warp!(result, weights, src::AbstractMatrix{T}, tfm; use_interp=false, supersample = 1) where T
    if (use_interp)
        warp_assign_interp!.(Ref(result), Ref(weights), Ref(src), axes(src, 1), transpose(axes(src, 2)), Ref(tfm), supersample)
    else
        warp_assign!.(Ref(result), Ref(weights), Ref(src), axes(src, 1), transpose(axes(src, 2)), Ref(tfm), supersample)
    end
    return result, weights
end

"""
    forward_warp(src, tfm, dest_size; supersample = 1)

Forward-mode warp: iterate over source, scatter to destination with accumulation.
"""
function forward_warp(src::AbstractMatrix{T}, tfm, dest_size; use_interp=false, supersample = 1) where T
    out_H, out_W = dest_size .* supersample
    result = zeros(eltype(src), out_H, out_W)
    weights = zeros(eltype(src), out_H, out_W)
    forward_warp!(result, weights, src, tfm; use_interp=use_interp, supersample)
end

"""
    drizzle_warp!(result, drizzle_mask, bayer_mosaic, tfm; use_interp=false, supersample = 2.0, bayer_pattern = "RGGB")

Performs the forward warping of in input bayer mosaic (`bayer_mosaic`) with the transformation as defined by `tfm`, but originally computed on the gridded data (i.e. the top left 4 pixels forming pixel 1).

# Parameters

* `result`:  A necessary output array into which the results are added. Outside pixels are ignored.
* `drizzle_mask`: An output array into which the value one is added at assigned pixel locations.
* `bayer_mosaic`: The input bayer-patter image to grid onto (add into) a color output
* `tfm`: The transformation, but calculated on the 2x2 binned data.
* `use_interp`: if `true` linar interpolation will be used on destination.
* `supersample`: The factor to supersample. The default of 2 means that the output size corresponds to the input size.
* `bayer_pattern`: The order of the pixels in the bayer pattern. Allowed tags are R,G and B.

"""
function drizzle_warp!(result, drizzle_mask, bayer_mosaic, inv_tfm; supersample = 2.0, use_interp=false, bayer_pattern = "RGGB")
    bayer_index = get_bayer_index(bayer_pattern)
    sindex_x = (1, 2, 1, 2)
    sindex_y = (1, 1, 2, 2)
    for bayer_pix in 1:4
        sx = sindex_x[bayer_pix] # Determines the offsets
        sy = sindex_y[bayer_pix]
        src_mat = @view bayer_mosaic[sx:2:end, sy:2:end]
        dst_mat = @view result[:,:,bayer_index[bayer_pix]]
        dst_mask_mat = @view drizzle_mask[:, :, bayer_index[bayer_pix]]
        my_src_shift = Translation(sx - 2, sy - 2)
        my_zoom = AffineMap([supersample 0; 0 supersample],[0, 0])
        tfm_both = compose(my_src_shift, compose(my_zoom, inv(inv_tfm)))
        forward_warp!(dst_mat, dst_mask_mat, src_mat, tfm_both, use_interp=use_interp)
    end 
    return result
end
