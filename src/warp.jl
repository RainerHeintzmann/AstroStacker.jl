# GPU-portable (KernelAbstractions.jl) scatter-accumulate warp kernels.
# Each work-item handles one source pixel (x,y): it maps it to (a) destination pixel(s) via the
# affine transform `tfm` and atomically accumulates into `result`/`weights`, since multiple source
# pixels (and, across drizzle_warp!'s Bayer loop, multiple sub-images) can land on the same
# destination pixel. Atomics make this correct on both the CPU backend and real GPU backends
# (CUDA/AMDGPU/Metal/oneAPI, via whichever such package the caller has loaded) without AstroStacker
# itself depending on any GPU backend package.

@kernel function warp_scatter_kernel!(result, weights, @Const(src), tfm, supersample)
    x, y = @index(Global, NTuple)
    @inbounds s = src[x, y]
    px = round(Int, (tfm.linear[1,1]*x + tfm.linear[1,2]*y + tfm.translation[1]) * supersample)
    py = round(Int, (tfm.linear[2,1]*x + tfm.linear[2,2]*y + tfm.translation[2]) * supersample)
    if checkbounds(Bool, result, px, py)
        KernelAbstractions.@atomic result[px, py] += s
        KernelAbstractions.@atomic weights[px, py] += one(eltype(weights))
    end
end

# the same but with bilinear splatting onto the 4 surrounding destination pixels
@kernel function warp_scatter_interp_kernel!(result, weights, @Const(src), tfm, supersample)
    x, y = @index(Global, NTuple)
    @inbounds s = src[x, y]
    fx = (tfm.linear[1,1]*x + tfm.linear[1,2]*y + tfm.translation[1]) * supersample
    fy = (tfm.linear[2,1]*x + tfm.linear[2,2]*y + tfm.translation[2]) * supersample
    px = floor(Int, fx)
    wx = fx - px
    py = floor(Int, fy)
    wy = fy - py

    if checkbounds(Bool, result, px, py)
        w = (1-wx)*(1-wy)
        KernelAbstractions.@atomic result[px, py] += w*s
        KernelAbstractions.@atomic weights[px, py] += w
    end
    if checkbounds(Bool, result, px+1, py)
        w = wx*(1-wy)
        KernelAbstractions.@atomic result[px+1, py] += w*s
        KernelAbstractions.@atomic weights[px+1, py] += w
    end
    if checkbounds(Bool, result, px, py+1)
        w = (1-wx)*wy
        KernelAbstractions.@atomic result[px, py+1] += w*s
        KernelAbstractions.@atomic weights[px, py+1] += w
    end
    if checkbounds(Bool, result, px+1, py+1)
        w = wx*wy
        KernelAbstractions.@atomic result[px+1, py+1] += w*s
        KernelAbstractions.@atomic weights[px+1, py+1] += w
    end
end

"""
    forward_warp!(result, weights, src, tfm; use_interp=false, supersample = 1, synchronize=true)

Forward-mode warp: iterate over source, scatter to destination with accumulation.

Runs on whichever backend `result`/`weights`/`src` live on (CPU `Array`s, or a GPU array such as a
`CuArray`/`ROCArray`/... provided the corresponding GPU package is loaded), via `KernelAbstractions.jl`.

`synchronize`: if `true` (default), blocks until the kernel has completed before returning, so
`result`/`weights` are immediately ready to read. Callers that launch several such kernels back-to-back
into independent regions of `result` (e.g. [`drizzle_warp!`](@ref), across its 4 Bayer sub-images) can
pass `false` and synchronize once themselves afterwards, avoiding one host/device round-trip per launch.
"""
function forward_warp!(result, weights, src::AbstractMatrix{T}, tfm; use_interp=false, supersample = 1, synchronize=true) where T
    backend = KernelAbstractions.get_backend(result)
    kernel! = use_interp ? warp_scatter_interp_kernel!(backend) : warp_scatter_kernel!(backend)
    kernel!(result, weights, src, tfm, supersample; ndrange=size(src))
    synchronize && KernelAbstractions.synchronize(backend)
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
    bayer_pixel_transform(inv_tfm, supersample, sx, sy)

Composes the per-Bayer-subpixel transform used by [`drizzle_warp!`](@ref) (source-grid shift, then
supersampling zoom, then the inverse of the frame-to-reference registration transform `inv_tfm`).

The result is forced into a static (`StaticArrays`-backed, isbits) `AffineMap`, since it is passed as an
argument to a `KernelAbstractions.jl` kernel in [`forward_warp!`](@ref) -- GPU kernel compilation requires
all arguments to be isbits, and a plain `Matrix`/`Vector`-backed `AffineMap` anywhere in the composition
chain (e.g. from a non-static `inv_tfm`) would otherwise silently poison the whole composed transform.
"""
function bayer_pixel_transform(inv_tfm, supersample, sx, sy)
    my_src_shift = Translation(sx - 2, sy - 2)
    my_zoom = AffineMap(SMatrix{2,2}(supersample, 0, 0, supersample), SVector(0.0, 0.0))
    tfm_both = compose(my_src_shift, compose(my_zoom, inv(inv_tfm)))
    return AffineMap(SMatrix{2,2}(tfm_both.linear), SVector{2}(tfm_both.translation))
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
        tfm_both = bayer_pixel_transform(inv_tfm, supersample, sx, sy)
        # the 4 sub-images write to disjoint regions of result/drizzle_mask, so the 4 kernel launches
        # need no synchronization between them -- only once, after all 4 are queued (below).
        forward_warp!(dst_mat, dst_mask_mat, src_mat, tfm_both, use_interp=use_interp, synchronize=false)
    end
    KernelAbstractions.synchronize(KernelAbstractions.get_backend(result))
    return result
end

"""
    do_drizzle_warp!(drizzle_mask, to_warp, inv_tfm, myaxes, )

an alternative to the `warp` function, to be provided to the alignment function instead of warp.

# Parameters:

- `drizzle_mask`: a mask, collecting which pixel where assigned.
- `result`: the assigned pixels.
- `use_interp`: whether to use interpolation (true) or not.
- `bayer_pattern`: a string of size 4 characters, indicating the order of colors. The default ("RGGB") corresponds to this pattern (starting from the top left corner of `input_stack`):
- `myaxes`: is ignored.
- `to_fast_mem`/`to_slow_mem`: functions to stage this single frame's source/destination onto "fast"
  (e.g. GPU) memory and back to `result`/`drizzle_mask`'s own ("slow") memory, so that only one frame at
  a time -- not the whole multi-frame stack -- needs to live on the fast device. Both default to
  `identity` (no staging; `result`/`drizzle_mask`/`to_warp` are used in place, as before).
- `fast_src_buf`/`fast_result_buf`/`fast_mask_buf`: optional pre-allocated fast-memory scratch buffers
  (one frame's worth each, already staged via `to_fast_mem`), reused via `copyto!` instead of allocating
  fresh ones on every call -- pass these (allocated once outside the per-frame loop) to avoid repeated
  (and, on a GPU, potentially expensive) allocation. If omitted, fresh buffers are allocated as before.
- `slow_src_buf`: an optional pre-allocated plain (concrete `Array`) CPU scratch buffer, the same shape as
  `to_warp`. Required alongside `fast_src_buf` -- `to_warp` is frequently a `SubArray` (or similar) view,
  and GPU array types generally only provide an efficient `copyto!` fast path from/to a concrete `Array`,
  not from/to an arbitrary view/lazy wrapper (which falls back to slow, and on a GPU array disallowed,
  scalar indexing on whichever side is the GPU array); `to_warp` is therefore always first `copyto!`'d
  into this concrete buffer, then that is `copyto!`'d into `fast_src_buf` -- and symmetrically on the way
  back out via `slow_result_buf`/`slow_mask_buf` below.
- `slow_result_buf`/`slow_mask_buf`: like `slow_src_buf`, but for the way back: `result`/`drizzle_mask`
  are themselves frequently `SubArray` views (into the caller's larger accumulator), so `fast_result`/
  `fast_mask` are first `copyto!`'d into these concrete buffers, then those into `result`/`drizzle_mask`.
  Required alongside `fast_result_buf`/`fast_mask_buf`.
"""
function do_drizzle_warp!(drizzle_mask, drizzle_supersampling, bayer_pattern, use_interp, result, to_warp, inv_tfm, myaxes;
                           to_fast_mem=identity, to_slow_mem=identity,
                           fast_src_buf=nothing, fast_result_buf=nothing, fast_mask_buf=nothing,
                           slow_src_buf=nothing, slow_result_buf=nothing, slow_mask_buf=nothing)
        isnothing(to_warp) && error("For drizzle you need to provide a drizzle_supersample! and a to_warp input, the Bayer-pattern mosaic input")
        if to_fast_mem === identity && to_slow_mem === identity
            # No staging requested at all: operate directly in place on whatever device result/
            # drizzle_mask/to_warp already live on -- zero extra buffers/copies, exactly matching
            # behavior from before fast/slow-mem tiering existed. This matters in particular for the
            # "pre-cast the whole stack to a CuArray, don't pass to_fast_mem/to_slow_mem" usage pattern:
            # it must stay fully GPU-resident with no round-trips through CPU staging buffers.
            drizzle_mask .= 0
            result .= 0
            drizzle_warp!(result, drizzle_mask, to_warp, inv_tfm; use_interp=use_interp, supersample=drizzle_supersampling, bayer_pattern)
            return result
        end
        fast_src = let
            if isnothing(fast_src_buf)
                to_fast_mem(to_warp)
            else
                copyto!(slow_src_buf, to_warp)
                copyto!(fast_src_buf, slow_src_buf)
            end
        end
        # result/drizzle_mask's current content is about to be overwritten anyway (they get zeroed right
        # below), so staging it onto fast memory via to_fast_mem would only transfer data we're going to
        # discard. Allocate fresh fast-memory buffers instead, on the same backend fast_src ended up on
        # (unless reusable ones were already provided).
        fast_result = something(fast_result_buf, similar(fast_src, eltype(result), size(result)))
        fast_mask = something(fast_mask_buf, similar(fast_src, eltype(drizzle_mask), size(drizzle_mask)))
        fast_result .= 0
        fast_mask .= 0
        drizzle_warp!(fast_result, fast_mask, fast_src, inv_tfm; use_interp=use_interp, supersample = drizzle_supersampling, bayer_pattern)
        if isnothing(fast_result_buf)
            result .= to_slow_mem(fast_result)
            drizzle_mask .= to_slow_mem(fast_mask)
        else
            # result/drizzle_mask are frequently SubArray views too, so -- symmetrically with the source
            # side above -- always funnel through a concrete Array on both ends of each copyto!.
            copyto!(slow_result_buf, fast_result)
            copyto!(result, slow_result_buf)
            copyto!(slow_mask_buf, fast_mask)
            copyto!(drizzle_mask, slow_mask_buf)
        end
        return result
end

function get_mono(data; use_drizzle, ref_col=(2,1), dim_color = 4)
    if (use_drizzle)
        return @view data[ref_col[1]:2:end, ref_col[2]:2:end]
    elseif (ndims(data)<3)
         return data
    else
         return @view data[:,:,min(size(data,dim_color),ref_col[1])]
    end
end
