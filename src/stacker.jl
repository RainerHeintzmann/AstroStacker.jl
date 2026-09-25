"""
    stack_many(input_stack; use_interp = false, use_drizzle=true, drizzle_supersampling = 2.0, min_sigma = 2.0,
                verbose = true, ref_slice = size(input_stack, 3)÷2 + 1, kwargs...)

Stacks many image frames (`input_stack`) stacked along the 3rd dimension into a single result image.
Returns a Tuple of the result image and a list of stacking parameters for each image.

# Parameters
* `input_stack`: input stack to align and sum in the stacking operation. This needs to be a bayer-pattern mosaic. This input stack should have the individual images stacked along dimension 3. 
    Internally first a binned version is calculated and then the transformation parameters are used to transform the original data.
* `use_drizzle`: if `true` the input_stack is interpreted as a bayer pattern and the drizzle algorithm  with the parameter below is applied.
* `use_interp`: if `true` linar interpolation will be used on destination.
* `drizzle_supersampling`: This is the supersampling factor in comparison to one original (red) color sampling. 
    The default of `2` means that the result size will be equal to the original size, but interpolation free.
    It is important to stack enough images such that no holes remain in the stacked image.
* `ref_col`: The index in X and Y as a tuple to use as the reference color channel for alignment only. default=(2,1), which is often the green channel.
* `ref_slice`: an integer indicating the slice to use as a reference image. (default: middle of the stack to minimize field rotation effects).
* `min_sigma`: minimum number of standard deviations a single pixel needs to be away from the mean of that pixel to be excluded. 
               If this number is set to zero, the outlier-exclusion algorithm will not be run.
* `min_fwhm`: minimum FWHM to accept for stars to be considered in the alignment
* `verbose`: prints diagnostic output, if `true`. (default: `true`)
* `bayer_pattern`: a string of size 4 characters, indicating the order of colors. The default ("RGGB") corresponds to this pattern (starting from the top left corner of `input_stack`):

    ```
    R G R G
    G B R B
    R G R G
    G B R B 
    ```
* `box_size`: the box size to use for identifying stars. You should try (15,15), which is not the default.

For other possible (optional) arguments, see the documentation of `align_frames` in the `Astroalign` package.

# GPU usage
`input_stack` is normally kept in "slow" memory (a plain CPU `Array`, or even a lazy/disk-backed array
such as `MultifileArrays.jl`'s `MultifileArray`, which reads a frame from disk only when it's touched) --
`all_results`/`all_masks` are allocated via `similar(input_stack, ...)` and so automatically inherit that
same slow storage. Star detection (used to determine each frame's alignment) is inherently CPU-bound (it
goes through `Astroalign`/`Photometry.jl`) and is transparently run on a small, separately-materialized
CPU copy of the relevant reference/frame data, regardless of `input_stack`'s own storage.

For the actual (per-pixel-expensive) drizzle/forward-warp accumulation, `to_fast_mem`/`to_slow_mem` let
you stage only the single frame currently being processed onto "fast" memory (e.g. a GPU), instead of
requiring the whole multi-frame stack to be resident there at once: pass e.g. `to_fast_mem=cu` (from
`CUDA.jl`) and `to_slow_mem=Array` to run that step on the GPU while `input_stack` stays lazily-loaded/on
the CPU. Both default to `identity` (no staging at all -- fine if `input_stack` is already a GPU array
itself, as in the old/simpler usage pattern, or if you just want to stay on the CPU).

`on_checkpoint`, if given, is called with a short descriptive label (e.g. `"after frame 3"`) at several
points during stacking, so you can trace where memory is actually being allocated -- e.g.
`on_checkpoint = label -> println(label, ": ", CUDA.memory_status())`.
"""
function stack_many(input_stack; use_drizzle=true, use_interp=false, drizzle_supersampling = 2.0, min_sigma = 2.0,
                verbose = true, ref_slice = size(input_stack,3)÷2 + 1, ref_col=(2,1), bayer_pattern = "RGGB",
                to_fast_mem=identity, to_slow_mem=identity, on_checkpoint=nothing, kwargs...)
    # calls on_checkpoint(label) if the caller supplied one -- e.g. `on_checkpoint = label ->
    # println(label, ": ", CUDA.memory_status())` -- to trace where memory (in particular GPU memory) is
    # actually being allocated, without AstroStacker.jl itself depending on any GPU package. A no-op if
    # `on_checkpoint` is left at its default `nothing`.
    checkpoint(label) = isnothing(on_checkpoint) ? nothing : on_checkpoint(label)
    if (!use_drizzle)
        drizzle_supersampling = 1
    end
    dim_color = 4 # see alsot the calculation of the destination size below
    dim_stack = 3
    # Sum over colors (for alignment only)
    # ref_mono = bin_mono(@view input_stack[:, :, ref_slice])[:, :, 1]
    # ref_mono = (use_drizzle) ? (@view input_stack[ref_col[1]:2:end, ref_col[2]:2:end, ref_slice]) : (@view input_stack[:,:,ref_slice])
    # Source-detection (Astroalign/Photometry.jl below) is inherently CPU/scalar-bound, and a plain
    # `input_stack[:,:,ref_slice,:]` materializing getindex on a GPU-array-backed `input_stack` (e.g.
    # wrapped in an OffsetArray by the FITS loader) can hit GPUArrays' "scalar indexing disallowed"
    # guard. Go through a lazy `selectdim` view first, then force ONE explicit bulk copy to a plain CPU
    # array -- this keeps the (large) `input_stack` itself untouched/GPU-resident for the warp/drizzle
    # step below, which now runs as a real (KernelAbstractions-based) GPU kernel; see warp.jl.
    ref_mono = get_mono(Array(selectdim(input_stack, dim_stack, ref_slice)); use_drizzle=use_drizzle, ref_col=ref_col)
    reduced_size = size(ref_mono)[1:2]

    Nimgs = size(input_stack, dim_stack)
    Ncol = 3
    if (!use_drizzle)        
        Ncol = size(input_stack, dim_color)
    end
    dst_size = round.(Int, ((reduced_size .* drizzle_supersampling)..., Nimgs, Ncol))
    all_params = []
    all_results = similar(input_stack, dst_size)
    all_masks = zeros(1,1,size(input_stack, dim_stack),1) # just a dummy to have something to iterate
    n = 1
    # ref_info = nothing

    warp_function = Astroalign.warp;

    # dst_size = round.(Int, ((reduced_size .* drizzle_supersampling)...,3))

    fast_src_buf = nothing
    fast_result_buf = nothing
    fast_mask_buf = nothing
    slow_src_buf = nothing
    slow_result_buf = nothing
    slow_mask_buf = nothing
    if !isnothing(drizzle_supersampling) && (drizzle_supersampling != 1)
        all_masks = similar(input_stack, eltype(all_results), dst_size)
        # Pre-allocate the fast-mem scratch buffers ONCE here (only when genuine staging is actually
        # requested -- to_fast_mem/to_slow_mem not both identity), reused for every frame below instead
        # of (re-)allocating (and, on a GPU backend, cudaMalloc-ing) fresh ones each frame. When both
        # hooks are identity (the default), do_drizzle_warp! operates directly in place with no staging
        # at all -- important for the "pre-cast the whole stack to a CuArray, don't pass to_fast_mem"
        # usage pattern, which must stay fully GPU-resident with zero extra buffers/copies.
        if !(to_fast_mem === identity && to_slow_mem === identity)
            # The slow_*_buf concrete-Array staging buffers exist because GPU array copyto! only has a
            # fast path from/to a concrete Array, not an arbitrary SubArray/lazy view (which
            # src/res_slice/mymask usually are).
            slow_src_buf = zeros(eltype(input_stack), size(selectdim(input_stack, dim_stack, 1)))
            fast_src_buf = to_fast_mem(slow_src_buf)
            slow_result_buf = zeros(eltype(all_results), size(selectdim(all_results, dim_stack, 1)))
            fast_result_buf = to_fast_mem(slow_result_buf)
            slow_mask_buf = zeros(eltype(all_masks), size(selectdim(all_masks, dim_stack, 1)))
            fast_mask_buf = to_fast_mem(slow_mask_buf)
        end
    end
    checkpoint("after allocating all_results/all_masks (+ scratch buffers)")

    myref_mono = ref_mono; # changes to the extracted coordinates after the firt run
    for (src, res_slice, mymask) in zip(eachslice(input_stack; dims = dim_stack), eachslice(all_results, dims = dim_stack), eachslice(all_masks, dims = dim_stack))
        # src_mono = bin_mono(src)[:, :, 1]; # Sum over colors
        # src_mono = (use_drizzle) ? (@view src[ref_col[1]:2:end, ref_col[2]:2:end, 1]) : src
        # forced to a plain CPU array for the same reason as ref_mono above; `src` itself (used for the
        # actual warp/drizzle accumulation below) stays untouched/GPU-resident.
        src_mono = Array(get_mono(src; use_drizzle=use_drizzle, ref_col=ref_col))

        if !isnothing(drizzle_supersampling) && (drizzle_supersampling != 1)
            warp_function(img_from, inv_tfm, myaxes) = do_drizzle_warp!(mymask, drizzle_supersampling, bayer_pattern, use_interp, res_slice, src, inv_tfm, myaxes;
                                                                         to_fast_mem, to_slow_mem, fast_src_buf, fast_result_buf, fast_mask_buf,
                                                                         slow_src_buf, slow_result_buf, slow_mask_buf)
        end

        tfm, params = find_transform(src_mono, myref_mono; kwargs...)
        myref_mono = params.phot_to # to speed up further rounds, sinc find_transform then ignores the photometry

        if (ndims(src) < 3)
            res_slice .= apply_transform(tfm, src_mono, ref_mono; warp_function = warp_function)
        else
            for (src_c, res_c) in zip(eachslice(src; dims = ndims(src)), eachslice(res_slice; dims = ndims(src)))
                res_c .= apply_transform(tfm, src_c, src_c; warp_function = warp_function)
            end
        end

        # if isempty(ref_info)
        #     @warn "ignoring slice $(n)"
        #     n += 1
        #     continue # Ignore this entry
        # end

        params = (params..., tfm = tfm) # store this for later use
        push!(all_params, params)

        if (verbose)
            a = atan(tfm.linear[1, 2], tfm.linear[1, 1]) * 180/pi
            println("stacking: $n, angle: $(round(a; sigdigits = 3)) deg, shift: $(round.(tfm.translation; sigdigits = 4))")
        end

        checkpoint("after frame $n")
        n += 1
    end

    result = nothing # Since it is returned
    checkpoint("before outlier removal / final combination")

    if (min_sigma > 0)
        if (use_drizzle)
            result = remove_outliers(all_results, all_masks; verbose=verbose, stack_dim=dim_stack, min_sigma=min_sigma, on_checkpoint)
        else
            if (size(all_results,dim_color)==1)
                result = remove_outliers(all_results;verbose=verbose, stack_dim=dim_stack, min_sigma=min_sigma, on_checkpoint)
            else
                res_size = ntuple(n -> (n==dim_stack) ? 1 : size(all_results, n), ndims(all_results))
                result = similar(all_results, res_size)
                for (all_c, res_c) in zip(eachslice(all_results; dims = dim_color),eachslice(result; dims = dim_color))
                    res_c .= remove_outliers(all_c; verbose=verbose, stack_dim=dim_stack, min_sigma=min_sigma, on_checkpoint)
                end
            end
        end
    else
        divisor = max.(1, sum(all_masks; dims = dim_stack))
        result = sum(all_results, dims = dim_stack) ./ divisor
    end
    checkpoint("after outlier removal / final combination")

    # Normalize drizzle result
    # result ./= max.(1, all_params[end][:drizzle_mask])
    return result, all_params
end

function remove_outliers(all_results; kwargs...)
    all_masks = .!isnan.(all_results)
    # Eliminate the NaNs
    # boolean-mask indexed assignment (all_results[.!all_masks] .= 0) is scalar iteration under the
    # hood and is disallowed on GPU arrays; ifelse.() is a plain broadcast and works on both CPU and GPU.
    all_results .= ifelse.(all_masks, all_results, zero(eltype(all_results)))
    return remove_outliers(all_results, all_masks; kwargs...)
end

function remove_outliers(all_results, all_masks; verbose = true, stack_dim = 3, min_sigma = 2.0, on_checkpoint=nothing)
        checkpoint(label) = isnothing(on_checkpoint) ? nothing : on_checkpoint(label)
        verbose && println("... summing results")
        divisor = max.(eltype(all_results)(1f-8), sum(all_masks; dims = stack_dim))
        result = sum(all_results; dims = stack_dim) ./ divisor
        checkpoint("remove_outliers: after initial result/divisor")

        verbose && println("... removing outliers")
        stddev_map = weighted_std(all_results, all_masks; dims = stack_dim)
        checkpoint("remove_outliers: after weighted_std")
        n=1
        # remove the stack_dim from the result and stddev_map
        crunched_dims = ntuple(n->(n!=stack_dim) ? (:) : 1, ndims(all_results))
        res_view = @view result[crunched_dims...]
        stddev_view = @view stddev_map[crunched_dims...]
        for (mask, masked_res) in zip(eachslice(all_masks; dims = stack_dim), eachslice(all_results; dims = stack_dim))
            # outliers = (all_masks .!= 0) .&& abs.(masked_res .- result) .> min_sigma .* weighted_std(masked_res, all_masks; dims = stack_dim)
            outliers = (mask .!= 0) .&& abs.(masked_res .- res_view .* mask) .> min_sigma .* stddev_view
            verbose && println("frame $(n) outliers found: $(sum(outliers)), $(round(100*sum(outliers)/length(outliers); sigdigits=3)) %")
            # boolean-mask indexed assignment is scalar iteration under the hood and is disallowed on
            # GPU arrays; ifelse.() is a plain broadcast and works on both CPU and GPU.
            masked_res .= ifelse.(outliers, zero(eltype(masked_res)), masked_res)
            mask .= ifelse.(outliers, zero(eltype(mask)), mask)
            checkpoint("remove_outliers: after outlier-loop frame $n")
            n += 1
        end
        divisor = max.(1, sum(all_masks; dims = stack_dim))
        final = sum(all_results; dims = stack_dim) ./ divisor
        checkpoint("remove_outliers: after final combination")
        return final
end
