"""
    do_drizzle_warp!(drizzle_mask, to_warp, inv_tfm, myaxes, )

an alternative to the `warp` function, to be provided to the alignment function instead of warp.

# Parameters:

- `drizzle_mask`: a mask, collecting which pixel where assigned.
- `result`: the assigned pixels.
- `use_interp`: whether to use interpolation (true) or not.
- `bayer_pattern`: a string of size 4 characters, indicating the order of colors. The default ("RGGB") corresponds to this pattern (starting from the top left corner of `input_stack`):
- `myaxes`: is ignored.
"""
function do_drizzle_warp!(drizzle_mask, drizzle_supersampling, bayer_pattern, use_interp, result, to_warp, inv_tfm, myaxes)
        isnothing(to_warp) && error("For drizzle you need to provide a drizzle_supersample! and a to_warp input, the Bayer-pattern mosaic input")
        # drizzle_mask = similar(to_warp, eltype(to_warp), dst_size)
        drizzle_mask .= 0
        # result = similar(to_warp, dst_size)
        result .= 0
        warped = drizzle_warp!(result, drizzle_mask, to_warp, inv_tfm; use_interp=use_interp, supersample = drizzle_supersampling, bayer_pattern)
        return warped
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

For other possible (optional) arguments, see the documentation of `align_frame`.
"""
function stack_many(input_stack; use_drizzle=true, use_interp=false, drizzle_supersampling = 2.0, min_sigma = 2.0,
                verbose = true, ref_slice = size(input_stack,3)÷2 + 1, ref_col=(2,1), bayer_pattern = "RGGB", kwargs...)
    if (!use_drizzle) 
        drizzle_supersampling = 1
    end 
    # Sum over colors (for alignment only)
    # ref_mono = bin_mono(@view input_stack[:, :, ref_slice])[:, :, 1]
    # ref_mono = (use_drizzle) ? (@view input_stack[ref_col[1]:2:end, ref_col[2]:2:end, ref_slice]) : (@view input_stack[:,:,ref_slice])
    ref_mono = get_mono(input_stack[:,:,ref_slice,:]; use_drizzle=use_drizzle, ref_col=ref_col)
    reduced_size = size(ref_mono)[1:2]

    dim_color = 4 # see alsot the calculation of the destination size below
    dim_stack = 3
    Nimgs = size(input_stack, dim_stack)
    Ncol = 3
    if (!use_drizzle)        
        Ncol = size(input_stack, dim_color)
    end
    @show dst_size = round.(Int, ((reduced_size .* drizzle_supersampling)..., Nimgs, Ncol))
    all_params = []
    all_results = similar(input_stack, dst_size)
    all_masks = zeros(1,1,size(input_stack, dim_stack),1) # just a dummy to have something to iterate
    n = 1
    # ref_info = nothing

    warp_function = Astroalign.warp;

    # dst_size = round.(Int, ((reduced_size .* drizzle_supersampling)...,3))

    if !isnothing(drizzle_supersampling) && (drizzle_supersampling != 1)
        all_masks = similar(input_stack, eltype(all_results), dst_size)
    end

    for (src, res_slice, mymask) in zip(eachslice(input_stack; dims = dim_stack), eachslice(all_results, dims = dim_stack), eachslice(all_masks, dims = dim_stack))
        # src_mono = bin_mono(src)[:, :, 1]; # Sum over colors
        # src_mono = (use_drizzle) ? (@view src[ref_col[1]:2:end, ref_col[2]:2:end, 1]) : src
        src_mono = get_mono(src; use_drizzle=use_drizzle, ref_col=ref_col)

        if !isnothing(drizzle_supersampling) && (drizzle_supersampling != 1)
            warp_function(img_from, inv_tfm, myaxes) = do_drizzle_warp!(mymask, drizzle_supersampling, bayer_pattern, use_interp, res_slice, src, inv_tfm, myaxes)
        end

        tfm, params = find_transform(src_mono, ref_mono; kwargs...)

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

        n += 1
    end

    result = nothing # Since it is returned

    if (min_sigma > 0)
        if (use_drizzle)
            @show size(all_results)
            @show size(all_masks)
            result = remove_outliers(all_results, all_masks; verbose=verbose, stack_dim=dim_stack, min_sigma=min_sigma)
        else
            if (size(all_results,dim_color)==1)
                result = remove_outliers(all_results;verbose=verbose, stack_dim=dim_stack, min_sigma=min_sigma)
            else
                res_size = ntuple(n -> (n==dim_stack) ? 1 : size(all_results, n), ndims(all_results))
                result = similar(all_results, res_size)
                for (all_c, res_c) in zip(eachslice(all_results; dims = dim_color),eachslice(result; dims = dim_color))
                    res_c .= remove_outliers(all_c; verbose=verbose, stack_dim=dim_stack, min_sigma=min_sigma)
                end
            end
        end
    else
        divisor = max.(1, sum(all_masks; dims = dim_stack))
        result = sum(all_results, dims = dim_stack) ./ divisor
    end

    # Normalize drizzle result
    # result ./= max.(1, all_params[end][:drizzle_mask])
    return result, all_params
end

function remove_outliers(all_results; kwargs...)
    all_masks = .!isnan.(all_results)
    # Eliminate the NaNs
    all_results[.!all_masks] .= 0
    return remove_outliers(all_results, all_masks; kwargs...)
end

function remove_outliers(all_results, all_masks; verbose = true, stack_dim = 3, min_sigma = 2.0)
        verbose && println("... summing results")
        divisor = max.(eltype(all_results)(1f-8), sum(all_masks; dims = stack_dim))
        result = sum(all_results; dims = stack_dim) ./ divisor

        verbose && println("... removing outliers")
        stddev_map = weighted_std(all_results, all_masks; dims = stack_dim)
        n=1
        # remove the stack_dim from the result and stddev_map 
        crunched_dims = ntuple(n->(n!=stack_dim) ? (:) : 1, 4)
        res_view = @view result[crunched_dims...]
        stddev_view = @view stddev_map[crunched_dims...]
        for (mask, masked_res) in zip(eachslice(all_masks; dims = stack_dim), eachslice(all_results; dims = stack_dim))
            # outliers = (all_masks .!= 0) .&& abs.(masked_res .- result) .> min_sigma .* weighted_std(masked_res, all_masks; dims = stack_dim)
            outliers = (mask .!= 0) .&& abs.(masked_res .- res_view .* mask) .> min_sigma .* stddev_view
            verbose && println("frame $(n) outliers found: $(sum(outliers)), $(round(100*sum(outliers)/length(outliers); sigdigits=3)) %")
            masked_res[outliers] .= 0
            mask[outliers] .= 0
            n += 1
        end
        divisor = max.(1, sum(all_masks; dims = stack_dim))
        return sum(all_results; dims = stack_dim) ./ divisor
end
