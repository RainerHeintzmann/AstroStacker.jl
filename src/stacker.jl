# """
#     stack_many(img_stack; verbose = true, ref_slice = size(img_stack,3)÷2+1, min_sigma=2.0, kwargs...)

# Stacks many image frames (`img_stack`) stacked along the 3rd dimension into a single result image.

# # Parameters

# * `img_stack`: input stack to align and sum in the stacking operation. This input stack should have the individual images stacked along dimension 3. It can have multiple colors, stacked along the 4th dimension.
# * `verbose`: prints diagnostic output, if `true`. (default: `true`)
# * `ref_slice`: an integer indicating the slice to use as a reference image. (default: middle of the stack to minimize field rotation effects).
# * `min_sigma`: minimum number of standard deviations a single pixel needs to be away from the mean of that pixel to be excluded. 
#                If this number is set to zero, the outlier-exclusion algorithm will not be run.

# For other (optional) arguments, see the documentation of `align_frame`.

# # Example

# ```julia
# using IndexFunArrays # For gaussian blob generation

# # Create Initial star coordinates
# sz = (100,100); mid_pos = [sz...] ./ 2
# N = 60
# star_pos = sz .* rand(2,N)
# star_amp = rand(N)
# star_shape = (0.3 .+ rand(2, N))

# # Create star image
# y1 = gaussian(sz; offset = star_pos, weight = star_amp, sigma = star_shape)
# rot_mat(alpha) = [cos(alpha) sin(alpha); -sin(alpha) cos(alpha)]

# # Rotate the stars by random angles
# frames = []
# coords = []
# N = 10

# for alpha in 33*rand(N)
#     shift_vec = 10.0 .* (rand(2).-0.5)
#     new_pos = rot_mat(alpha * pi/180) * (star_pos .- mid_pos) .+ mid_pos .+ shift_vec
#     frame = gaussian(sz; offset = new_pos, weight = star_amp, sigma = star_shape);
#     push!(frames, frame); push!(coords, new_pos)
# end

# frames = cat(frames...; dims=3)
# ```
# """
# function stack_many(img_stack; verbose = true, ref_slice = size(img_stack,3)÷2 + 1, min_sigma = 2.0, kwargs...)
#     num_cols = size(img_stack, 4)
#     Nimgs = size(img_stack, 3)
#     dst_size = (size(img_stack)[1:2]..., num_cols, Nimgs)
#     all_params = []
#     all_results = similar(img_stack, dst_size)
#     all_results .= 0

#     # Keep the alignment images in mono
#     # Do NOT use a @view in the line below, as this leads to errors!
#     ref_img = sum(img_stack[:, :, ref_slice:ref_slice, :]; dims = 4)[:, :, 1, 1]
#     # ref_img = @view img_stack[:,:,ref_slice,:]
#     ref_info = nothing
#     n=1

#     for (src_color, res_slice) in zip(eachslice(img_stack; dims = 3), eachslice(all_results; dims = 4))
#         # sum over colors. Note that eachslice removes dimension 3
#         src_mono = sum(src_color; dims = 3)[:, :, 1, 1]
#         myres, _, ref_info, params = align_frame(ref_img, src_mono; to_warp = src_color, kwargs...)

#         if isempty(ref_info)
#             @warn "ignoring slice $(n)"
#             n += 1
#             continue # ignore this entry
#         end

#         push!(all_params, params)
#         res_slice .= myres

#         if (verbose)
#             tfm = params[:tfm]
#             a = atan(tfm.linear[1,2], tfm.linear[1,1]) * 180/pi
#             println("stacking frame $n, angle: $a deg, shift: $(tfm.translation)")
#         end

#         n += 1
#     end

#     result = nothing
#     stack_dim = 4

#     if (min_sigma > 0)
#         # @show size(all_results)
#         # @show min_sigma
#         result = remove_outliers(all_results; verbose, stack_dim, min_sigma)
#     else
#         all_masks = .!isnan.(all_results)
#         # Eliminate the NaNs
#         all_results[.!all_masks] .= 0
#         result = sum(all_results; dims=stack_dim) ./ max.(1, sum(all_masks; dims = stack_dim))
#     end

#     return result, all_params
# end

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

"""
    stack_many(input_stack; use_interp = false, use_drizzle=true, drizzle_supersampling = 2.0, min_sigma = 2.0,
                verbose = true, ref_slice = size(input_stack, 3)÷2 + 1, kwargs...)

Stacks many image frames (`input_stack`) stacked along the 3rd dimension into a single result image.

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
    ref_mono = (use_drizzle) ? (@view input_stack[ref_col[1]:2:end, ref_col[2]:2:end, ref_slice]) : (@view input_stack[:,:,ref_slice])

    reduced_size = size(ref_mono)[1:2]

    Nimgs = size(input_stack, 3)
    NZ = 3
    if (!use_drizzle)
        NZ = 1
    end
    dst_size = round.(Int, ((reduced_size .* drizzle_supersampling)..., NZ, Nimgs))
    all_params = []
    all_results = similar(input_stack, dst_size)
    all_masks = nothing
    n = 1
    # ref_info = nothing

    final_warp_function = Astroalign.warp;

    # dst_size = round.(Int, ((reduced_size .* drizzle_supersampling)...,3))

    if !isnothing(drizzle_supersampling) && (drizzle_supersampling != 1)
        all_masks = similar(input_stack, eltype(all_results), dst_size)
    end

    for (src, res_slice, mymask) in zip(eachslice(input_stack; dims = 3), eachslice(all_results, dims = 4), eachslice(all_masks, dims = 4))
        # src_mono = bin_mono(src)[:, :, 1]; # Sum over colors
        src_mono = (use_drizzle) ? (@view src[ref_col[1]:2:end, ref_col[2]:2:end, 1]) : src

        if !isnothing(drizzle_supersampling) && (drizzle_supersampling != 1)
            final_warp_function(img_from, inv_tfm, myaxes) = do_drizzle_warp!(mymask, drizzle_supersampling, bayer_pattern, use_interp, res_slice, src, inv_tfm, myaxes)
        end

        # ref_info, 
        myres, params = align_frame(src_mono, ref_mono;
            # use_interp=use_interp, drizzle_supersampling = drizzle_supersampling,
            # to_warp = src,
            # ref_info = ref_info,
            # verbose = verbose,
            final_warp_function = final_warp_function,
            kwargs...
        )

        # if isempty(ref_info)
        #     @warn "ignoring slice $(n)"
        #     n += 1
        #     continue # Ignore this entry
        # end

        push!(all_params, params)
        res_slice .= myres

        if (verbose)
            tfm = params[:tfm]
            a = atan(tfm.linear[1, 2], tfm.linear[1, 1]) * 180/pi
            println("stacking: $n, angle: $(round(a; sigdigits = 3)) deg, shift: $(round.(tfm.translation; sigdigits = 4))")
        end

        n += 1
    end

    result = nothing # Since it is returned
    stack_dim = 4

    if (min_sigma > 0)
        if (use_drizzle)
            result = remove_outliers(all_results, all_masks; verbose, stack_dim, min_sigma)
        else
            result = remove_outliers(all_results; verbose, stack_dim, min_sigma)
        end
    else
        divisor = max.(1, sum(all_masks; dims = stack_dim))
        result = sum(all_results, dims = stack_dim) ./ divisor
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

function remove_outliers(all_results, all_masks; verbose = true, stack_dim = 4, min_sigma = 2.0)
        verbose && println("... summing results")
        divisor = max.(eltype(all_results)(1f-8), sum(all_masks; dims = stack_dim))
        result = sum(all_results; dims = stack_dim) ./ divisor

        verbose && println("... removing outliers")
        stddev_map = weighted_std(all_results, all_masks; dims = stack_dim)
        n=1
        for (mask, masked_res) in zip(eachslice(all_masks; dims = stack_dim), eachslice(all_results; dims = stack_dim))
            # outliers = (all_masks .!= 0) .&& abs.(masked_res .- result) .> min_sigma .* weighted_std(masked_res, all_masks; dims = stack_dim)
            outliers = (mask .!= 0) .&& abs.(masked_res .- result.*mask) .> min_sigma .* stddev_map
            verbose && println("frame $(n) outliers found: $(sum(outliers)), $(round(100*sum(outliers)/length(outliers); sigdigits=3)) %")
            masked_res[outliers] .= 0
            mask[outliers] .= 0
            n += 1
        end
        divisor = max.(1, sum(all_masks; dims = stack_dim))
        return sum(all_results; dims = stack_dim) ./ divisor
end
