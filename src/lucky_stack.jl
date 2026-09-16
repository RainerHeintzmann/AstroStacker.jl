"""
    stack_many_lucky(input_stack; ref_slice=size(input_stack,3)÷2+1, grid_size=(10,10), shift_fun=FindShift.find_shift_iter, quality_power=2.0, verbose=true, align_kwargs...)

Stacks a series of mono 2D frames (`input_stack`, with individual frames stacked along the 3rd
dimension) using patch-wise thin-plate-spline (TPS) non-rigid registration (via
`FindShift.align_images`) followed by local-quality-weighted blending ("lucky imaging"): at each pixel,
frames are weighted according to a smoothly-interpolated local sharpness (Laplacian-variance) map, so
that the sharpest locally-available data dominates the result even when different regions of different
frames are sharp due to spatially-varying (e.g. atmospheric-seeing) blur.

This complements [`stack_many`](@ref), which only performs a single global rigid/similarity registration
per frame; `stack_many_lucky` is intended for extended, textured targets (e.g. the Moon or planets) where
a single per-frame transform cannot capture spatially-varying distortion and sharpness.

# Arguments
+ `input_stack`: input frames stacked along the 3rd dimension (mono, i.e. `(x,y,frame)`).
+ `ref_slice`: index of the frame used as the fixed registration reference.
+ `grid_size`: the TPS/quality patch grid size, passed on to `FindShift.align_images`.
+ `shift_fun`: the per-patch shift-estimation function passed to `FindShift.align_images`/`find_deformations`
    (`FindShift.find_shift_iter` by default, an FFT/Optim-based estimator; `FindShift.find_shift_lk`, a
    Lucas-Kanade based alternative, can be passed in instead).
+ `quality_power`: exponent applied to the (non-negative) quality values before weighting; higher values
    make the blend more aggressively favor the sharpest frame/region ("harder" lucky-imaging selection).
+ `verbose`: print per-frame diagnostic information if `true`.
+ `align_kwargs`: further keyword arguments forwarded to `FindShift.align_images` (e.g. `patch_size`, `tolerance`, `band_pass_freq`).

Returns a `NamedTuple` of `(result, aligned, warps, quality)`, each (except `result`) given in the
original frame order of `input_stack`.
"""
function stack_many_lucky(input_stack::AbstractArray{T,3}; ref_slice=size(input_stack, 3) ÷ 2 + 1, grid_size=(10, 10),
    shift_fun=FindShift.find_shift_iter, quality_power=2.0, verbose=true, align_kwargs...) where {T}
    Nimgs = size(input_stack, 3)
    order = [ref_slice, (n for n in 1:Nimgs if n != ref_slice)...]
    frames = [Float64.(input_stack[:, :, n]) for n in order]

    aligned, warps = FindShift.align_images(frames; grid_size=grid_size, shift_fun=shift_fun, align_kwargs...)

    nodes = FindShift.get_default_markers(aligned[1], grid_size)
    img_size = size(aligned[1])

    qualities = Vector{Matrix{Float64}}(undef, length(aligned))
    weights = Vector{Matrix{Float64}}(undef, length(aligned))
    for (n, frame) in enumerate(aligned)
        q = patch_quality_grid(frame, grid_size)
        qualities[n] = q
        weights[n] = quality_weight_map(q, nodes, img_size; power=quality_power)
        verbose && println("lucky-stacking: frame $(n), mean patch quality: $(round(mean(q); sigdigits=3))")
    end

    acc = zeros(Float64, img_size)
    wsum = zeros(Float64, img_size)
    for (frame, w) in zip(aligned, weights)
        acc .+= w .* frame
        wsum .+= w
    end
    result = acc ./ max.(wsum, eps(Float64))

    # undo the reordering, so the per-frame outputs line up with the original frame order of input_stack
    inv_order = sortperm(order)
    return (result=result, aligned=aligned[inv_order], warps=warps[inv_order], quality=qualities[inv_order])
end
