"""
    stack_many_fft(input_stack; ref_slice=nothing, min_sigma=2.0, verbose=true, align_kwargs...)

Stacks a series of mono 2D frames (`input_stack`, with individual frames stacked along the 3rd
dimension) using fast, translation-only FFT-correlation registration (`FindShift.align_stack`),
rather than the full patch-wise thin-plate-spline non-rigid registration used by
[`stack_many_lucky`](@ref).

This is a much cheaper option than [`stack_many_lucky`](@ref), suitable when frames only jitter
(no rotation/scale, no spatially-varying local distortion) -- e.g. for a quick preview, or to build a
more stable reference frame for [`stack_many_lucky`](@ref) than a single raw input frame.

# Arguments
+ `input_stack`: input frames stacked along the 3rd dimension (mono, i.e. `(x,y,frame)`).
+ `ref_slice`: index of the reference frame to align to (default: the middle frame, see `FindShift.align_stack`'s `refno`).
+ `min_sigma`: minimum number of standard deviations a pixel needs to be away from the mean of that
    pixel to be excluded as an outlier (sigma-clipping, as in [`stack_many`](@ref)). Set to `0` to
    disable outlier rejection (plain mean).
+ `verbose`: print per-frame shift information if `true`.
+ `align_kwargs`: further keyword arguments forwarded to `FindShift.align_stack` (e.g. `damp`, `max_freq`, `method`).

Returns a `NamedTuple` of `(result, aligned, shifts)`.
"""
function stack_many_fft(input_stack::AbstractArray{T,3}; ref_slice=nothing, min_sigma=2.0, verbose=true, align_kwargs...) where {T}
    aligned, shifts = FindShift.align_stack(input_stack; refno=ref_slice, align_kwargs...)

    if verbose
        for (n, s) in enumerate(shifts)
            println("stacking (fft): frame $n, shift: $(round.(s[1:2]; sigdigits=4))")
        end
    end

    result = let
        if min_sigma > 0
            remove_outliers(aligned; verbose=verbose, stack_dim=3, min_sigma=min_sigma)
        else
            sum(aligned; dims=3) ./ size(aligned, 3)
        end
    end

    return (result=dropdims(result, dims=3), aligned=aligned, shifts=shifts)
end
