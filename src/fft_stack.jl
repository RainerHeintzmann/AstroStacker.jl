"""
    fft_find_transform(src_mono, ref_mono; shift_fun=FindShift.find_shift_iter, kwargs...)

An `Astroalign.find_transform`-compatible transform estimator based on fast, translation-only FFT
cross-correlation ([`FindShift.find_shift_iter`](@ref)/[`FindShift.find_shift_lk`](@ref)) instead of star
detection + RANSAC. This is what [`stack_many_fft`](@ref) passes as `find_transform_fun` to
[`stack_many`](@ref), so the whole existing pipeline (Bayer/drizzle handling, GPU tiering, buffer reuse,
outlier rejection) is reused unchanged -- only the per-frame registration step itself differs.

Since there is no reusable photometry table here (unlike `Astroalign.find_transform`'s `phot_to`), the
returned `params.phot_to` is just `ref_mono` itself, unchanged -- i.e. every frame is registered against
the same fixed reference, which is exactly what's wanted for a pure-translation estimator.

Any other keyword arguments (e.g. Astroalign-specific ones like `box_size`/`ap_radius`/`f`/`min_fwhm`/`N_max`
that `stack_many`'s `kwargs...` may still be carrying over from a copy-pasted call) are accepted and ignored,
so the same call site works whichever `find_transform_fun` is plugged in.
"""
function fft_find_transform(src_mono, ref_mono; shift_fun=FindShift.find_shift_iter, kwargs...)
    Δx = shift_fun(ref_mono, src_mono)
    tfm = AffineMap(SMatrix{2,2}(1.0, 0.0, 0.0, 1.0), SVector(Δx[1], Δx[2]))
    return tfm, (phot_to=ref_mono,)
end

"""
    stack_many_fft(input_stack; shift_fun=FindShift.find_shift_iter, kwargs...)

Stacks a series of frames using fast, translation-only FFT-correlation registration
([`fft_find_transform`](@ref)) in place of `Astroalign.find_transform`'s star-detection + RANSAC rigid/
similarity fit -- otherwise identical to, and sharing the full implementation of, [`stack_many`](@ref):
the same Bayer/drizzle handling, GPU fast/slow-memory tiering (`to_fast_mem`/`to_slow_mem`), buffer reuse,
sigma-clip outlier rejection, and `on_checkpoint` instrumentation all apply unchanged. See
[`stack_many`](@ref)'s docstring for the full list of supported keyword arguments.

This is much cheaper than star-based registration (no source detection/RANSAC), and works on frames with
no discrete point sources at all (e.g. the Moon/planets) -- but it only fits a pure per-frame translation,
no rotation. Good for a quick preview, for frames that only jitter, or to build a more stable reference
frame for [`stack_many_lucky`](@ref) than a single raw input frame.

# Arguments
+ `input_stack`: input frames stacked along the 3rd dimension, as in [`stack_many`](@ref) (a Bayer-pattern
    mosaic by default -- see `use_drizzle`/`bayer_pattern` there).
+ `shift_fun`: the per-frame shift estimator used by [`fft_find_transform`](@ref) -- `FindShift.find_shift_iter`
    (FFT/Optim-based, default) or `FindShift.find_shift_lk` (Lucas-Kanade based alternative).
+ `kwargs...`: forwarded to [`stack_many`](@ref) (e.g. `use_drizzle`, `bayer_pattern`, `drizzle_supersampling`,
    `min_sigma`, `to_fast_mem`, `to_slow_mem`, `on_checkpoint`, `ref_slice`, `verbose`).

Returns the same `(result, all_params)` as [`stack_many`](@ref).
"""
function stack_many_fft(input_stack; shift_fun=FindShift.find_shift_iter, kwargs...)
    find_transform_fun(src_mono, ref_mono; kw...) = fft_find_transform(src_mono, ref_mono; shift_fun, kw...)
    return stack_many(input_stack; find_transform_fun, kwargs...)
end

"""
    apply_shift(input_stack, params)

creates a shifted stack using the shifts stored in `params`, for diagnostic/debug reasons only.
"""
function apply_shift(input_stack, params)
    res = []
    for (slice, p) in zip(eachslice(input_stack, dims=3), params)
        @show size(slice)
        push!(res, FindShift.shift(slice, Tuple(p.tfm.translation)))
    end
    return cat(res..., dims=3)
end
