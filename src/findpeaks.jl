"""
    get_sources(img; box_size = nothing, nsigma = 1, N_max = 10)

Extract candidate sources in `img` according to [`Photometry.Detection.extract_sources`](@extref). By default, `img` is first sigma clipped and then background subtracted before the candidate sources are extracted. `box_size` is passed to [`Photometry.Background.estimate_background`](@extref), and `nsigma` is passed to [`Photometry.Detection.extract_sources`](@extref). See the [Photometry.jl](@extref) documentation for more.

TODO: Pass more options to clipping, background estimating, and extraction methods in [Photometry.jl](@extref).
"""
function get_sources(img; box_size = _compute_box_size(img), nsigma = 1, N_max = 10)
    # Background subtract `img`
    clipped = sigma_clip(img, 1, fill = NaN)
    bkg, bkg_rms = estimate_background(clipped, box_size)
    subt = img .- bkg[axes(img)...]

    return (
        # Sort detected sources from brightest to darkest
        first(extract_sources(PeakMesh(; box_size, nsigma), subt, bkg, true), N_max),
        # And also return the inputs, handy for debugging and data viz
        subt,
        bkg,
    )
end

"""
    com_psf(img_ap)

determins peak parameters via a fast, non-iterative center-of-mass approach.

returned is a named Tuple:
(; psf_params, psf_model, psf_data)
with the `psf_params` also being a named tuple 
(; x, y, fwhm)
"""
function com_psf(img_ap; rel_thresh=0.1f0)
    psf_data = collect(Float32, img_ap)
    psf_data ./= maximum(psf_data)

    x_coords = axes(psf_data, 1)
    y_coords = axes(psf_data, 2)'
    # use the average of the edge of the ROI as background estimate
    bg = (sum(psf_data[1,:])+sum(psf_data[end,:])+sum(psf_data[:,1])+sum(psf_data[:,end])) / 
        (2*size(psf_data,1)+2*size(psf_data,2))

    # clamp at zero before com calculation
    psf_data = max.(0, psf_data .- bg .- rel_thresh)
    sum_psf = sum(psf_data)
    x = sum(x_coords .* psf_data) / sum_psf
    y = sum(y_coords .* psf_data) / sum_psf

    var_x = sum(abs2.(x_coords.-x) .* psf_data) / sum_psf
    var_y = sum(abs2.(y_coords.-y) .* psf_data) / sum_psf
    fac = 2*sqrt(2*log(2))
    fwhm = (sqrt(var_x)*fac, sqrt(var_y)*fac)

    psf_params = (; x, y, fwhm)
    psf_model = "com"
    return (; psf_params, psf_model, psf_data)
end
