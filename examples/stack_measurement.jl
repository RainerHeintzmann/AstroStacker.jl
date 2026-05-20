# =============================================================================
# stack_measurement.jl - AstroStacker Example Script
# =============================================================================
# This script demonstrates how to use AstroStacker to stack astronomical images
# from multiple exposures. It shows three different stacking approaches:
#   1. Color camera with Bayer pattern (RGGB)
#   2. Monochromatic single channel
#   3. Pre-debayered color images
# The script also generates diagnostic plots showing alignment parameters
# and displays the final stacked images.
# =============================================================================

# Load required packages
using AstroStacker          # Core stacking functionality
# using Pluto                # Interactive notebooks (commented out)
using Plots                 # Plotting and visualization
using Images                # Image display utilities
# using View5D               # Interactive 5D viewer (commented out)
using Statistics: mean, median  # Statistical functions
using NDTools: select_region    # N-dimensional array tools
using AstroImages           # FITS file format support
using FileIO                # Generic file loading interface
using MultifileArrays: load_series  # Load image series from files

# -----------------------------------------------------------------------------
# DATA DOWNLOAD SECTION
# Run this block (set true) to download example astronomical image data
# The dataset (~500 MB) contains exposures from a Dwarf III telescope
# -----------------------------------------------------------------------------
if (false)
    using Downloads
    
    # Create temporary directory for download
    mydir = mktempdir()
    
    # Example data URL from University of Jena cloud storage
    zip_url = "https://cloud.uni-jena.de/s/9dqYsosX2DxTT5G/download"
    
    # Download ZIP archive containing FITS image files
    Downloads.download(zip_url, joinpath(mydir, "data.zip"))
    
    # Extract ZIP file using tar (works on Windows with modern tar utility)
    # This creates the examples/example_data directory structure
    run(`tar -xf $(joinpath(mydir, "data.zip"))`)
end


function main()

        # -----------------------------------------------------------------------------
        # FILE PATH CONFIGURATION
        # Define paths to image data, dark frames, and flat fields
        # Toggle between example data and custom backup paths
        # -----------------------------------------------------------------------------
        if (true)
                # Use example_data folder with calibration frames in the same directory
                folder = "example_data\\"
                darkfolder = folder      # Dark frames located in same folder
                flatfolder = folder      # Flat frames located in same folder
                
                # Master dark frame for 15 second exposures
                file_dark15 = raw"dark_exp_15.000000_gain_60_bin_1_12C_stack_9.fits"
                
                # Flat field calibration image
                file_flat = raw"flat_gain_2_bin_1_ir_0.fits"
                
                # Pattern matching for M35 star cluster exposures (15s, date 2026-02-15)
                files = raw"M 35_15s60_Astro_20260215-*_16C.fits"
        else # other local file paths.
                # Alternative: Custom paths from backup location
                #folder = raw"C:\NoBackup\dwarf\Astronomy\DWARF_RAW_TELE_M 35_EXP_15_GAIN_60_2026-02-15-20-07-10-140\\"
                folder = raw"C:\NoBackup\dwarf\Astronomy\DWARF_RAW_TELE_M 101_EXP_60_GAIN_60_2026-03-22-23-55-44-199\\"
                
                # Separate calibration folders organized by type
                darkfolder = raw"C:\NoBackup\dwarf\Astronomy\CALI_FRAME\dark\cam_0\\"
                flatfolder = raw"C:\NoBackup\dwarf\Astronomy\CALI_FRAME\flat\cam_0\\"
                
                # Master dark for 60 second exposures at 12C
                file_dark15 = raw"dark_exp_60.000000_gain_60_bin_1_14C_stack_10.fits"
                # file_dark15 = raw"dark_exp_15.000000_gain_60_bin_1_12C_stack_9.fits"
                
                # Flat field calibration
                file_flat = raw"flat_gain_2_bin_1_ir_0.fits"
                
                # Pattern matching for M101 galaxy exposures (60s, date 2026-03-23)
                # files = raw"M 35_15s60_Astro_20260215-*_16C.fits"
                files = raw"M 101_60s60_Astro_20260323-*_14C.fits"
        end


        # -----------------------------------------------------------------------------
        # LOAD CALIBRATION FRAMES AND IMAGE DATA
        # Load master dark and flat field, then apply calibration to science images
        # -----------------------------------------------------------------------------
        # data = load_series(load, files);

        # Load master dark frame (averaged dark current subtraction reference)
        dark = load(joinpath(darkfolder, file_dark15))

        # Load flat field (dust motes and vignetting correction reference)
        flat = load(joinpath(flatfolder, file_flat))

        # Workaround for Windows path handling issues with FITS loading
        # Save current directory, change to data folder, load images, restore
        curpath = pwd()
        cd(folder)  # Required due to Windows path handling bug
        data = load_series(load, files)
        cd(curpath)

        # Apply dark subtraction and flat field division to all images
        # Note: Conversion to Float32 is required - Float64 causes FITS loading issues
        @time data = correct_dark_flat(data, dark, flat)

        # -----------------------------------------------------------------------------
        # STACKING PARAMETERS CONFIGURATION
        # Define parameters for source detection, alignment, and stacking
        # -----------------------------------------------------------------------------
        box_size = (15, 15)           # Detection box size in pixels for peak finding
        ap_radius = 0.6 * first(box_size)  # Aperture radius for photometry measurements
        min_fwhm = 1.5                # Minimum full width at half maximum for star detection
        N_max = 30                    # Maximum number of reference sources to use
        use_interp = true             # Use interpolation for sub-pixel alignment

        # Source detection function selection
        # f = AstroStacker.Astroalign.PSF()  # Alternative: PSF-fitting approach
        f = com_psf                   # Use common centroid method for alignment

        # ----- stack color camera images which follow a bayer pattern "RGGB" -------
        bayer_pattern = "RGGB"
        @time stacked_d, all_params_d = stack_many(data; use_interp=use_interp, use_drizzle=true, f=f, N_max=N_max,
                box_size=box_size, ap_radius=ap_radius, min_sigma = 2.5, nsigma = 1, min_fwhm = min_fwhm, bayer_pattern = bayer_pattern, drizzle_supersampling = 2.0);

        # dist_limit = 2, 
        # using ProfileCanvas
        # ProfileCanvas.@profview stacked_d, all_params_d = stack_many(data; use_interp=use_interp, use_drizzle=true, f=f, N_max=N_max,
        #         box_size, ap_radius, min_sigma = 2.5, nsigma = 1, min_fwhm = min_fwhm, drizzle_supersampling = 2.0)

        # defining a new axis na
        # newaxis = [CartesianIndex()]
        # @vt stacked_d # [:,:,newaxis,:]
        # @vt sum(data, dims=4)

        # all_med_fwhms_x, all_med_fwhms_y, 
        all_stars_used, all_shift_x, all_shift_y, all_rotation, all_med_fwhms_x, all_med_fwhms_y = collect_info(all_params_d);
        plot(all_shift_x, title="Shifts", xlabel="frame #", ylabel="shift / pixel", label="X");plot!(all_shift_y, label="Y")
        plot(all_med_fwhms_x, title="FWHMs", xlabel="frame #", ylabel="shift / pixel", label="X");plot!(all_med_fwhms_y, label="Y")
        plot(all_rotation .* (180/pi), title="rotation", xlabel="frame #", ylabel="angle / deg", label="X")

        prepare_for_viewer(v) = sqrt.(max.(0, v .- median(v, dims=(1,2)))[:,:,1,:]) #reshape( (size(v)[1:2]...,1,size(v,3))
        prepare_for_display(v, m=1.0) = m .*colorview(RGB, permutedims(prepare_for_viewer(v), (3,2,1)))
        mono_for_display(v, m=1.0) = m .*Gray.(permutedims(prepare_for_viewer(v)[:,:,1], (2,1)))
        plot(prepare_for_display(stacked_d, 0.08)) 
        # @vt prepare_for_viewer(stacked_d)

        # --------- Now let's stack a monochromatic single channel only -------------
        single_channel = data[2:2:end, 1:2:end, :];
        @time stacked_s, all_params_s = stack_many(single_channel; use_drizzle=false, f=f, N_max=N_max,
                box_size, ap_radius, min_sigma = 1.5, nsigma = 1, min_fwhm = min_fwhm);

        all_stars_used, all_shift_x, all_shift_y, all_rotation, all_med_fwhms_x, all_med_fwhms_y = collect_info(all_params_s);
        plot(all_shift_x, title="Shifts", xlabel="frame #", ylabel="shift / pixel", label="X");plot!(all_shift_y, label="Y")
        plot(all_med_fwhms_x, title="FWHMs", xlabel="frame #", ylabel="shift / pixel", label="X");plot!(all_med_fwhms_y, label="Y")

        # prepare_for_display(v, m=10.0) = m .*colorview(RGB,permutedims(prepare_for_viewer(v), (4, 2, 1, 3))[:,:,:,1])
        # display the result as an RGB image
        plot(mono_for_display(stacked_s, 0.08)) 


        # -------- and now a color image, which is already de-bayered ------------
        # bin first and process the binned color data
        all_color = bin_rgb(data); # color is in the 4th dimension
        box_size = (9, 9)
        ap_radius = 0.6 * first(box_size);
        stacked_c, all_params_c = stack_many(all_color; use_drizzle=false, 
                f=f, N_max=N_max, box_size=box_size, ap_radius=ap_radius, min_sigma = 2.5, nsigma = 1, min_fwhm = min_fwhm);

        # @vt prepare_for_viewer(stacked_c)
        plot(prepare_for_display(stacked_c, 0.08)) 

        # bin and sum colors first and process the binned monochrome data
        all_binned_m = bin_mono(data)[:,:,:,1];
        box_size = (9, 9)
        ap_radius = 0.6 * first(box_size);
        @time stacked_m, all_params_m = stack_many(all_binned_m; use_drizzle=false, f=f, N_max=N_max,
                box_size, ap_radius, min_sigma = 1.5, nsigma = 1, min_fwhm = min_fwhm);

        plot(mono_for_display(stacked_m, 0.08)) 
        # heatmap(sqrt.(clamp.(stacked_m[:,:,1,1], 200, 250)))
        @vt prepare_for_viewer(stacked_m)
end

function better_speed()
        ##### speed improvements ?
        using ImageFiltering
        data = rand(2048, 2048)
        box_size = (7, 7)
        @time data_max = mapwindow(maximum, data, box_size, border = Fill(zero(eltype(data))))

        using Photometry
        using TypedTables

        using BenchmarkTools
        pm = PeakMesh(box_size, 3.0)
        @btime s =  extract_sources($pm, $data); #  384 ms
        # @btime sm =  extract_sources3($pm, $data); # 58 ms
end
