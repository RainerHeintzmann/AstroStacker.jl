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
# using CUDA

# -----------------------------------------------------------------------------
# DATA DOWNLOAD SECTION
# Run this block (set true) to download example astronomical image data
# The dataset (~500 MB) contains exposures from a Dwarf III telescope
# -----------------------------------------------------------------------------

function main()

        # -----------------------------------------------------------------------------
        # FILE PATH CONFIGURATION
        # Define paths to image data, dark frames, and flat fields
        # Toggle between example data and custom backup paths
        # -----------------------------------------------------------------------------
        folder = raw"C:\NoBackup\dwarf\Astronomy\DWARF_RAW_TELE_Moon_EXP_0.004_GAIN_10_2026-09-26-19-28-35-816\\"
        
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
        files = raw"Moon_0.004s10_VIS_20260926-*_33C.fits"

        # -----------------------------------------------------------------------------
        # LOAD CALIBRATION FRAMES AND IMAGE DATA
        # Load master dark and flat field, then apply calibration to science images
        # -----------------------------------------------------------------------------
        # data = load_series(load, files);

        # Load master dark frame (averaged dark current subtraction reference)
        dark = load(joinpath(darkfolder, file_dark15));

        # Load flat field (dust motes and vignetting correction reference)
        flat = load(joinpath(flatfolder, file_flat));

        # Workaround for Windows path handling issues with FITS loading
        # Save current directory, change to data folder, load images, restore
        curpath = pwd()
        cd(folder)  # Required due to Windows path handling bug
        data = load_series(load, files)
        cd(curpath)

        # Apply dark subtraction and flat field division to all images
        # Note: Conversion to Float32 is required - Float64 causes FITS loading issues
        # @time data = correct_dark_flat(data, dark, flat);
        # @time data = correct_dark_flat(data, dark);

        grid_size = (10,10)
        quality_power = 2.0
        data = collect(data)
        single_channel = data[2:2:end, 1:2:end, :];

        bayer_pattern = "RGGB"
        @time fft_stacked_bayer, all_params_fft_bayer = stack_many_fft(Float32.(data); bayer_pattern=bayer_pattern, drizzle_supersampling=2.0);

        fft_stacked_mono, all_params_fft_mono = stack_many_fft(single_channel; use_drizzle=false);
        shifted = AstroStacker.apply_shift(Float32.(single_channel), all_params_fft_mono)
        @vt single_channel shifted
        
        plot()
        @time fft_stacked_mono, all_params_fft_mono = stack_many_fft(Float32.(single_channel); drizzle_supersampling=2.0);
        # result_fft, aligned_fft, shifts_fft = stack_many_fft(single_channel)

        result, aligned, warps, quality = stack_many_lucky(single_channel; grid_size=grid_size, quality_power=quality_power, verbose=true)


        # dist_limit = 2, 
        # using ProfileCanvas
        # ProfileCanvas.@profview stacked_d, all_params_d = stack_many(data; use_interp=use_interp, use_drizzle=true, f=f, N_max=N_max,
        #         box_size, ap_radius, min_sigma = 2.5, nsigma = 1, min_fwhm = min_fwhm, drizzle_supersampling = 2.0)

        # defining a new axis na
        # newaxis = [CartesianIndex()]
        # @vt stacked_d # [:,:,newaxis,:]
        # @vt sum(data, dims=4)

        # all_med_fwhms_x, all_med_fwhms_y, 
end
