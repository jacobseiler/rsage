#!/usr/bin/env python
"""
This file is the driver script that creates the plots for the Seiler et. al
(2018) paper.  The associated files ``GalaxyData.py`` and ``GalaxyPlots.py``
control the reading in of galaxy data and the creation of the plots.  Files
``ReionData.py`` and ``ReionPlots.py`` perform similar roles for the
reionization plots.

We also include a number of extra plots that can be turned on/off as desired.
Please refer to the ``README`` in this directory (``output``) for more
information on how to use this for your own data.

Author: Jacob Seiler
Version: 0.4
"""

from __future__ import print_function

import GalaxyData as galdata
import ReionData as reiondata
import MiscData as misc
import AllVars as av
import PlotScripts as ps

import numpy as np 
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

output_format = "png"


if __name__ == "__main__":

    av.Set_Constants()

    # Set label sizes, colours etc.
    ps.rsage_paper_plot_params()

    # Directory name where we want the plots to be saved.
    output_directory = "./2018_planck"

    # Format all plots are saved as.
    output_format = "png"

    # Plotting is driven entirely through specifying the .ini files. 
    # For this reason, the directories specified by the .ini files (e.g.,
    # `OutputDir`)  **MUST** be absolute paths, **NOT** relative.

    # First we specify the .ini files for the SAGE galaxy evolution.
    gal_ini_model1="/fred/oz004/jseiler/kali/self_consistent_output/rsage_constant/ini_files/const_0.20_SAGE.ini"
    gal_ini_model2="/fred/oz004/jseiler/kali/self_consistent_output/rsage_MHneg/ini_files/MHneg_1e7_1e12_0.99_0.05_SAGE.ini"   
    gal_ini_model3="/fred/oz004/jseiler/kali/self_consistent_output/rsage_MHpos/ini_files/MHpos_1e8_1e12_0.01_0.40_SAGE.ini"
    gal_ini_model4="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.30_beta0.00_SAGE.ini"
    gal_ini_model5="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha1.00_beta1.50_delta1.00_SAGE.ini"

    gal_ini_model6="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.30_beta0.12_SAGE.ini"
    gal_ini_model7="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.20_beta0.19_SAGE.ini"

    gal_ini_model8="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha4.50_beta0.50_delta1.00_SAGE.ini"
    gal_ini_model9="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha0.63_beta1.50_delta1.00_SAGE.ini"

    # Then the .ini files for cifog reionization.
    reion_ini_model1="/fred/oz004/jseiler/kali/self_consistent_output/rsage_constant/ini_files/const_0.20_cifog.ini"
    reion_ini_model2="/fred/oz004/jseiler/kali/self_consistent_output/rsage_MHneg/ini_files/MHneg_1e7_1e12_0.99_0.05_cifog.ini"   
    reion_ini_model3="/fred/oz004/jseiler/kali/self_consistent_output/rsage_MHpos/ini_files/MHpos_1e8_1e12_0.01_0.40_cifog.ini"
    reion_ini_model4="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.30_beta0.00_cifog.ini"
    reion_ini_model5="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha1.00_beta1.50_delta1.00_cifog.ini"

    reion_ini_model6="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.30_beta0.12_cifog.ini"
    reion_ini_model7="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.20_beta0.19_cifog.ini"

    reion_ini_model8="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha4.50_beta0.50_delta1.00_cifog.ini"
    reion_ini_model9="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha0.63_beta1.50_delta1.00_cifog.ini"


    # Higher tau
    gal_ini_model10="/fred/oz004/jseiler/kali/self_consistent_output/rsage_constant/ini_files/const_0.3_SAGE.ini"
    gal_ini_model11="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.45_beta0.00_SAGE.ini"
    gal_ini_model12="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha1.00_beta0.80_delta1.00_SAGE.ini"

    reion_ini_model10="/fred/oz004/jseiler/kali/self_consistent_output/rsage_constant/ini_files/const_0.3_cifog.ini"
    reion_ini_model11="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.45_beta0.00_cifog.ini"
    reion_ini_model12="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha1.00_beta0.80_delta1.00_cifog.ini"

    # All .ini files included in this array will be plotted.
    paper_models = 0
    const_fej_SFR = 0
    no_const = 0
    higher_tau = 0
    const = 0
    sliding_fesc = 0
    combined_higher_tau = 0
    UV = 1

    if paper_models:

        output_directory = "./paper_models"
        gal_ini_files = [gal_ini_model1,
                         gal_ini_model2, 
                         gal_ini_model3, 
                         gal_ini_model4, 
                         gal_ini_model5]

        reion_ini_files = [reion_ini_model1,
                           reion_ini_model2,
                           reion_ini_model3,
                           reion_ini_model4,
                           reion_ini_model5]

        model_tags = [r"$\mathbf{f_\mathrm{esc} = 0.20}$",
                      r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}^{-1}}$",
                      r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}}$",
                      r"$\mathbf{f_\mathrm{esc} \: \propto \: f_\mathrm{ej}}$",
                      r"$\mathbf{f_\mathrm{esc} \: \propto \: SFR}$"]


    if const_fej_SFR:
        output_directory = "./30month"
        gal_ini_files = [gal_ini_model1,
                         gal_ini_model4, 
                         gal_ini_model5]

        reion_ini_files = [reion_ini_model1,
                           reion_ini_model4,
                           reion_ini_model5]

        model_tags = [r"$\mathbf{f_\mathrm{esc} = 0.20}$",
                      r"$\mathbf{f_\mathrm{esc} \: \propto \: SN}$",
                      r"$\mathbf{f_\mathrm{esc} \: \propto \: SFR}$"]

    if no_const:
        output_directory = "./no_const"
        gal_ini_files = [gal_ini_model2, 
                         gal_ini_model3, 
                         gal_ini_model4, 
                         gal_ini_model5]

        reion_ini_files = [reion_ini_model2,
                           reion_ini_model3,
                           reion_ini_model4,
                           reion_ini_model5]

        model_tags = [r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}^{-1}}$",
                      r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}}$",
                      r"$\mathbf{f_\mathrm{esc} \: \propto \: f_\mathrm{ej}}$",
                      r"$\mathbf{f_\mathrm{esc} \: \propto \: SFR}$"]

    if higher_tau:
        output_directory = "./higher_tau"
        gal_ini_files = [gal_ini_model10, 
                         gal_ini_model11, 
                         gal_ini_model12]

        reion_ini_files = [reion_ini_model10,
                           reion_ini_model11,
                           reion_ini_model12]

        model_tags = [r"$\mathbf{f_\mathrm{esc} = 0.30}$",
                      r"$\mathbf{f_\mathrm{esc} \: \propto \: f_\mathrm{ej}}$",
                      r"$\mathbf{f_\mathrm{esc} \: \propto \: SFR}$"]


    if const:

        output_directory = "./const"
        gal_ini_files = [gal_ini_model1]
        reion_ini_files = [reion_ini_model1]

        model_tags = [r"$\mathbf{Fiducial \: \texttt{RSAGE}}$"]

    if sliding_fesc:

        output_directory = "./sliding_fesc"

        gal_ini_1="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.40_beta0.00_SAGE.ini"
        gal_ini_2="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.30_beta0.05_SAGE.ini"
        gal_ini_3="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.10_beta0.20_SAGE.ini"
        gal_ini_4="/fred/oz004/jseiler/kali/self_consistent_output/rsage_constant/ini_files/const_0.25_SAGE.ini"
        gal_ini_5="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha0.20_beta4.50_delta1.00_SAGE.ini"
        gal_ini_6="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha0.30_beta2.30_delta1.00_SAGE.ini"
        gal_ini_7="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha0.63_beta1.50_delta1.00_SAGE.ini"
        gal_ini_8="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha1.00_beta1.00_delta1.00_SAGE.ini"
        gal_ini_9="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha4.50_beta0.50_delta1.00_SAGE.ini"

        reion_ini_1="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.40_beta0.00_cifog.ini"
        reion_ini_2="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.30_beta0.05_cifog.ini"
        reion_ini_3="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.10_beta0.20_cifog.ini"
        reion_ini_4="/fred/oz004/jseiler/kali/self_consistent_output/rsage_constant/ini_files/const_0.25_cifog.ini"      
        reion_ini_5="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha0.20_beta4.50_delta1.00_cifog.ini"
        reion_ini_6="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha0.30_beta2.30_delta1.00_cifog.ini"
        reion_ini_7="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha0.63_beta1.50_delta1.00_cifog.ini"
        reion_ini_8="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha1.00_beta1.00_delta1.00_cifog.ini"
        reion_ini_9="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha4.50_beta0.50_delta1.00_cifog.ini"

        gal_ini_files = [gal_ini_1,
                         gal_ini_2,
                         gal_ini_3,
                         gal_ini_4,
                         gal_ini_5,
                         gal_ini_6,
                         gal_ini_7,
                         gal_ini_8,
                         gal_ini_9]

        reion_ini_files = [reion_ini_1,
                         reion_ini_2,
                         reion_ini_3,
                         reion_ini_4,
                         reion_ini_5,
                         reion_ini_6,
                         reion_ini_7,
                         reion_ini_8,
                         reion_ini_9]

        model_tags = [r"fej1",
                      r"fej2",
                      r"fej3",
                      r"Constant",
                      r"SFR1",
                      r"SFR2",
                      r"SFR3",
                      r"SFR4",
                      r"SFR5"]


    if combined_higher_tau:
        output_directory = "./combined_higher_tau_400hrs"

        gal_ini_files = [gal_ini_model1,
                         gal_ini_model4, 
                         gal_ini_model5,
                         gal_ini_model10, 
                         gal_ini_model11, 
                         gal_ini_model12]

        reion_ini_files = [reion_ini_model1,
                           reion_ini_model4,
                           reion_ini_model5,
                           reion_ini_model10,
                           reion_ini_model11,
                           reion_ini_model12]

        model_tags = [r"$\mathbf{f_\mathrm{esc} = 0.30}$",
                      r"$\mathbf{f_\mathrm{esc} \: \propto \: f_\mathrm{ej}}$",
                      r"$\mathbf{f_\mathrm{esc} \: \propto \: SFR}$"]

    if UV:
        output_directory = "./UV"

        gal_ini="/fred/oz004/jseiler/kali/self_consistent_output/rsage_UV/ini_files/test_UV_SAGE.ini"
        reion_ini="/fred/oz004/jseiler/kali/self_consistent_output/rsage_UV/ini_files/test_UV_cifog.ini"

        gal_ini_files = [gal_ini]
        reion_ini_files = [reion_ini]

        model_tags = [r"$Test$"]

    # ============================================================= #
    # Switches to control what plots to make. 0 to skip, 1 to plot. #
    # ============================================================= #

    # ==================== #
    # Galaxy Plot Toggles. # 
    # ==================== #
    gal_nion   = 0
    mstar_fesc = 0
    mstar_fej  = 0
    mstar_SFR  = 0
    SMF        = 1
    UVLF       = 1

    # ==================== #
    # Galaxy Plot Options. # 
    # ==================== #
    # For some plots, there is the option of plotting one model at different
    # snapshots (within one panel) or plotting all models at one snapshot
    # (within one panel).  
    plot_snaps_for_models = [[76],  # For each panel, plot one model 
                             [76],  # at specified snapshots.
                             [76],   
                             [76],   
                             [76],   
                             [76]]  
    #plot_snaps_for_models = None
    plot_models_at_snaps = [76] # For each panel, plot one specified snapshot 
                                # for all models.

    plot_single_panel = 1  # Instead of plotting multiple panels, plot all on a
                           # single panel.

    # For the stellar mass function, we may only want to plot at specific
    # redshifts for each model.
    SMF_plot_z = [[6.0, 7.0, 8.0]]

    # And same for the UV luminosity function.
    UVLF_plot_z = [[6.0, 7.0, 8.0]]

    # If you only want to use a specific range of files, specify here. If you
    # want to use all the files from a model, set these to ``None``.
    first_file = None 
    last_file = None
   
    # =============================================== #
    # Galaxy Plotting Time (Shouldn't need to touch). # 
    # =============================================== #
    galaxy_plots = {"nion" :                  gal_nion,
                    "mstar_fesc" :            mstar_fesc,
                    "SMF" :                   SMF,
                    "mstar_fej" :             mstar_fej,
                    "mstar_SFR" :             mstar_SFR,
                    "UVLF" :                  UVLF,
                    "plot_snaps_for_models" : plot_snaps_for_models,
                    "plot_models_at_snaps" :  plot_models_at_snaps,
                    "plot_single_panel" :     plot_single_panel,
                    "SMF_plot_z" :            SMF_plot_z,
                    "UVLF_plot_z" :           UVLF_plot_z,
                    "first_file" :            first_file,
                    "last_file" :             last_file}

    # Check if any galaxy plots need to be done.
    for field in galaxy_plots.keys():
        if galaxy_plots[field] == 1:
            galdata.plot_galaxy_properties(rank, size, comm, gal_ini_files,
                                           model_tags, galaxy_plots,
                                           output_directory, output_format)
            break
    exit()

    # =========================== #
    # Reionization Plot Toggles . # 
    # =========================== #
    single_slice         = 0
    single_ps            = 0
    history              = 0
    reion_nion           = 0
    ps_fixed_XHI         = 1
    optical_depth        = 0
    ps_scales            = 1
    ps_scales_beta       = 1
    slices_fixed_XHI     = 1
    bubble_size          = 0
    zreion_dens_cross    = 0
    dens_ion_contours    = 0
    dens_zreion_contours = 0

    # ========================== #
    # Reionization Plot Options. # 
    # ========================== #
    # `px_fixed_XHI` and `slices_fixed_XHI` are done at fixed neutral
    # fractions. `ps_scales` also has these marked on the plot. 
    #fixed_XHI_values = [0.90, 0.75, 0.50, 0.25, 0.10]
    fixed_XHI_values = [0.50]

    # `ps_scales` has the option to instead mark redshift instead of fixed
    # neutral fractions. Set this to `None` to use the fixed neutral fractions
    # above.
    ps_scales_z = [11.0, 10.0, 9.0, 8.0, 7.0]
    #ps_scales_z = None 

    # Neutral fractions that define the star/mid/end of reionization.
    duration_definition = [0.90, 0.50, 0.01]

    # What k value (h/Mpc) defines 'small' and 'large' scales?                                             
    small_scale_def = 1.0
    large_scale_def = 0.3

    # The uncertainties of various instruments.
    # Using Cath Trott's Marvelous Spreadsheet, assuming 400 hours integration
    # time with 10MHz bandwidth.

    # Select which instrument to use. Only accepted are "MWA", "SKA", "HERA", 
    # "All" or None. Selecting None will not plot any errors. 
    use_instrument = "SKA"

    # First define the number of hours we will be integrating for.
    integration_time = 200

    if use_instrument == "MWA":
        small_scale_err = 4.52e4   # 256 dish configuration at 1.0 h/Mpc.
        large_scale_err = 1.18e3  # 256 dish configuration at 0.3 h/Mpc.

    elif use_instrument == "SKA":
        small_scale_err = 1.37  # 1.0 h/Mpc. 
        large_scale_err = 1.15e-2  # 0.3 h/Mpc. 

    elif use_instrument == "HERA":
        # WARNING: HERA ONLY PROBES 0.0757 h/Mpc to 0.5604 h/Mpc.
        print("Warning!  HERA only probes 0.0757 h/Mpc to 0.5604 h/Mpc.")
        small_scale_err = 2.25  # 0.56 h/Mpc. 
        large_scale_err = 2.65e-2  # 0.3 h/Mpc.
    elif use_instrument is None:
        small_scale_err = None
        large_scale_err = None
    else:
        print("The only accepted values for 'use_instrument' are 'MWA', "
              "'SKA', 'HERA' or None.")
        raise ValueError

    # The error scales linearly with integration time. So we scale all our
    # values to an integration time of 400 hours.
    if use_instrument:
        small_scale_err *= (400.0 / integration_time)
        large_scale_err *= (400.0 / integration_time)

    # For plots of the slices of the ionization field, what grid index should
    # we start our cuts and what should the thickness (in grid cells) be?  All
    # models are normalized to have the same **spatial** dimensions as model 0. 
    cut_slice = 40
    cut_thickness = 1

    # Finally, if we want to sweep a parameter space worth of models and
    # construct contours of constant tau and duration. To grab all the .ini
    # files from a directory with specified alpha values, specify a directory
    # as `ini_dir` and the alpha/beta ranges of interest. 
    contours = 0
    if contours:

        # These ranges are **INCLUSIVE**
        alpha_low = 0.0
        alpha_high = 0.60
        alpha_step = 0.10     
        alpha_vals = np.arange(alpha_low, alpha_high + alpha_step, alpha_step)

        beta_low = 0.0
        beta_high = 0.30
        beta_step = 0.05     
        beta_vals = np.arange(beta_low, beta_high + beta_step, beta_step)

        ini_dir = "/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/"
        gal_ini_files, reion_ini_files = misc.get_ini_from_dir(ini_dir,
                                                               alpha_vals=alpha_vals,
                                                               beta_vals=beta_vals)
        alpha_beta_limits = [[alpha_low, alpha_high, alpha_step],
                             [beta_low, beta_high, beta_step]]
    else:
        alpha_beta_limits = None

    # ===================================================== #
    # Reionization Plotting Time (Shouldn't need to touch). # 
    # ===================================================== #
    reion_plots = {"single_slice" :         single_slice,
                   "single_ps" :            single_ps,
                   "history" :              history,
                   "nion" :                 reion_nion,
                   "ps_fixed_XHI" :         ps_fixed_XHI,
                   "contours" :             contours, 
                   "optical_depth" :        optical_depth,
                   "ps_scales" :            ps_scales,
                   "ps_scales_beta" :       ps_scales_beta,
                   "bubble_size" :          bubble_size,
                   "slices_fixed_XHI" :     slices_fixed_XHI,
                   "zreion_dens_cross" :    zreion_dens_cross,
                   "dens_ion_contours" :    dens_ion_contours,
                   "dens_zreion_contours" : dens_zreion_contours,
                   "fixed_XHI_values" :     fixed_XHI_values,
                   "ps_scales_z" :          ps_scales_z,
                   "alpha_beta_limits" :    alpha_beta_limits,
                   "duration_definition" :  duration_definition,
                   "small_scale_def" :      small_scale_def,
                   "large_scale_def" :      large_scale_def,
                   "small_scale_err" :      small_scale_err,
                   "large_scale_err" :      large_scale_err,
                   "cut_slice" :            cut_slice,
                   "cut_thickness" :        cut_thickness}

    # Check if any reionization plots need to be done.
    for field in reion_plots.keys():
        if reion_plots[field] == 1:
            reiondata.plot_reion_properties(rank, size, comm, reion_ini_files,
                                            gal_ini_files, model_tags,
                                            reion_plots, output_directory,
                                            output_format)
            break
