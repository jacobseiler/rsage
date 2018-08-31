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
Version: 0.2
"""

from __future__ import print_function

import GalaxyData as galdata
import ReionData as reiondata
import MiscData as misc
import AllVars as av

import numpy as np 
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

output_format = "png"


if __name__ == "__main__":

    av.Set_Constants()

    output_directory = "./31st_Aug"

    # Plotting is driven entirely through specifying the ``.ini`` files. 
    # For this reason, the directories in the ``.ini`` files **MUST** be
    # absolute paths, **NOT** relative. 
    gal_ini_model1="/fred/oz004/jseiler/kali/self_consistent_output/rsage_constant/ini_files/const_0.3_SAGE.ini"
    gal_ini_model2="/fred/oz004/jseiler/kali/self_consistent_output/rsage_MHneg/ini_files/MHneg_1e8_1e12_0.90_0.05_SAGE.ini"   
    gal_ini_model3="/fred/oz004/jseiler/kali/self_consistent_output/rsage_MHpos/ini_files/MHpos_1e8_1e12_0.01_0.50_SAGE.ini"
    gal_ini_model4="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.40_beta0.05_SAGE.ini"
    gal_ini_model5="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha1.00_beta1.00_delta1.00_SAGE.ini"
    gal_ini_model6="/home/jseiler/rsage/ini_files/kali_noreion_SAGE.ini"
    gal_ini_model7="/fred/oz004/jseiler/kali/self_consistent_output/rsage_myMHneg/ini_files/MHneg_1e8_1e12_0.95_0.05_SAGE.ini"

    reion_ini_model1="/fred/oz004/jseiler/kali/self_consistent_output/rsage_constant/ini_files/const_0.3_cifog.ini"
    reion_ini_model2="/fred/oz004/jseiler/kali/self_consistent_output/rsage_MHneg/ini_files/MHneg_1e8_1e12_0.90_0.05_cifog.ini"   
    reion_ini_model3="/fred/oz004/jseiler/kali/self_consistent_output/rsage_MHpos/ini_files/MHpos_1e8_1e12_0.01_0.50_cifog.ini"
    reion_ini_model4="/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.40_beta0.05_cifog.ini"
    reion_ini_model5="/fred/oz004/jseiler/kali/self_consistent_output/rsage_SFR/ini_files/SFR_alpha1.00_beta1.00_delta1.00_cifog.ini"
    reion_ini_model6="/home/jseiler/rsage/ini_files/kali_noreion_cifog.ini"
    reion_ini_model7="/fred/oz004/jseiler/kali/self_consistent_output/rsage_myMHneg/ini_files/MHneg_1e8_1e12_0.95_0.05_cifog.ini"

    # All ``.ini`` files included in this array will be plotted.
    gal_ini_files = [gal_ini_model1,
                     gal_ini_model2, 
                     gal_ini_model3, 
                     gal_ini_model4, 
                     gal_ini_model5]

    #gal_ini_files = [gal_ini_model1]

    reion_ini_files = [reion_ini_model1,
                       reion_ini_model2, 
                       reion_ini_model3, 
                       reion_ini_model4, 
                       reion_ini_model5]

    #reion_ini_files = [reion_ini_model1]

    ini_dir = "/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/"

    # These ranges are **INCLUSIVE**
    alpha_low = 0.0
    alpha_high = 0.60
    alpha_step = 0.10     
    alpha_vals = np.arange(alpha_low, alpha_high + alpha_step, alpha_step)

    # These ranges are **INCLUSIVE**
    beta_low = 0.0
    beta_high = 0.30
    beta_step = 0.05     
    beta_vals = np.arange(beta_low, beta_high + beta_step, beta_step)

    #gal_ini_files, reion_ini_files = misc.get_ini_from_dir(ini_dir,
    #                                                       alpha_vals=alpha_vals,
    #                                                       beta_vals=beta_vals)

    #gal_ini_files = gal_ini_files[0:3]
    #reion_ini_files = reion_ini_files[0:3]
    # These are the labels that will appear on the axis legends for each model.
    model_tags = [r"$\mathbf{f_\mathrm{esc} = 0.30}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}^{-1}}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: f_\mathrm{ej}}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: SFR}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: M_\mathrm{H}^{-1} \: Mine}$"]

    '''
    model_tags = [r"$\mathbf{f_\mathrm{esc} \: \propto \: f_\mathrm{ej} \: \alpha = 0.00, \beta = 0.05}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: f_\mathrm{ej} \: \alpha = 0.10, \beta = 0.05}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: f_\mathrm{ej} \: \alpha = 0.20, \beta = 0.05}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: f_\mathrm{ej} \: \alpha = 0.30, \beta = 0.05}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: f_\mathrm{ej} \: \alpha = 0.40, \beta = 0.05}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: f_\mathrm{ej} \: \alpha = 0.50, \beta = 0.05}$",
                  r"$\mathbf{f_\mathrm{esc} \: \propto \: f_\mathrm{ej} \: \alpha = 0.60, \beta = 0.05}$"]
    '''
    # ============================================================= #
    # Switches to control what plots to make. 0 to skip, 1 to plot. #
    # ============================================================= #

    ## Galaxy Plots ##
    gal_nion = 0
    mstar_fesc = 0
    SMF = 0
    mstar_fej = 0
    mstar_SFR = 0
    mstar_infall = 0

    # For some plots, there is the option of plotting one model at different
    # snapshots (within one panel) or plotting all models at one snapshot
    # (within one panel).  
    plot_snaps_for_models = [[33, 50, 76, 93],  # For each panel, plot one model 
                             [33, 50, 76, 93],  # at specified snapshots.
                             [33, 50, 76, 93],   
                             [33, 50, 76, 93],   
                             [33, 50, 76, 93],   
                             [33, 50, 76, 93]]  

    plot_models_at_snaps = None  # For each panel, plot one specified snapshot 
                                 # for all models.

    ## Reionization Plots ##
    history = 0
    reion_nion = 0
    ps_fixed_XHI = 0
    duration_contours = 0
    optical_depth = 0
    ps_scales = 1

    fixed_XHI_values = [0.90, 0.75, 0.50, 0.25, 0.10]
    duration_contours_limits = [[alpha_low, alpha_high, alpha_step],  # This is the min/max/step
                                [beta_low, beta_high, beta_step]] # for the alpha/beta values. 
    duration_definition = [0.90, 0.50, 0.01]  # Neutral fraction that defines the
                                              # start/mid/end of reionization. 

    small_scale_def = 1.0
    large_scale_def = 0.2

    # ======================= #
    # Don't touch below here. # 
    # ======================= #
    galaxy_plots = {"nion" : gal_nion,
                    "mstar_fesc" : mstar_fesc,
                    "SMF" : SMF,
                    "mstar_fej" : mstar_fej,
                    "mstar_SFR" : mstar_SFR,
                    "mstar_infall" : mstar_infall,
                    "plot_snaps_for_models" : plot_snaps_for_models,
                    "plot_models_at_snaps" : plot_models_at_snaps}

    # Check if any galaxy plots need to be done.
    for field in galaxy_plots.keys():
        if galaxy_plots[field] == 1:
            galdata.plot_galaxy_properties(rank, size, comm, gal_ini_files,
                                           model_tags, galaxy_plots, output_directory)
            break

    reion_plots = {"history" : history,
                   "nion" : reion_nion,
                   "ps_fixed_XHI" : ps_fixed_XHI,
                   "duration_contours" : duration_contours,
                   "optical_depth" : optical_depth,
                   "ps_scales" : ps_scales,
                   "fixed_XHI_values" : fixed_XHI_values,
                   "duration_contours_limits" : duration_contours_limits,
                   "duration_definition" : duration_definition,
                   "small_scale_def" : small_scale_def,
                   "large_scale_def" : large_scale_def}

    # Check if any reionization plots need to be done.
    for field in reion_plots.keys():
        if reion_plots[field] == 1:
            reiondata.plot_reion_properties(rank, size, comm, reion_ini_files,
                                            gal_ini_files, model_tags,
                                            reion_plots, output_directory)
            break
