#!/usr/bin/env python
from __future__ import print_function

import os
import matplotlib
matplotlib.use('Agg')

import pylab as plt
import numpy as np
from numpy.fft import fftn, ifftn

import scipy.integrate as integrate
from scipy import stats

import PlotScripts
import ReadScripts
import AllVars


def calculate_HI_frac(XHII, density):
    """
    Calculates the mass-weighted fraction of ionized hydrogen for a given 
    ionization grid.

    Parameters
    ---------
    XHII: 3-Dimensional Array of floats. Required. 
        Grid that contains the fraction of ionized hydrogen (XHII) in each 
        cell. 

    density: 3-Dimensional Array of floats. Required. 
        Grid that contains the overdensity (rho/<rho>) of dark matter in each
        cell. 
 
    Returns
    -------
    HI: Float.
        Fraction of ionized hydrogen. 

    Units
    -----
    XHII and HI are unitless.
    Density is unitless (overdensity, rho/<rho>).    
    """

    HI = 1.0 - np.sum(XHII * density / np.sum(density))
       
    print("")
    print("Mass averaged HI fraction is {0:.4f}".format(HI))
   
    return HI 


def determine_close_idx(fname_HII, fname_density, SnapList, GridSize, 
                        precision, target_XHI_fraction, model_tags):

    XHII_fraction = np.zeros_like(SnapList, dtype=np.float32)

    for model_number in range(len(fname_HII)):
        for snapnum in range(len(SnapList[model_number])):
            HII_fname = "{0}_{1:03d}".format(fname_HII[model_number], 
                                              SnapList[model_number][snapnum])

            HII = ReadScripts.read_binary_grid(HII_fname,
                                               GridSize[model_number],
                                               precision[model_number])

            density_fname = "{0}{1:03d}.dens.dat".format(fname_density[model_number], 
                                                         SnapList[model_number][snapnum])
 
            density = ReadScripts.read_binary_grid(density_fname,
                                                   GridSize[model_number],
                                                   precision[model_number])

            HI_frac = calculate_HI_frac(HII, density)
            XHII_fraction[model_number][snapnum] = HI_frac

    for model_number in range(len(fname_HII)):
        print("Model {0}".format(model_tags[model_number])) 
        for val in target_XHI_fraction: 
            idx = (np.abs(XHII_fraction[model_number] - val)).argmin()
            print("HI Fract {0}: Nearest Idx {1} with value {2}".format(val, 
                                                                        idx, 
                                                                        XHII_fraction[model_number][idx]))
