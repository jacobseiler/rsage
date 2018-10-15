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

    SnapList = []
    for model_number in range(len(fname_HII)):
        SnapList.append([])
        print("Model {0}".format(model_tags[model_number])) 
        for val in target_XHI_fraction: 
            idx = (np.abs(XHII_fraction[model_number] - val)).argmin()
            print("HI Fract {0}: Nearest Idx {1} with value {2}".format(val, 
                                                                        idx, 
                                                                        XHII_fraction[model_number][idx]))
            SnapList[model_number].append(idx)

    return SnapList 


def determine_MH_fesc_constants(MH_low, MH_high, fesc_low, fesc_high):
    
    log_A = (np.log10(fesc_high) - (np.log10(fesc_low)*np.log10(MH_high)/np.log10(MH_low))) * pow(1 - (np.log10(MH_high) / np.log10(MH_low)), -1)
    B = (np.log10(fesc_low) - log_A) / np.log10(MH_low)
    A = pow(10, log_A)

    return A, B


def plot_Anne_MH(MH_low, MH_high, fesc_low, fesc_high, pos_scaling, ax1):

    halomass = np.arange(7.0, 18.0, 0.01)
    halomass = pow(10, halomass)

    if pos_scaling:
        fesc = 1.0 - pow((1.0 - fesc_low) * (1.0-fesc_low)/(1.0-fesc_high), -np.log10(halomass/MH_low)/np.log10(MH_high/MH_low))
        fesc[fesc < fesc_low] = fesc_low
        fesc[fesc > fesc_high] = fesc_high
    else:
        fesc = pow(fesc_low * fesc_low/fesc_high, -np.log10(halomass/MH_low)/np.log10(MH_high/MH_low))
        fesc[fesc > fesc_low] = fesc_low
        fesc[fesc < fesc_high] = fesc_high

    ax1.plot(np.log10(halomass),
             fesc, ls = '-', color = 'k',
             label = "Anne")

    ax1.set_xlabel("Halo Mass [Msun]")
    ax1.set_ylabel("fesc")

    ax1.set_ylim([0.0, 1.1])

    return ax1


def plot_my_MH(MH_low, MH_high, fesc_low, fesc_high, ax1):

    halomass = np.arange(7.0, 18.0, 0.01)
    halomass = pow(10, halomass)

    alpha, beta = determine_MH_fesc_constants(MH_low, MH_high,
                                              fesc_low, fesc_high)

    print("Alpha = {0} Beta = {1}".format(alpha, beta))
    fesc = alpha*pow(halomass,beta)

    if fesc_low > fesc_high:
        fesc[fesc > fesc_low] = fesc_low
        fesc[fesc < fesc_high] = fesc_high
    else:
        fesc[fesc < fesc_low] = fesc_low
        fesc[fesc > fesc_high] = fesc_high

    ax1.plot(np.log10(halomass),
             fesc, ls = '--', color = 'r',
             label = "Mine")

    ax1.set_xlabel("Halo Mass [Msun]")
    ax1.set_ylabel("fesc")

    ax1.set_ylim([0.0, 1.1])

    return ax1

def plot_MHs(MH_low, MH_high, fesc_low, fesc_high):

    fig1 = plt.figure(figsize = (8,8))
    ax1 = fig1.add_subplot(111)

    if fesc_high > fesc_low:
        pos_scaling = 1
    else:
        pos_scaling = 0
    ax1 = plot_Anne_MH(MH_low, MH_high, fesc_low, fesc_high, pos_scaling, ax1)
    ax1 = plot_my_MH(MH_low, MH_high, fesc_low, fesc_high, ax1)

    leg = ax1.legend(loc='upper right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(10)

    outputFile1 = "./fescMH_{0:.2e}_{1:.2e}_{2}_{3}.png".format(MH_low,
                                                                MH_high,
                                                                fesc_low,
                                                                fesc_high)
    fig1.savefig(outputFile1, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile1))
    plt.close(fig1)


def plot_SFR_fesc(alpha, beta, delta):


    fig1 = plt.figure(figsize = (8,8))
    ax1 = fig1.add_subplot(111)

    SFR = np.arange(-5, 2, 0.01)

    for alpha_val, beta_val, delta_val in zip(alpha, beta, delta):
        fesc = delta_val / (1.0 + np.exp(-alpha_val*(SFR-beta_val)))

        label = r"$\alpha = " + str(alpha_val) + r", \beta = " + str(beta_val) +\
                r", \delta = " + str(delta_val) + "$"
        
        print(label)

        ax1.plot(SFR, fesc, label=label)
    
    leg = ax1.legend(loc='upper left', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(10)

    outputFile1 = "./fesc_SFR.png"

    fig1.savefig(outputFile1, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile1))
    plt.close(fig1)

if __name__ == "__main__":

    MH_low = 1.0e8
    MH_high = 1.0e12
    fesc_high = 0.05

    #for fesc_low in [0.95]:
    #    plot_MHs(MH_low, MH_high, fesc_low, fesc_high)

    alpha = [0.2, 0.3, 0.63, 1.0, 4.50]
    beta = [4.5, 2.3, 1.5, 1.0, 0.5]
    delta = [1.0, 1.0, 1.0, 1.0, 1.0]
    plot_SFR_fesc(alpha, beta, delta)
