#!/usr/bin/env python
from __future__ import print_function

import os
import matplotlib
matplotlib.use('Agg')

import pylab as plt
import numpy as np
import scipy.integrate as integrate

import PlotScripts
import ReadScripts
import AllVars

PlotScripts.Set_Params_Plot()

output_format = ".png"

def calculate_tau(base_XHII_fname, base_density_fname, snaps):

    def integrand(z):
        H = AllVars.Hubble_Param(z, AllVars.Hubble_h, AllVars.Omega_m) / (AllVars.pc_to_m * 1.0e6 / 1.0e3) # Hubble Parameter in Mpc / s / Mpc. 
        return (((1 + z)**2) / H)

    tau_04 = integrate.quad(integrand, 0, 4)[0] # Tau for z = 0 to 4.
    tau_04 *= (1 + 2*AllVars.Y/(4 * (1-AllVars.Y)))

    tau_46 = integrate.quad(integrand, 4, AllVars.SnapZ[snaps[-1]])[0] # Tau for z = 4 to 6.
    tau_46 *= (1 + AllVars.Y/(4* (1-AllVars.Y)))

    tau_06 = tau_04 + tau_46

    tau = np.empty(len(snaps))
 
    for count, snap in enumerate(reversed(snaps)):
        XHII_fname = "{0}_XHII_{1:03d}".format(base_XHII_fname, snap)
        XHII = ReadScripts.read_binary_grid(XHII_fname, 256, 2)
       
        density_fname = "{0}snap{1:03d}.dens.dat".format(base_density_fname, snap)
        density = ReadScripts.read_binary_grid(density_fname, 256, 2)

        mass_fraction = np.sum(XHII * density / np.sum(density))
	
        H = AllVars.Hubble_Param(AllVars.SnapZ[snap], AllVars.Hubble_h, AllVars.Omega_m) / (AllVars.pc_to_m * 1.0e6 / 1.0e3) # Hubble Parameter in Mpc / s / Mpc.	
        numerator = ((1 + AllVars.SnapZ[snap]) **2) * mass_fraction 

        idx = len(snaps)-count-1
        if count == 0:
            tau[idx] = tau_06 + (( numerator / H) * (AllVars.SnapZ[snap] - AllVars.SnapZ[snaps[-1]]) * (1 + AllVars.Y/(4 * (1 - AllVars.Y)))) 
        else:
            tau[idx] = tau[idx+1] + (( numerator / H) * (AllVars.SnapZ[snap] - AllVars.SnapZ[snap+1]) * (1 + AllVars.Y/(4 * (1-AllVars.Y)))) 


    tau *= AllVars.n_HI(0, AllVars.Hubble_h, AllVars.Omega_b, AllVars.Y) * AllVars.c_in_ms * AllVars.Sigmat

    return tau[0]


def plot_contours(base_XHII_fname, base_density_fname, alpha, beta, snaps,
                  output_tag="tau_contours"):


    tau = np.zeros((len(alpha),len(beta)))

    for alpha_count, alpha_val in enumerate(alpha):
        for beta_count, beta_val in enumerate(beta):
            fname = "{0}_alpha{1}_beta{2}".format(base_fname, alpha_val,
                                                  beta_val)
            tau[alpha_count, beta_count] = calculate_tau(fname, 
                                                         base_density_fname, 
                                                         snaps)
            print("alpha {0}\tbeta {1}\tTau {2}".format(alpha_val, beta_val,
                                                        tau[alpha_count, beta_count]))

    fig1 = plt.figure(111)
    ax1 = fig1.add_subplot(111)

    X, Y = np.meshgrid(beta, alpha)

    CS = ax1.contour(Y, X, tau) 
    ax1.clabel(CS, inline=1, fontsize=10)

    ax1.set_xlabel(r"$\mathbf{\alpha}$", size = PlotScripts.global_fontsize)
    ax1.set_ylabel(r"$\mathbf{\beta}$", size = PlotScripts.global_fontsize)

    outputFile = "./{0}{1}".format(output_tag, output_format) 
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))
    plt.close()

if __name__ == '__main__':

    AllVars.Set_Params_Kali()
    AllVars.Set_Constants()

    alpha = [0.2, 0.3, 0.4, 0.5, 0.6]
    beta = [0.0, 0.05, 0.10, 0.15]
    density_fname="/fred/oz004/jseiler/kali/density_fields/1024_subsampled_256/"
    base_fname="/fred/oz004/jseiler/kali/self_consistent_output/shifted_fej/grids/cifog/shifted_fej"
    plot_contours(base_fname, density_fname, alpha, beta, np.arange(28,98),
                  output_tag="tau_contours_shifted")
