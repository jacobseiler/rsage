#!/usr/bin/env python from __future__ import print_function

import os
import matplotlib
matplotlib.use('Agg')

import pylab as plt
import numpy as np
import scipy.integrate as integrate

import PlotScripts as ps
import ReadScripts as rs
import AllVars as av

ps.Set_Params_Plot()

output_format = ".png"

def calculate_tau(base_XHII_fname, base_density_fname, snaps):

    def integrand(z):
        H = av.Hubble_Param(z, av.Hubble_h, av.Omega_m) / (av.pc_to_m * 1.0e6 / 1.0e3) # Hubble Parameter in Mpc / s / Mpc. 
        return (((1 + z)**2) / H)

    tau_04 = integrate.quad(integrand, 0, 4)[0] # Tau for z = 0 to 4.
    tau_04 *= (1 + 2*av.Y/(4 * (1-av.Y)))

    tau_46 = integrate.quad(integrand, 4, av.SnapZ[snaps[-1]])[0] # Tau for z = 4 to 6.
    tau_46 *= (1 + av.Y/(4* (1-av.Y)))

    tau_06 = tau_04 + tau_46

    tau = np.empty(len(snaps))
 
    for count, snap in enumerate(reversed(snaps)):
        XHII_fname = "{0}_XHII_{1:03d}".format(base_XHII_fname, snap)
        XHII = rs.read_binary_grid(XHII_fname, 256, 2)
       
        density_fname = "{0}{1:03d}.dens.dat".format(base_density_fname, snap)
        density = rs.read_binary_grid(density_fname, 256, 2)

        mass_fraction = np.sum(XHII * density / np.sum(density))
	
        H = av.Hubble_Param(av.SnapZ[snap], av.Hubble_h, av.Omega_m) / (av.pc_to_m * 1.0e6 / 1.0e3) # Hubble Parameter in Mpc / s / Mpc.	
        numerator = ((1 + av.SnapZ[snap]) **2) * mass_fraction 

        idx = len(snaps)-count-1
        if count == 0:
            tau[idx] = tau_06 + (( numerator / H) * (av.SnapZ[snap] - av.SnapZ[snaps[-1]]) * (1 + av.Y/(4 * (1 - av.Y)))) 
        else:
            tau[idx] = tau[idx+1] + (( numerator / H) * (av.SnapZ[snap] - av.SnapZ[snap+1]) * (1 + av.Y/(4 * (1-av.Y)))) 


    tau *= av.n_HI(0, av.Hubble_h, av.Omega_b, av.Y) * av.c_in_ms * av.Sigmat

    return tau[0]


def fit_tau(alpha, beta, tau=None, tau_file=None):

    from scipy.optimize import curve_fit

    def func(X, a, b, c):
        alpha = X[0, :]
        beta = X[1, :]
        return alpha*a + beta*b +c

    alpha_arr = []
    beta_arr = []

    for alpha_val in alpha:
        for beta_val in beta:
            alpha_arr.append(alpha_val)
            beta_arr.append(beta_val)

    if tau_file:
        tau = np.loadtxt(tau_file)
        tau = tau.ravel()

    alpha = np.array(alpha_arr)
    beta = np.array(beta_arr)
    tau = np.array(tau)

    X = np.array([alpha, beta])

    popt, pcov = curve_fit(func, X, tau)
    print("Function is of the form alpha*a + beta*b +c")
    print("Best fit values are {0}".format(popt))
    print("Covariance is {0}".format(pcov))
    print("Error on fits is {0}".format(np.sqrt(np.diag(pcov))))

    return


def calculate_all_taus(XHII_basename, density_basename, alpha, beta, snaps):

    tau = np.zeros((len(alpha),len(beta)))

    for alpha_count, alpha_val in enumerate(alpha):
        for beta_count, beta_val in enumerate(beta):
            if beta_val == 0.20:
                fname = "{0}_alpha{1:.2f}_beta{2:.2f}".format(XHII_basename, alpha_val,
                                                              beta_val)
            elif (alpha_val == 0.30 and beta_val == 0.1):
                fname = "{0}_alpha{1:.2f}_beta{2:.2f}".format(XHII_basename, alpha_val,
                                                          beta_val)
            else:
                fname = "{0}_alpha{1:.2f}_beta{2}".format(XHII_basename, alpha_val,
                                                          beta_val)
            tau[alpha_count, beta_count] = calculate_tau(fname, 
                                                         density_basename, 
                                                         snaps)
            print("alpha {0}\tbeta {1}\tTau {2}".format(alpha_val, beta_val,
                                                        tau[alpha_count, beta_count]))

    header = "alpha {0} beta {1}".format(alpha, beta)
    np.savetxt("./tau_values.txt", tau, header=header)

    return


def plot_tau(alpha, beta, output_tag, tau=None, tau_file=None):

    fig1 = plt.figure(111)
    ax1 = fig1.add_subplot(111)

    if tau_file:
        tau = np.loadtxt(tau_file)        

    X, Y = np.meshgrid(beta, alpha)

    CS = ax1.contour(Y, X, tau) 
    ax1.clabel(CS, inline=1, fontsize=10)

    ax1.set_xlabel(r"$\mathbf{\alpha}$", size = ps.global_fontsize)
    ax1.set_ylabel(r"$\mathbf{\beta}$", size = ps.global_fontsize)

    outputFile = "./{0}{1}".format(output_tag, output_format) 
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))
    plt.close()


if __name__ == '__main__':

    av.Set_Params_Kali()
    av.Set_Constants()

    alpha = [0.2, 0.3, 0.4, 0.6]
    beta = [0.0, 0.05, 0.10, 0.15, 0.20]

    fit_tau(alpha, beta, tau_file="./tau_values.txt")
    plot_tau(alpha, beta, "new_contours", tau_file="./tau_values.txt")
    exit()    
    
    cifog_ini_name = "/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/ini_files/fej_alpha0.40_beta0.0_cifog.ini"

    cifog_params, _ = rs.read_cifog_ini(cifog_ini_name)

    first_snap = int(cifog_params["SimulationLowSnap"])
    last_snap = int(cifog_params["SimulationHighSnap"])
    snaps = np.arange(first_snap+1, last_snap)
    density_basename = cifog_params["inputIgmDensityFile"]

    XHII_basename = "/fred/oz004/jseiler/kali/self_consistent_output/rsage_fej/grids/cifog/fej"

    calculate_all_taus(XHII_basename, density_basename, alpha, beta, snaps)
