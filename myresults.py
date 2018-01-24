#!/usr/bin/env python
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')

import os
import heapq
import numpy as np
import pylab as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from numpy import *
from random import sample, seed, randint
from os.path import getsize as getFileSize
import math
import random
import csv
from cycler import cycler
from io import StringIO
#np.set_printoptions(threshold=np.nan)
from collections import Counter
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import AxesGrid
from astropy import units as u
from astropy import cosmology

import matplotlib.ticker as mtick
import PlotScripts
import ReadScripts
import AllVars

from mpi4py import MPI
from tqdm import tqdm

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

AllVars.Set_Constants()
AllVars.Set_Params_Mysim()
PlotScripts.Set_Params_Plot()

cosmo = cosmology.FlatLambdaCDM(H0 = AllVars.Hubble_h*100, Om0 = AllVars.Omega_m) 
t_BigBang = cosmo.lookback_time(100000).value # Lookback time to the Big Bang in Gyr.

output_format = ".png"

# For the Tiamat extended results there is a weird hump when calculating the escape fraction.
# This hump occurs at a halo mass of approximately 10.3. 
# The calculation of fesc skips this hump range (defined from kink_low to kink_high)
kink_low = 10.3
kink_high = 10.30000001

m_low = 7.0 # We only sum the photons coming from halos within the mass range m_low < Halo Mass < m_high
m_high = 15.0

m_gal_low = 4.8
m_gal_high = 12.0

m_low_SAGE = pow(10, m_low)/1.0e10 * AllVars.Hubble_h
m_high_SAGE = pow(10, m_high)/1.0e10 * AllVars.Hubble_h

bin_width = 0.2
NB = int((m_high - m_low) / bin_width)
NB_gal = int((m_gal_high - m_gal_low) / bin_width)

def raise_list_power(my_list, n):
    return [pow(x, n) for x in my_list]
def raise_power_list(my_list, n):
    return [pow(n, x) for x in my_list]

def calculate_beta(MUV, z):
    ''' 
    Calculation of the dust attenuation parameter Beta. Fit values are from Bouwens (2015) ApJ 793, 115.
    For z = 5 and 6, Bouwens uses a piece-wise linear relationship and a linear relationship for higher redshift. ##

    Parameters
    ----------
        MUV : `float'
        A value of the absolute magnitude in the UV (generally M1600) in the AB magnitude system.

    z : `float' 
        Redshift the attenuation is calculated at.

    Returns
    ------
    beta : `float'
        Value of the UV continuum paramaeter beta. 
    '''

    if (z >= 4.5 and z < 5.5): # z = 5 fits.
        if (MUV > -18.8):
            dB = -0.08
        else:
            dB = -0.17
        B = -2.05
        offset = 18.8
    elif (z >= 5.5 and z < 6.5): # z = 6 fits.
        if (MUV > -18.8):
            dB = -0.08
        else:
            dB = -0.24
        B = -2.22
        offset = 18.8

    elif (z >= 6.5 and z < 7.5): # z = 7 fits.
        dB = -0.20
        B = -2.05
        offset = 19.5
    elif (z >= 7.5 and z < 8.5): # z = 8 fits.
        dB = -0.15
        B = -2.13
        offset = 19.5
    elif (z >= 8.5 and z < 9.5): # z = 9 fits.
        dB = -0.16
        B = -2.19
        offset = 19.5
    elif (z >= 9.5 and z < 10.5): # z = 10 fits.
        dB = -0.16
        B = -2.16
        offset = 19.5

    beta = dB * (MUV + offset) + B

    return beta
        
def multiply(array):
    '''
    Performs element wise multiplication.

    Parameters
    ----------
    array : `~numpy.darray'
        The array to be multiplied.

    Returns
    -------
    total : `float'
        Total of the elements multiplied together.
    '''

    total = 1
    for i in range(0, len(array)):
        total *= array[i]
    return total

##

def Sum_Log(array):
    '''
    Performs an element wise sum of an array who's elements are in log-space.

    Parameters
    ----------
    array : array
    Array with elements in log-space.

    Returns
    ------
    sum_total : float
    Value of the elements taken to the power of 10 and summed.
    Units
    -----
    All units are kept the same as the inputs.
    '''

    sum_total = 0.0
    for i in range(0, len(array)):
        sum_total += 10**array[i]

    return sum_total

##

def Std_Log(array, mean):
    '''
    Calculates the standard deviation of an array with elements in log-space. 

    Parameters
    ----------
    array : array
    Array with elements in log-space.
    mean : float
    Mean of the array (not in log).

    Returns
    ------
    std : float
    Standard deviation of the input array taken to the power of 10. 
    Units
    -----
    All units are kept the same as the inputs.
    '''

    sum_total = 0.0
    for i in range(0, len(array)):
        sum_total += (10**array[i] - mean)**2

    sum_total *= 1.0/len(array)

    std = np.sqrt(sum_total)
    return std

##

def calculate_pooled_stats(mean_pool, std_pool, mean_local, std_local, N_local):
    '''
    Calculates the pooled mean and standard deviation from multiple processors and appends it to an input array.
    Formulae taken from https://en.wikipedia.org/wiki/Pooled_variance
    As we only care about these stats on the rank 0 process, we make use of junk inputs/outputs for other ranks.

    NOTE: Since the input data may be an array (e.g. pooling the mean/std for a stellar mass function).

    Parameters
    ----------
    mean_pool, std_pool : array of floats.
        Arrays that contain the current pooled means/standard deviation (for rank 0) or just a junk input (for other ranks).
    mean_local, mean_std : float or array of floats.
        The non-pooled mean and standard deviation unique for each process.
    N_local : int or array of ints. 
        Number of data points used to calculate the mean/standard deviation that is going to be added to the pool.

    Returns
    -------
    mean_pool, std_pool : array of floats.
        Original array with the new pooled mean/standard deviation appended (for rank 0) or the new pooled mean/standard deviation only (for other ranks).

    Units
    -----
    All units are the same as the input.
    All inputs MUST BE real-space (not log-space).
    '''

    if isinstance(mean_local, list) == True:    
        if len(mean_local) != len(std_local):
            print("len(mean_local) = {0} \t len(std_local) = {1}".format(len(mean_local), len(std_local)))
            raise ValueError("Lengths of mean_local and std_local should be equal")
   
    if ((type(mean_local).__module__ == np.__name__) == True or (isinstance(mean_local, list) == True)): # Checks to see if we are dealing with arrays. 
    
        N_times_mean_local = np.multiply(N_local, mean_local)
        N_times_var_local = np.multiply(N_local, np.multiply(std_local, std_local))
        
        N_local = np.array(N_local).astype(float)
        N_times_mean_local = np.array(N_times_mean_local).astype(np.float32)

        if rank == 0: # Only rank 0 holds the final arrays so only it requires proper definitions.
            N_times_mean_pool = np.zeros_like(N_times_mean_local) 
            N_pool = np.zeros_like(N_local)
            N_times_var_pool = np.zeros_like(N_times_var_local)

            N_times_mean_pool = N_times_mean_pool.astype(np.float64) # Recast everything to double precision then use MPI.DOUBLE.
            N_pool = N_pool.astype(np.float64)
            N_times_var_pool = N_times_var_pool.astype(np.float64)
        else:
            N_times_mean_pool = None
            N_pool = None
            N_times_var_pool = None

        comm.Barrier()

        N_times_mean_local = N_times_mean_local.astype(np.float64)
        N_local = N_local.astype(np.float64)
        N_times_var_local = N_times_var_local.astype(np.float64)

        comm.Reduce([N_times_mean_local, MPI.DOUBLE], [N_times_mean_pool, MPI.DOUBLE], op = MPI.SUM, root = 0) # Sum the arrays across processors.
        comm.Reduce([N_local, MPI.DOUBLE],[N_pool, MPI.DOUBLE], op = MPI.SUM, root = 0)   
        comm.Reduce([N_times_var_local, MPI.DOUBLE], [N_times_var_pool, MPI.DOUBLE], op = MPI.SUM, root = 0)
        
    else:
    
        N_times_mean_local = N_local * mean_local
        N_times_var_local = N_local * std_local * std_local

        N_times_mean_pool = comm.reduce(N_times_mean_local, op = MPI.SUM, root = 0)
        N_pool = comm.reduce(N_local, op = MPI.SUM, root = 0)
        N_times_var_pool = comm.reduce(N_times_var_local, op = MPI.SUM, root = 0)
    
    if rank == 0:

        mean_pool_function = np.zeros((len(N_pool)))
        std_pool_function = np.zeros((len(N_pool)))

        for i in range(0, len(N_pool)):
            if N_pool[i] == 0:
                mean_pool_function[i] = 0.0
            else:
                mean_pool_function[i] = np.divide(N_times_mean_pool[i], N_pool[i])
            if N_pool[i] < 3:
                std_pool_function[i] = 0.0
            else:
                std_pool_function[i] = np.sqrt(np.divide(N_times_var_pool[i], N_pool[i]))
       
        mean_pool.append(mean_pool_function)
        std_pool.append(std_pool_function)

        return mean_pool, std_pool 
    else:
    
        return mean_pool, std_pool 
##


def StellarMassFunction(SnapList, SMF, simulation_norm, FirstFile, LastFile, NumFile, ResolutionLimit_mean, model_tags, observations, output_tag):
    '''
    Calculates the stellar mass function for given galaxies with the option to overplot observations by Song et al. (2013) at z = 6, 7, 8 and/or Baldry et al. (2008) at z = 0.1. 
    Parallel compatible.
    NOTE: The plotting assumes the redshifts we are plotting at are (roughly) the same for each model. 

    Parameters
    ---------
    SnapList : Nested 'array-like`, SnapList[model_number0] = [snapshot0_model0, ..., snapshotN_model0], with length equal to the number of models.
        Snapshots that we plot the stellar mass function at for each model.
    SMF : Nested 2-dimensional array, SMF[model_number0][snapshot0]  = [bin0galaxies, ..., binNgalaxies], with length equal to the number of bins (NB_gal). 
        The count of galaxies within each stellar mass bin.  Bounds are given by 'm_gal_low' and 'm_gal_high' in bins given by 'bin_width'. 
    simulation_norm : array with length equal to the number of models.
        Denotes which simulation each model uses.  
        0 : MySim
        1 : Mini-Millennium
        2 : Tiamat (down to z = 5)
        3 : Extended Tiamat (down to z = 1.6ish).
        4 : Britton's Simulation
        5 : Kali
    FirstFile, LastFile, NumFile : array of integers with length equal to the number of models.
        The file numbers for each model that were read in (defined by the range between [FirstFile, LastFile] inclusive) and the TOTAL number of files for this model (we may only be plotting a subset of the volume). 
    ResolutionLimit_mean : array of floats with the same shape as SMF.
        This is the mean stellar mass for a halo with len (number of N-body simulation particles) between 'stellar_mass_halolen_lower' and 'stellar_mass_halolen_upper'. 
    model_tags : array of strings with length equal to the number of models.
        Strings that contain the tag for each model.  Will be placed on the plot.
    observations : int
        Denotes whether we want to overplot observational results. 
        0 : Don't plot anything. 
        1 : Plot Song et al. (2016) at z = 6, 7, 8. 
        2 : Plot Baldry et al. (2008) at z = 0.1.
        3 : Plot both of these.
    output_tag : string
        Name of the file that will be generated. File will be saved in the current directory with the output format defined by the 'output_format' variable at the beggining of the file.

    Returns
    -------
    No returns.
    Generates and saves the plot (named via output_tag).

    Units
    -----
    Stellar Mass is in units of log10(Msun).
    '''

    ## Empty array initialization ##
    title = []
    normalization_array = []
    redshift_labels = []
    counts_array = []
    bin_middle_array = []

    for model_number in range(0, len(SnapList)):
        counts_array.append([])
        bin_middle_array.append([])
        redshift_labels.append([])

    ####
    for model_number in range(0, len(SnapList)): # Does this for each of the models. 

        ## Normalization for each model. ##
        if (simulation_norm[model_number] == 0):
            AllVars.Set_Params_Mysim()
        elif (simulation_norm[model_number] == 1):
            AllVars.Set_Params_MiniMill()
        elif (simulation_norm[model_number] == 2):
            AllVars.Set_Params_Tiamat()
        elif (simulation_norm[model_number] == 3):
            AllVars.Set_Params_Tiamat_extended()
        elif (simulation_norm[model_number] == 4):
            AllVars.Set_Params_Britton()       
        elif(simulation_norm[model_number] == 5):
            AllVars.Set_Params_Kali()
 
        box_factor = (LastFile[model_number] - FirstFile[model_number] + 1.0)/(NumFile[model_number]) # This factor allows us to take a sub-volume of the box and scale the results to represent the entire box.
        print("We are creating the stellar mass function using {0:.4f} of the box's volume.".format(box_factor))
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * bin_width * box_factor 
        normalization_array.append(norm)

        ####
    
        for snapshot_idx in range(0, len(SnapList[model_number])): # Loops for each snapshot in each model.
            tmp = 'z = %.2f' %(AllVars.SnapZ[SnapList[model_number][snapshot_idx]]) # Assigns a redshift label. 
            redshift_labels[model_number].append(tmp)

            ## We perform the plotting on Rank 0 so only this rank requires the final counts array. ##
            if rank == 0:
                counts_total = np.zeros_like(SMF[model_number][snapshot_idx])
            else:
                counts_total = None

            comm.Reduce([SMF[model_number][snapshot_idx], MPI.FLOAT], [counts_total, MPI.FLOAT], op = MPI.SUM, root = 0) # Sum all the stellar mass and pass to Rank 0.

            if rank == 0:
                counts_array[model_number].append(counts_total)
                bin_middle_array[model_number].append(np.arange(m_gal_low, m_gal_high+bin_width, bin_width)[:-1] + bin_width * 0.5)
            ####
    

    ## Plotting ##

    if rank == 0: # Plot only on rank 0.
        f = plt.figure()  
        ax = plt.subplot(111)  

        for model_number in range(0, len(SnapList)):
            for snapshot_idx in range(0, len(SnapList[model_number])):
                if model_number == 0: # We assume the redshifts for each model are the same, we only want to put a legend label for each redshift once.
                    title = redshift_labels[model_number][snapshot_idx]
                else:
                    title = ''
    
                plt.plot(bin_middle_array[model_number][snapshot_idx], counts_array[model_number][snapshot_idx] / normalization_array[model_number], color = PlotScripts.colors[snapshot_idx], linestyle = PlotScripts.linestyles[model_number], rasterized = True, label = title, linewidth = PlotScripts.global_linewidth) 

        print(np.min(np.log10(ResolutionLimit_mean)))
    
        ax.axvline(np.min(np.log10(ResolutionLimit_mean)), color = 'k', linewidth = PlotScripts.global_linewidth, linestyle = '--')    
        ax.text(np.min(np.log10(ResolutionLimit_mean)) + 0.1, 1e-3, "Resolution Limit", color = 'k')
 
        for model_number in range(0, len(SnapList)): # Place legend labels for each of the models. NOTE: Placed after previous loop for proper formatting of labels. 
            plt.plot(1e100, 1e100, color = 'k', linestyle = PlotScripts.linestyles[model_number], label = model_tags[model_number], rasterized=True, linewidth = PlotScripts.global_linewidth)
    
        ## Adjusting axis labels/limits. ##

        plt.yscale('log', nonposy='clip')
        plt.axis([4, 11.5, 1e-6, 1e-0])

        ax.set_xlabel(r'$\log_{10}\ m_{\mathrm{*}} \:[M_{\odot}]$', fontsize = PlotScripts.global_fontsize)
        ax.set_ylabel(r'$\Phi\ [\mathrm{Mpc}^{-3}\: \mathrm{dex}^{-1}]$', fontsize = PlotScripts.global_fontsize)
        ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))
        delta = 0.05
        caps = 5
        ####

        if ((observations == 1 or observations == 3) and rank == 0): # If we wanted to plot Song.
        ## Observational data definition. Columns are redshift, mean value, lower error, upper error.##
            Gonzalez_z6 = np.array([[7.77 + delta, -2.0956, -1.8596, -2.3539],
                        [8.27 + delta, -2.1742, -1.9494, -2.4101],
                        [8.77 + delta, -2.5674, -2.3876, -2.7921],
                        [9.27 + delta, -2.8483, -2.6573, -3.0843],
                        [9.77 + delta, -3.5787, -3.3764, -3.8258],
                        [10.27 + delta, -4.3202, -4.0281, -4.5674]], dtype = np.float32)

            Gonzalez_z7 = np.array([[7.75, -2.1828, -1.7463, -2.5858],
                        [8.26, -2.25, -1.8694, -2.2631],
                        [8.77, -2.7425, -2.3731, -3.1231],
                        [9.27, -3.0672, -2.6753, -3.4142],
                        [9.76, -3.8731, -3.4831, -4.2537]], dtype = np.float32)

            Song_z6 = np.array([[7.25 - delta, -1.47, -1.47 + 0.35, -1.47 - 0.23],
                        [7.75 - delta, -1.81, -1.81 + 0.23, -1.81 - 0.28],
                        [8.25 - delta, -2.26, -2.26 + 0.21, -2.26 - 0.16],
                        [8.75 - delta, -2.65, -2.65 + 0.15, -2.65 - 0.15],
                        [9.25 - delta, -3.14, -3.14 + 0.12, -3.14 - 0.11],
                        [9.75 - delta, -3.69, -3.69 + 0.12, -3.69 - 0.13],
                        [10.25 - delta, -4.27, -4.27 + 0.38, -4.27 - 0.86]], dtype = np.float32)

            Song_z7 = np.array([[7.25, -1.63, -1.63 + 0.54, -1.63 - 0.54],
                        [7.75, -2.07, -2.07 + 0.45, -2.07 - 0.41],
                        [8.25, -2.49, -2.49 + 0.38, -2.49 - 0.32],
                        [8.75, -2.96, -2.96 + 0.32, -2.96 - 0.30],
                        [9.25, -3.47, -3.47 + 0.32, -3.47 - 0.35],
                        [9.75, -4.11, -4.11 + 0.41, -4.11 - 0.57],
                        [10.25, -4.61, -4.61 + 0.72, -4.61 - 0.82],
                        [10.75, -5.24, -5.24 + 0.90, -5.25 - 0.57]], dtype = np.float32)

            Song_z8 = np.array([[7.25, -1.73, -1.73 + 1.01, -1.73 - 0.84],
                        [7.75, -2.28, -2.28 + 0.84, -2.28 - 0.64],
                        [8.25, -2.88, -2.88 + 0.75, -2.88 - 0.57],
                        [8.75, -3.45, -3.45 + 0.57, -3.45 - 0.60],
                        [9.25, -4.21, -4.21 + 0.63, -4.21 - 0.78],
                        [9.75, -5.31, -5.31 + 1.01, -5.31 - 1.64]], dtype = np.float32)
            ####        

            ## Gonzalaez (2011) Plotting ##                   
            #plt.errorbar(Gonzalez_z6[:,0], 10**Gonzalez_z6[:,1], yerr= (10**Gonzalez_z6[:,3], 10**Gonzalez_z6[:,2]) , alpha=0.75, lw=1.0, marker='o', ls='none', label = 'Gonzalez 2011, z = 6', color = 'cyan')
            #plt.errorbar(Gonzalez_z7[:,0], 10**Gonzalez_z7[:,1], yerr= (10**Gonzalez_z7[:,3], 10**Gonzalez_z7[:,2]) , alpha=0.75, lw=1.0, marker='o', ls='none'    , label = 'Gonzalez 2011, z = 7', color = 'magenta')
            ####

            ## Song (2016) Plotting ##
            plt.errorbar(Song_z6[:,0], 10**Song_z6[:,1], yerr= (10**Song_z6[:,1] - 10**Song_z6[:,3], 10**Song_z6[:,2] - 10**Song_z6[:,1]), xerr = 0.25, capsize = caps, elinewidth = PlotScripts.global_errorwidth, alpha = 1.0, lw=2.0, marker='o', ls='none', label = 'Song 2015, z = 6', color = 'green', rasterized=True)
            plt.errorbar(Song_z7[:,0], 10**Song_z7[:,1], yerr= (10**Song_z7[:,1] - 10**Song_z7[:,3], 10**Song_z7[:,2] - 10**Song_z7[:,1]), xerr = 0.25, capsize = caps, alpha=0.75, elinewidth = PlotScripts.global_errorwidth, lw=1.0, marker='o', ls='none', label = 'Song 2015, z = 7', color = 'blue', rasterized=True)
            plt.errorbar(Song_z8[:,0], 10**Song_z8[:,1], yerr= (10**Song_z8[:,1] - 10**Song_z8[:,3], 10**Song_z8[:,2] - 10**Song_z8[:,1]), xerr = 0.25, capsize = caps, alpha=0.75, elinewidth = PlotScripts.global_errorwidth, lw=1.0, marker='o', ls='none', label = 'Song 2015, z = 8', color = '#bd0026', rasterized=True)
            ####

        if ((observations == 2 or observations == 3) and rank == 0): # If we wanted to plot Baldry.
            ## Data definition. Columns are redshift, lower bound, upper bound.##   
            Baldry = np.array([
                [7.05, 1.3531e-01, 6.0741e-02],
                [7.15, 1.3474e-01, 6.0109e-02],
                [7.25, 2.0971e-01, 7.7965e-02],
                [7.35, 1.7161e-01, 3.1841e-02],
                [7.45, 2.1648e-01, 5.7832e-02],
                [7.55, 2.1645e-01, 3.9988e-02],
                [7.65, 2.0837e-01, 4.8713e-02],
                [7.75, 2.0402e-01, 7.0061e-02],
                [7.85, 1.5536e-01, 3.9182e-02],
                [7.95, 1.5232e-01, 2.6824e-02],
                [8.05, 1.5067e-01, 4.8824e-02],
                [8.15, 1.3032e-01, 2.1892e-02],
                [8.25, 1.2545e-01, 3.5526e-02],
                [8.35, 9.8472e-02, 2.7181e-02],
                [8.45, 8.7194e-02, 2.8345e-02],
                [8.55, 7.0758e-02, 2.0808e-02],
                [8.65, 5.8190e-02, 1.3359e-02],
                [8.75, 5.6057e-02, 1.3512e-02],
                [8.85, 5.1380e-02, 1.2815e-02],
                [8.95, 4.4206e-02, 9.6866e-03],
                [9.05, 4.1149e-02, 1.0169e-02],
                [9.15, 3.4959e-02, 6.7898e-03],
                [9.25, 3.3111e-02, 8.3704e-03],
                [9.35, 3.0138e-02, 4.7741e-03],
                [9.45, 2.6692e-02, 5.5029e-03],
                [9.55, 2.4656e-02, 4.4359e-03],
                [9.65, 2.2885e-02, 3.7915e-03],
                [9.75, 2.1849e-02, 3.9812e-03],
                [9.85, 2.0383e-02, 3.2930e-03],
                [9.95, 1.9929e-02, 2.9370e-03],
                [10.05, 1.8865e-02, 2.4624e-03],
                [10.15, 1.8136e-02, 2.5208e-03],
                [10.25, 1.7657e-02, 2.4217e-03],
                [10.35, 1.6616e-02, 2.2784e-03],
                [10.45, 1.6114e-02, 2.1783e-03],
                [10.55, 1.4366e-02, 1.8819e-03],
                [10.65, 1.2588e-02, 1.8249e-03],
                [10.75, 1.1372e-02, 1.4436e-03],
                [10.85, 9.1213e-03, 1.5816e-03],
                [10.95, 6.1125e-03, 9.6735e-04],
                [11.05, 4.3923e-03, 9.6254e-04],
                [11.15, 2.5463e-03, 5.0038e-04],
                [11.25, 1.4298e-03, 4.2816e-04],
                [11.35, 6.4867e-04, 1.6439e-04],
                [11.45, 2.8294e-04, 9.9799e-05],
                [11.55, 1.0617e-04, 4.9085e-05],
                [11.65, 3.2702e-05, 2.4546e-05],
                [11.75, 1.2571e-05, 1.2571e-05],
                [11.85, 8.4589e-06, 8.4589e-06],
                [11.95, 7.4764e-06, 7.4764e-06],
                ], dtype=np.float32)
            

            Baldry_xval = np.log10(10 ** Baldry[:, 0]  /AllVars.Hubble_h/AllVars.Hubble_h)
            Baldry_xval = Baldry_xval - 0.26  # convert back to Chabrier IMF
            Baldry_yvalU = (Baldry[:, 1]+Baldry[:, 2]) * AllVars.Hubble_h*AllVars.Hubble_h*AllVars.Hubble_h
            Baldry_yvalL = (Baldry[:, 1]-Baldry[:, 2]) * AllVars.Hubble_h*AllVars.Hubble_h*AllVars.Hubble_h

            plt.fill_between(Baldry_xval, Baldry_yvalU, Baldry_yvalL, 
                facecolor='purple', alpha=0.25, label='Baldry et al. 2008 (z=0.1)')

            ####

        leg = plt.legend(loc='lower left', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize(PlotScripts.global_legendsize)

        outputFile = './%s%s' %(output_tag, output_format) 
        plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
        print('Saved file to {0}'.format(outputFile))
        plt.close()

##

def plot_fesc(SnapList, mean_z_fesc, std_z_fesc, N_fesc, model_tags, output_tag):
    '''
    Plots the escape fraction as a function of redshift for the given galaxies. 
    Parallel compatible.
    Accepts 2D arrays of the escape fraction at each redshift for each model. 

    Parameters
    ---------
    SnapList : Nested array, SnapList[model_number0] = [snapshot0_model0, ..., snapshotN_model0], with length equal to the number of models.
        Snapshots for each model. 
    mean_z_fesc, std_z_fesc, N_fesc : Nested 2-dimensional array, mean_z_fesc[model_number0] = [z0_meanfesc, ..., zN_meanfesc], with length equal to the number of models 
        Mean/Standard deviation for fesc at each redshift. N_fesc is the number of data points in each bin. 
    model_tags : array of strings with length equal to the number of models.
        Strings that contain the tag for each model.  Will be placed on the plot.
    output_tag : string
        Name of the file that will be generated.

    Returns
    -------
    No returns.
    Generates and saves the plot (named via output_tag).   
    '''

    print("Plotting fesc as a function of redshift.")

    ## Array initialization ##
    pooled_mean_fesc = []
    pooled_std_fesc = []

    for model_number in range(0, len(SnapList)): # Loop for each model. 
    
        pooled_mean_fesc, pooled_std_fesc = calculate_pooled_stats(pooled_mean_fesc, pooled_std_fesc, mean_z_fesc[model_number], std_z_fesc[model_number], N_fesc[model_number]) # Calculates the pooled mean/standard deviation for this snapshot.  Only rank 0 receives a proper value here; the other ranks don't need this information. 
    
    if (rank == 0):
        ax1 = plt.subplot(111)

        for model_number in range(0, len(SnapList)):
    
            ## Calculate lookback time for each snapshot ##
            t = np.empty(len(SnapList[model_number]))
            for snapshot_idx in range(0, len(SnapList[model_number])):  
                t[snapshot_idx] = (t_BigBang - cosmo.lookback_time(AllVars.SnapZ[SnapList[model_number][snapshot_idx]]).value) * 1.0e3   
                    
            mean = pooled_mean_fesc[model_number]
            std = pooled_std_fesc[model_number]   

            ax1.plot(t, mean, color = PlotScripts.colors[model_number], linestyle = PlotScripts.linestyles[model_number], label = model_tags[model_number], linewidth = PlotScripts.global_linewidth)  
            ax1.fill_between(t, np.subtract(mean,std), np.add(mean,std), color = PlotScripts.colors[model_number], alpha = 0.25)

        ax1.xaxis.set_minor_locator(mtick.MultipleLocator(PlotScripts.time_tickinterval))
        ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.025))
        ax1.set_xlim(PlotScripts.time_xlim)

        ## Create a second axis at the top that contains the corresponding redshifts. ##
        ## The redshift defined in the variable 'z_plot' will be displayed. ##
        ax2 = ax1.twiny()

        t_plot = (t_BigBang - cosmo.lookback_time(PlotScripts.z_plot).value) * 1.0e3 # Corresponding time values on the bottom.
        z_labels = ["$%d$" % x for x in PlotScripts.z_plot] # Properly Latex-ize the labels.

        ax2.set_xlabel(r"$z$", size = PlotScripts.global_labelsize) 
        ax2.set_xlim(PlotScripts.time_xlim)
        ax2.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
        ax2.set_xticklabels(z_labels) # But label them as redshifts.

        ax1.set_ylim([0.0, 1.0])
        ax1.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = PlotScripts.global_labelsize) 
        ax1.set_ylabel(r'$f_\mathrm{esc}$', fontsize = PlotScripts.global_fontsize) 

        leg = ax1.legend(loc=1, numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize(PlotScripts.global_legendsize)

        plt.tight_layout()
        outputFile = './{0}{1}'.format(output_tag, output_format)
        plt.savefig(outputFile)  # Save the figure
        print('Saved file to {0}'.format(outputFile))
        plt.close()

##

def plot_fesc_galaxy(SnapList, PlotSnapList, simulation_norm, mean_galaxy_fesc, std_galaxy_fesc, N_galaxy_fesc, mean_halo_fesc, std_halo_fesc, N_halo_fesc, ResolutionLimit_mean, model_tags, output_tag, save_data, save_tags):
    """
    Plots the escape fraction as a function of stellar/halo mass.
    Parallel compatible.
    Accepts 3D arrays of the escape fraction binned into Stellar Mass bins to plot the escape fraction for multiple models. 
    Mass units are 1e10 Msun (no h).

    Parameters
    ---------
    SnapList : Nested array, SnapList[model_number0] = [snapshot0_model0, ..., snapshotN_model0], with length equal to the number of models.
        Snapshots for each model. 
    simulation_norm : array with length equal to the number of models.
        Denotes which simulation each model uses.  
        0 : MySim
        1 : Mini-Millennium
        2 : Tiamat (down to z = 5)
        3 : Extended Tiamat (down to z = 1.6ish).
        4 : Britton's Simulation
        5 : Kali
    mean_galaxy_fesc, std_galaxy_fesc, N_galaxy_fesc : Nested 3-dimensional array, mean_galaxy_fesc[model_number0][snapshot0]  = [bin0_meanfesc, ..., binN_meanfesc], with length equal to the number of models. 
        Mean/Standard deviation for fesc in each stellar mass bin, for each [model_number] and [snapshot_number]. N_galaxy_fesc is the number of galaxies placed into each mass bin.
    mean_halo_fesc, std_halo_fesc, N_halo_fesc  Nested 3-dimensional array, mean_halo_fesc[model_number0][snapshot0]  = [bin0_meanfesc, ..., binN_meanfesc], with length equal to the number of models. 
        Identical to previous except using the halo virial mass for the binning rather than stellar mass. 
    ResolutionLimit_mean : array of floats with the same shape as mean_galaxy_fesc.
        This is the mean stellar mass for a halo with len (number of N-body simulation particles) between 'stellar_mass_halolen_lower' and 'stellar_mass_halolen_upper'. 
    model_tags : array of strings with length equal to the number of models.
        Strings that contain the tag for each model.  Will be placed on the plot.
    output_tag : string
        Name of the file that will be generated.

    Returns
    -------
    No returns.
    Generates and saves the plot (named via output_tag).  

    Units
    -----

    Mass units are 1e10 Msun (no h). 
    """

    print("Plotting fesc as a function of stellar mass.")

    ## Array initialization ##
    title = []
    redshift_labels = []

    mean_fesc_stellar_array = []
    std_fesc_stellar_array = []

    mean_fesc_halo_array = []
    std_fesc_halo_array = []

    bin_middle_stellar_array = []
    bin_middle_halo_array = []

    for model_number in range(0, len(SnapList)):
        redshift_labels.append([])

        mean_fesc_stellar_array.append([])
        std_fesc_stellar_array.append([])

        mean_fesc_halo_array.append([])
        std_fesc_halo_array.append([])

        bin_middle_stellar_array.append([])
        bin_middle_halo_array.append([])

    for model_number in range(0, len(SnapList)): 

        ## Normalization for each model. ##
        if (simulation_norm[model_number] == 0):
            AllVars.Set_Params_Mysim()
        elif (simulation_norm[model_number] == 1):
            AllVars.Set_Params_MiniMill()
        elif (simulation_norm[model_number] == 2):
            AllVars.Set_Params_Tiamat()
        elif (simulation_norm[model_number] == 3):
            AllVars.Set_Params_Tiamat_extended()
        elif (simulation_norm[model_number] == 4):
            AllVars.Set_Params_Britton()       
        elif(simulation_norm[model_number] == 5):
            AllVars.Set_Params_Kali()


        for snapshot_idx in range(0, len(SnapList[model_number])):
            tmp = 'z = %.2f' %(AllVars.SnapZ[SnapList[model_number][snapshot_idx]])
            redshift_labels[model_number].append(tmp)
 
            mean_fesc_stellar_array[model_number], std_fesc_stellar_array[model_number] = calculate_pooled_stats(mean_fesc_stellar_array[model_number], std_fesc_stellar_array[model_number], mean_galaxy_fesc[model_number][snapshot_idx], std_galaxy_fesc[model_number][snapshot_idx], N_galaxy_fesc[model_number][snapshot_idx]) 

            mean_fesc_halo_array[model_number], std_fesc_halo_array[model_number] = calculate_pooled_stats(mean_fesc_halo_array[model_number], std_fesc_halo_array[model_number], mean_halo_fesc[model_number][snapshot_idx], std_halo_fesc[model_number][snapshot_idx], N_halo_fesc[model_number][snapshot_idx]) 

            bin_middle_stellar_array[model_number].append(np.arange(m_gal_low, m_gal_high+bin_width, bin_width)[:-1] + bin_width * 0.5)
            bin_middle_halo_array[model_number].append(np.arange(m_low, m_high+bin_width, bin_width)[:-1] + bin_width * 0.5)
 
    if rank == 0:
        
        fig = plt.figure()  
        ax1 = fig.add_subplot(111)  
        ax2 = ax1.twinx()

        fig2 = plt.figure()
        ax3 = fig2.add_subplot(111)
        ax4 = ax3.twinx()
       
        for model_number in range(0, len(SnapList)):

            ## Normalization for each model. ##
            if (simulation_norm[model_number] == 0):
                AllVars.Set_Params_Mysim()
            elif (simulation_norm[model_number] == 1):
                AllVars.Set_Params_MiniMill()
            elif (simulation_norm[model_number] == 2):
                AllVars.Set_Params_Tiamat()
            elif (simulation_norm[model_number] == 3):
                AllVars.Set_Params_Tiamat_extended()
            elif (simulation_norm[model_number] == 4):
                AllVars.Set_Params_Britton()       
            elif(simulation_norm[model_number] == 5):
                AllVars.Set_Params_Kali()

            if (save_data == 1):
                fname = "/lustre/projects/p004_swin/jseiler/kali_analysis/halo_bins_{0}.txt".format(save_tags[model_number])
                np.savetxt(fname, bin_middle_halo_array[model_number])

                fname = "/lustre/projects/p004_swin/jseiler/kali_analysis/mean_fesc_halo_{0}.txt".format(save_tags[model_number])
                np.savetxt(fname, mean_fesc_halo_array[model_number])

                fname = "/lustre/projects/p004_swin/jseiler/kali_analysis/N_halo_{0}.txt".format(save_tags[model_number])
                np.savetxt(fname, N_halo_fesc[model_number]) 


            plot_count = 0
            for snapshot_idx in range(0, len(SnapList[model_number])):
                
                if (SnapList[model_number][snapshot_idx] == PlotSnapList[model_number][plot_count]):

                    if (model_number == 0):
                        label = redshift_labels[model_number][snapshot_idx]
                    else:
                        label = ""

                    ## Plots as a function of stellar mass ##
                    w = np.where((N_galaxy_fesc[model_number][snapshot_idx] < 1))[0] # If there are no galaxies in the bin we don't want to plot. 
                    N_galaxy_fesc[model_number][snapshot_idx][w] = np.nan 
                    mean_fesc_stellar_array[model_number][snapshot_idx][w] = np.nan

                    ax1.plot(bin_middle_stellar_array[model_number][snapshot_idx], mean_fesc_stellar_array[model_number][snapshot_idx], color = PlotScripts.colors[plot_count], linestyle = PlotScripts.linestyles[model_number], rasterized = True, label = label, linewidth = PlotScripts.global_linewidth) # Plots the escape fraction on the left.
                    ax2.plot(bin_middle_stellar_array[model_number][snapshot_idx], N_galaxy_fesc[model_number][snapshot_idx], color = PlotScripts.colors[plot_count], linestyle = '-.', rasterized = True, linewidth = PlotScripts.global_linewidth) # And the number of galaxies in the bin on the right.

                    print("Resolution limit for model {0} at snapshot {1} is {2}".format(model_number, snapshot_idx, np.log10(ResolutionLimit_mean[model_number][snapshot_idx])))
                    if plot_count == 2 and model_number == 0:
                        ax1.axvline(np.log10(ResolutionLimit_mean[model_number][snapshot_idx]), color = 'k', linewidth = PlotScripts.global_linewidth, linestyle = '--')    

                    ## Plots as a function of halo mass ##
                    w = np.where((N_halo_fesc[model_number][snapshot_idx] < 1))[0] # If there are no galaxies in the bin we don't want to plot.                     
                    N_halo_fesc[model_number][snapshot_idx][w] = np.nan 
                    mean_fesc_halo_array[model_number][snapshot_idx][w] = np.nan

                    '''
                    if (model_number == 0):
                        print(bin_middle_halo_array[model_number][snapshot_idx])
                        print(mean_fesc_halo_array[model_number][snapshot_idx])
                    '''
                    ax3.plot(bin_middle_halo_array[model_number][snapshot_idx], mean_fesc_halo_array[model_number][snapshot_idx], color = PlotScripts.colors[plot_count], linestyle = PlotScripts.linestyles[model_number], rasterized = True, label = label, linewidth = PlotScripts.global_linewidth) # Plots the escape fraction on the left.
                    ax4.plot(bin_middle_halo_array[model_number][snapshot_idx], N_halo_fesc[model_number][snapshot_idx], color = PlotScripts.colors[plot_count], linestyle = '-.', rasterized = True, linewidth = PlotScripts.global_linewidth) # And the number of halos in the bin on the right.

                    plot_count += 1                
                    if (plot_count == len(PlotSnapList[model_number])):
                        break

        for model_number in range(0, len(SnapList)): # Just plot some garbage to get the legend labels correct.
            ax1.plot(np.nan, np.nan, color = 'k', linestyle = PlotScripts.linestyles[model_number], rasterized = True, label = model_tags[model_number], linewidth = PlotScripts.global_linewidth)
            ax3.plot(np.nan, np.nan, color = 'k', linestyle = PlotScripts.linestyles[model_number], rasterized = True, label = model_tags[model_number], linewidth = PlotScripts.global_linewidth)

        ## Stellar Mass plots ##
 
        ax1.axhline(0.2, 0, 100, color ='k', linewidth = PlotScripts.global_linewidth, linestyle = '--')

        ax1.set_xlabel(r'$\log_{10}\ M_*\ [M_{\odot}]$', size = PlotScripts.global_fontsize) 
        ax1.set_ylabel(r'$\mathrm{f_\mathrm{esc}}$', size = PlotScripts.global_fontsize)
        ax2.set_ylabel(r'$N_\mathrm{gal}$', size = PlotScripts.global_fontsize)
        ax1.set_xlim([4.0, 12])
        ax1.set_ylim([-0.05, 1.0])   

        ax1.xaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.025))
        
        ax2.set_yscale('log', nonposy='clip')

        leg = ax1.legend(loc=1, numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')

        ## Halo mass plots ##

        ax3.axhline(0.2, 0, 100, color ='k', linewidth = PlotScripts.global_linewidth, linestyle = '--')
        ax3.axvline(np.log10(32.0*AllVars.PartMass / AllVars.Hubble_h), color = 'k', linewidth = PlotScripts.global_linewidth, linestyle = '--')   
 
        ax3.set_xlabel(r'$\log_{10}\ M_\mathrm{vir}\ [M_{\odot}]$', size = PlotScripts.global_fontsize) 
        ax3.set_ylabel(r'$\mathrm{f_\mathrm{esc}}$', size = PlotScripts.global_fontsize)
        ax4.set_ylabel(r'$N_\mathrm{Halo}$', size = PlotScripts.global_fontsize)
        ax3.set_xlim([7.0, 14])
        ax3.set_ylim([-0.05, 1.0])   

        ax3.set_xticks(np.arange(7.0, 14.0))  
        ax3.xaxis.set_minor_locator(mtick.MultipleLocator(0.25))
        ax3.yaxis.set_minor_locator(mtick.MultipleLocator(0.025))
        
        ax4.set_yscale('log', nonposy='clip')

        leg = ax3.legend(loc=1, numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')

        ## Output ##

        outputFile = './%s%s' %(output_tag, output_format)
        fig.savefig(outputFile, bbox_inches='tight')  # Save the figure
        print('Saved file to {0}'.format(outputFile))

        outputFile = './%s_Halo%s' %(output_tag, output_format)
        fig2.savefig(outputFile, bbox_inches='tight')  # Save the figure
        print('Saved file to {0}'.format(outputFile))

        plt.close(fig)
        plt.close(fig2)

##

def plot_ejectedfraction(SnapList, mean_mvir_ejected, std_mvir_ejected, N_ejected, model_tags, output_tag): 
    '''
    Plots the ejected fraction as a function of the halo mass. 
    Parallel compatible.
    Accepts a 3D array of the ejected fraction so we can plot for multiple models and redshifts. 

    Parameters
    ---------
    SnapList : Nested array, SnapList[model_number0] = [snapshot0_model0, ..., snapshotN_model0], with length equal to the number of models.
    Snapshots for each model. 
    mean_mvir_ejected, std_mvir_ejected, N_ejected : Nested 3-dimensional array, mean_mvir_ejected[model_number0][snapshot0]  = [bin0_meanejected, ..., binN_meanejected], with length equal to the number of models. 
    Mean/Standard deviation for the escape fraction binned into Halo Mass bins. N_ejected is the number of data points in each bin. Bounds are given by 'm_low' and 'm_high' in bins given by 'bin_width'.   
    model_tags : array of strings with length equal to the number of models.
    Strings that contain the tag for each model.  Will be placed on the plot.
    output_tag : string
    Name of the file that will be generated.

    Returns
    -------
    No returns.
    Generates and saves the plot (named via output_tag).

    Units
    -----
    Halo Mass is in units of log10(Msun). 
    '''

    print("Plotting the Ejected Fraction as a function of halo mass.")

    ## Array initialization. ##
    title = []
    redshift_labels = []

    mean_ejected_array = []
    std_ejected_array = []

    mean_halomass_array = []
    std_halomass_array = []

    bin_middle_array = []

    for model_number in range(0, len(SnapList)):
        redshift_labels.append([])

        mean_ejected_array.append([])
        std_ejected_array.append([])

        mean_halomass_array.append([])
        std_halomass_array.append([])

        bin_middle_array.append([])
    
    bin_width = 0.1
 
    for model_number in range(0, len(SnapList)): 
        for snapshot_idx in range(0, len(SnapList[model_number])):
            print("Doing Snapshot {0}".format(SnapList[model_number][snapshot_idx]))
            tmp = 'z = %.2f' %(AllVars.SnapZ[SnapList[model_number][snapshot_idx]])
            redshift_labels[model_number].append(tmp)
            
            mean_ejected_array[model_number], std_ejected_array[model_number] = calculate_pooled_stats(mean_ejected_array[model_number], std_ejected_array[model_number], mean_mvir_ejected[model_number][snapshot_idx], std_mvir_ejected[model_number][snapshot_idx], N_ejected[model_number][snapshot_idx]) # Calculates the pooled mean/standard deviation for this snapshot.  Only rank 0 receives a proper value here; the other ranks don't need this information. 

            bin_middle_array[model_number].append(np.arange(m_low, m_high+bin_width, bin_width)[:-1] + bin_width * 0.5)
    
    if rank == 0:
        f = plt.figure()  
        ax1 = plt.subplot(111)  

        for model_number in range(0, len(SnapList)):
            for snapshot_idx in range(0, len(SnapList[model_number])): 
                ax1.plot(bin_middle_array[model_number][snapshot_idx], mean_ejected_array[model_number][snapshot_idx], color = PlotScripts.colors[snapshot_idx], linestyle = PlotScripts.linestyles[model_number], rasterized = True, label = redshift_labels[model_number][snapshot_idx], linewidth = PlotScripts.global_linewidth) 
                

        for model_number in range(0, len(SnapList)): # Just plot some garbage to get the legend labels correct.
            ax1.plot(np.nan, np.nan, color = 'k', linestyle = PlotScripts.linestyles[model_number], rasterized = True, label = model_tags[model_number], linewidth = PlotScripts.global_linewidth)

        ax1.set_xlabel(r'$\log_{10}\ M_{\mathrm{vir}}\ [M_{\odot}]$', size = PlotScripts.global_fontsize) 
        ax1.set_ylabel(r'$\mathrm{Ejected \: Fraction}$', size = PlotScripts.global_fontsize)
        ax1.set_xlim([8.0, 12])
        ax1.set_ylim([-0.05, 1.0])   

        ax1.xaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.025))
        
        leg = ax1.legend(loc=1, numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')

        outputFile = './%s%s' %(output_tag, output_format)
        plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
        print('Saved file to {0}'.format(outputFile))
        plt.close()

##


def plot_mvir_fesc(SnapList, mass_central, fesc, model_tags, output_tag): 

    title = []
    redshift_labels = []

    mean_fesc_array = []
    std_fesc_array = []

    mean_halomass_array = []
    std_halomass_array = []

    bin_middle_array = []

    for model_number in range(0, len(SnapList)):
        redshift_labels.append([])

    mean_fesc_array.append([])
    std_fesc_array.append([])

    mean_halomass_array.append([])
    std_halomass_array.append([])

    bin_middle_array.append([])
    print("Plotting fesc against Mvir") 
    
    binwidth = 0.1
    Frequency = 1  
 
    for model_number in range(0, len(SnapList)): 
        for snapshot_idx in range(0, len(SnapList[model_number])):
            print("Doing Snapshot {0}".format(SnapList[model_number][snapshot_idx]))
            tmp = 'z = %.2f' %(AllVars.SnapZ[SnapList[model_number][snapshot_idx]])
            redshift_labels[model_number].append(tmp)

            minimum_mass = np.floor(min(mass_central[model_number][snapshot_idx])) - 10*binwidth
            maximum_mass = np.floor(max(mass_central[model_number][snapshot_idx])) + 10*binwidth

            minimum_mass = 6.0
            maximum_mass = 12.0

            binning_minimum = comm.allreduce(minimum_mass, op = MPI.MIN)
            binning_maximum = comm.allreduce(maximum_mass, op = MPI.MAX)
            
            halomass_nonlog = [10**x for x in mass_central[model_number][snapshot_idx]]
            (mean_fesc, std_fesc, N, bin_middle) = AllVars.Calculate_2D_Mean(mass_central[model_number][snapshot_idx], fesc[model_number][snapshot_idx], binwidth, binning_minimum, binning_maximum)

            mean_fesc_array[model_number], std_fesc_array[model_number] = calculate_pooled_stats(mean_fesc_array[model_number], std_fesc_array[model_number], mean_fesc, std_fesc, N)
            mean_halomass_array[model_number], std_halomass_array[model_number] = calculate_pooled_stats(mean_halomass_array[model_number], std_halomass_array[model_number], np.mean(halomass_nonlog), np.std(halomass_nonlog), len(mass_central[model_number][snapshot_idx]))

            ## If want to do mean/etc of halo mass need to update script. ##
            bin_middle_array[model_number].append(bin_middle)
        
        mean_halomass_array[model_number] = np.log10(mean_halomass_array[model_number]) 
        
    if rank == 0:
        f = plt.figure()  
        ax1 = plt.subplot(111)  

        for model_number in range(0, len(SnapList)):
            for snapshot_idx in range(0, len(SnapList[model_number])):
                if model_number == 0:
                    title = redshift_labels[model_number][snapshot_idx]
                else:
                    title = ''
                
                mean = mean_fesc_array[model_number][snapshot_idx]
                std = std_fesc_array[model_number][snapshot_idx] 
                bin_middle = bin_middle_array[model_number][snapshot_idx]

                ax1.plot(bin_middle, mean, color = colors[snapshot_idx], linestyle = linestyles[model_number], rasterized = True, label = title)
                #ax1.scatter(mean_halomass_array[model_number][snapshot_idx], np.mean(~np.isnan(mean)), color = colors[snapshot_idx], marker = 'o', rasterized = True, s = 40, lw = 3)  
                if (len(SnapList) == 1):
                        ax1.fill_between(bin_middle, np.subtract(mean,std), np.add(mean,std), color = colors[snapshot_idx], alpha = 0.25)

        ax1.set_xlabel(r'$\log_{10}\ M_{\mathrm{vir}}\ [M_{\odot}]$', size = PlotScripts.global_fontsize) 
        ax1.set_ylabel(r'$f_\mathrm{esc}$', size = PlotScripts.global_fontsize) 
        #ax1.set_xlim([8.5, 12])
        #ax1.set_ylim([0.0, 1.0])   

        ax1.xaxis.set_minor_locator(mtick.MultipleLocator(0.1))
#       ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
    
#       ax1.set_yscale('log', nonposy='clip')
#       for model_number in range(0, len(SnapList)):
#       ax1.plot(1e100, 1e100, color = 'k', ls = linestyles[model_number], label = model_tags[model_number], rasterized=True)
    
    
        leg = ax1.legend(loc='upper left', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')

        outputFile = './' + output_tag + output_format
        plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
        print('Saved file to'.format(outputFile))
        plt.close()
##

def plot_mvir_Ngamma(SnapList, mean_mvir_Ngamma, std_mvir_Ngamma, N_Ngamma, model_tags, output_tag,fesc_prescription=None, fesc_normalization=None, fitpath=None): 
    '''
    Plots the number of ionizing photons (pure ngamma times fesc) as a function of halo mass. 
    Parallel compatible.
    The input data has been binned as a function of halo virial mass (Mvir), with the bins defined at the top of the file (m_low, m_high, bin_width). 
    Accepts 3D arrays to plot ngamma for multiple models. 

    Parameters
    ----------
    SnapList : Nested array, SnapList[model_number0] = [snapshot0_model0, ..., snapshotN_model0], with length equal to the number of models.
    Snapshots for each model. 
    mean_mvir_Ngamma, std_mvir_Ngamma, N_Ngamma : Nested 2-dimensional array, mean_mvir_Ngamma[model_number0][snapshot0]  = [bin0_meanNgamma, ..., binN_meanNgamma], with length equal to the number of bins. 
    Mean/Standard deviation/number of data points in each halo mass (Mvir) bin.
    The number of photons is in units of 1.0e50 s^-1.   
    model_tags : array of strings with length equal to the number of models.
    Strings that contain the tag for each model.  Will be placed on the plot.
    output_tag : string
    Name of the file that will be generated.
    fesc_prescription : int (optional)
    If this parameter is defined, we will save the Mvir-Ngamma results in a text file (not needed if not saving).
    Number that controls what escape fraction prescription was used to generate the escape fractions.
    0 : Constant, fesc = Constant.
    1 : Scaling with Halo Mass, fesc = A*Mh^B.
    2 : Scaling with ejected fraction, fesc = fej*A + B.
    fesc_normalization : float (if fesc_prescription == 0) or `numpy.darray' with length 2 (if fesc_prescription == 1 or == 2) (optional).
    If this parameter is defined, we will save the Mvir-Ngamma results in a text file (not needed if not saving).
    Parameter not needed if you're not saving the Mvir-Ngamma results.
    If fesc_prescription == 0, gives the constant value for the escape fraction.
    If fesc_prescription == 1 or == 2, gives A and B with the form [A, B].
    fitpath : string (optional) 
    If this parameter is defined, we will save the Mvir-Ngamma results in a text file (not needed if not saving).
    Defines the base path for where we are saving the results.

    Returns
    -------
    No returns.
    Generates and saves the plot (named via output_tag).

    Units
    -----
    Ngamma is in units of 1.0e50 s^-1. 
    '''

    print("Plotting ngamma*fesc against the halo mass") 

    ## Array initialization. ##
    title = []
    redshift_labels = []

    mean_ngammafesc_array = []
    std_ngammafesc_array = []

    mean_halomass_array = []
    std_halomass_array = []

    bin_middle_array = []

    for model_number in range(0, len(SnapList)):
        redshift_labels.append([])

    mean_ngammafesc_array.append([])
    std_ngammafesc_array.append([])

    mean_halomass_array.append([])
    std_halomass_array.append([])

    bin_middle_array.append([])
    
    for model_number in range(0, len(SnapList)): 
        for snapshot_idx in range(0, len(SnapList[model_number])):
            print("Doing Snapshot {0}".format(SnapList[model_number][snapshot_idx]))
            tmp = 'z = %.2f' %(AllVars.SnapZ[SnapList[model_number][snapshot_idx]])
            redshift_labels[model_number].append(tmp)

            N = N_Ngamma[model_number][snapshot_idx]
            
            mean_ngammafesc_array[model_number], std_ngammafesc_array[model_number] = calculate_pooled_stats(mean_ngammafesc_array[model_number], std_ngammafesc_array[model_number], mean_mvir_Ngamma[model_number][snapshot_idx], std_mvir_Ngamma[model_number][snapshot_idx], N) # Collate the values from all processors.   
            bin_middle_array[model_number].append(np.arange(m_low, m_high+bin_width, bin_width)[:-1] + bin_width * 0.5) 
    
    if rank == 0:
        f = plt.figure()  
        ax1 = plt.subplot(111)  

        for model_number in range(0, len(SnapList)):
            count = 0
            for snapshot_idx in range(0, len(SnapList[model_number])):
                if model_number == 0:
                    title = redshift_labels[model_number][snapshot_idx]
                else:
                    title = ''

                mean = np.zeros((len(mean_ngammafesc_array[model_number][snapshot_idx])), dtype = np.float32)
                std = np.zeros((len(mean_ngammafesc_array[model_number][snapshot_idx])), dtype=np.float32)

                for i in range(0, len(mean)):
                    if(mean_ngammafesc_array[model_number][snapshot_idx][i] < 1e-10):
                        mean[i] = np.nan
                        std[i] = np.nan
                    else:   
                        mean[i] = np.log10(mean_ngammafesc_array[model_number][snapshot_idx][i] * 1.0e50) # Remember that the input data is in units of 1.0e50 s^-1.
                        std[i] = 0.434 * std_ngammafesc_array[model_number][snapshot_idx][i] / mean_ngammafesc_array[model_number][snapshot_idx][i] # We're plotting in log space so the standard deviation is 0.434*log10(std)/log10(mean).

                bin_middle = bin_middle_array[model_number][snapshot_idx]
    
                if (count < 4): # Only plot at most 5 lines.
                    ax1.plot(bin_middle, mean, color = PlotScripts.colors[snapshot_idx], linestyle = PlotScripts.linestyles[model_number], rasterized = True, label = title, linewidth = PlotScripts.global_linewidth)  
                    count += 1
                ## In this block we save the Mvir-Ngamma results to a file. ##
                if (fesc_prescription == None or fesc_normalization == None or fitpath == None):
                    raise ValueError("You've specified you want to save the Mvir-Ngamma results but haven't provided an escape fraction prescription, normalization and base path name")
                # Note: All the checks that escape fraction normalization was written correctly were performed in 'calculate_fesc()', hence it will be correct by this point and we don't need to double check. 

                if (fesc_prescription[model_number] == 0): # Slightly different naming scheme for the constant case (it only has a float for fesc_normalization).
                    fname = "%s/fesc%d_%.3f_z%.3f.txt" %(fitpath, fesc_prescription[model_number], fesc_normalization[model_number], AllVars.SnapZ[SnapList[model_number][snapshot_idx]])
                elif (fesc_prescription[model_number] == 1 or fesc_prescription[model_number] == 2):    
                    fname = "%s/fesc%d_A%.3eB%.3f_z%.3f.txt" %(fitpath, fesc_prescription[model_number], fesc_normalization[model_number][0], fesc_normalization[model_number][1], AllVars.SnapZ[SnapList[model_number][snapshot_idx]])
                f = open(fname, "w+")
                if not os.access(fname, os.W_OK):
                    print("The filename is {0}".format(fname))
                    raise ValueError("Can't write to this file.")
        
                for i in range(0, len(bin_middle)):
                    f.write("%.4f %.4f %.4f %d\n" %(bin_middle[i], mean[i], std[i], N_Ngamma[model_number][snapshot_idx][i]))
                f.close() 
                print("Wrote successfully to file {0}".format(fname))
            ##

        for model_number in range(0, len(SnapList)): # Just plot some garbage to get the legend labels correct.
            ax1.plot(np.nan, np.nan, color = 'k', linestyle = PlotScripts.linestyles[model_number], rasterized = True, label = model_tags[model_number], linewidth = PlotScripts.global_linewidth)
    
        ax1.set_xlabel(r'$\log_{10}\ M_{\mathrm{vir}}\ [M_{\odot}]$', size = PlotScripts.global_fontsize) 
        ax1.set_ylabel(r'$\log_{10}\ \dot{N}_\gamma \: f_\mathrm{esc} \: [\mathrm{s}^{-1}]$', size = PlotScripts.global_fontsize) 
        ax1.set_xlim([8.5, 12])

        ax1.xaxis.set_minor_locator(mtick.MultipleLocator(0.1)) 
    
        leg = ax1.legend(loc='upper left', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')

        outputFile = './' + output_tag + output_format
        plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
        print('Saved file to'.format(outputFile))
    
        plt.close()


def bin_Simfast_halos(RedshiftList, SnapList, halopath, fitpath, fesc_prescription, fesc_normalization, GridSize, output_tag):
   
    for model_number in range(0, len(fesc_prescription)):
        for halo_z_idx in range(0, len(RedshiftList)):
            snapshot_idx = min(range(len(SnapList)), key=lambda i: abs(SnapList[i]-RedshiftList[halo_z_idx])) # This finds the index of the simulation redshift that most closely matches the Halo redshift.
            print("Binning Halo redshift {0}".format(RedshiftList[halo_z_idx]))
            print("For the Halo redshift {0:.3f} the nearest simulation redshift is {1:.3f}".format(RedshiftList[halo_z_idx], SnapList[snapshot_idx])) 
            if (fesc_prescription[model_number] == 0):
                fname = "%s/fesc%d_%.3f_z%.3f.txt" %(fitpath, fesc_prescription[model_number], fesc_normalization[model_number], AllVars.SnapZ[snapshot_idx]) 
            elif (fesc_prescription[model_number] == 1 or fesc_prescription[model_number] == 2):
                fname = "%s/fesc%d_A%.3eB%.3f_z%.3f.txt" %(fitpath, fesc_prescription[model_number], fesc_normalization[model_number][0], fesc_normalization[model_number][1], AllVars.SnapZ[snapshot_idx])

            print("Reading in file {0}".format(fname))
            ## Here we read in the results from the Mvir-Ngamma binning. ##
            f = open(fname, 'r')
            fit_mvir, fit_mean, fit_std, fit_N = np.loadtxt(f, unpack = True)
            f.close()

            ## Here we read in the halos created by Simfast21 ##
            # The data file has the structure:
                # long int N_halos
                # Then an entry for each halo:
                    # float Mass
                    # float x, y, z positions.
                # NOTE: The x,y,z positions are the grid indices but are still floats (because Simfast21 is weird like that).

            Halodesc_full = [ 
                 ('Halo_Mass', np.float32),
                 ('Halo_x', np.float32),
                 ('Halo_y', np.float32),
                 ('Halo_z', np.float32)
                ]

            names = [Halodesc_full[i][0] for i in range(len(Halodesc_full))]
            formats = [Halodesc_full[i][1] for i in range(len(Halodesc_full))]
            Halo_Desc = np.dtype({'names':names, 'formats':formats}, align=True)

            fname = "%s/halonl_z%.3f_N%d_L100.0.dat.catalog" %(halopath, RedshiftList[halo_z_idx], GridSize)
            f = open(fname, 'rb')
            N_Halos = np.fromfile(f, count = 1, dtype = np.long)    
            Halos = np.fromfile(f, count = N_Halos, dtype = Halo_Desc)  

            binned_nion = np.zeros((GridSize*GridSize*GridSize), dtype = float32) # This grid will contain the ionizing photons that results from the binning.
            binned_Halo_Mass = np.digitize(np.log10(Halos['Halo_Mass']), fit_mvir) # Places the Simfast21 halos into the correct halo mass bins defined by the Mvir-Ngamma results.
            binned_Halo_Mass[binned_Halo_Mass == len(fit_mvir)] = len(fit_mvir) - 1 # Fixes up the edge case.

            ## Fore each Halo we now assign it an ionizing flux. ##
            # This flux is determined by drawing a random number from a normal distribution with mean and standard deviation given by the Mvir-Ngamma results.
            # NOTE: Remember the Mvir-Ngamma results are in units of log10(s^-1).
            fit_nan = 0
            for i in range(0, N_Halos):
                if(np.isnan(fit_mean[binned_Halo_Mass[i]]) == True or np.isnan(fit_std[binned_Halo_Mass[i]]) == True): # This halo had mass that was not covered by the Mvir-Ngamma fits.
                    fit_nan += 1
                    continue
                nion_halo = np.random.normal(fit_mean[binned_Halo_Mass[i]], fit_std[binned_Halo_Mass[i]])

                ## Because of how Simfast21 does their binning, we have some cases where the Halos are technically outside the box. Just fix them up. ##
                x_grid = int(Halos['Halo_x'][i])
                if x_grid >= GridSize:
                    x_grid = GridSize - 1
                if x_grid < 0:
                    x_grid = 0

                y_grid = int(Halos['Halo_y'][i])
                if y_grid >= GridSize:
                    y_grid = GridSize - 1
                if y_grid < 0:
                    y_grid = 0
                    
                z_grid = int(Halos['Halo_z'][i])
                if z_grid >= GridSize:
                    z_grid = GridSize - 1
                if z_grid < 0:
                    z_grid = 0
                
                idx = x_grid * GridSize*GridSize + y_grid * GridSize + z_grid
                binned_nion[idx]  += pow(10, nion_halo)/1.0e50 
#            print"We had %d halos (out of %d, so %.4f fraction) that had halo mass that was not covered by the Mvir-Ngamma results." %(fit_nan, N_Halos, float(fit_nan)/float(N_Halos))
 #           print "There were %d cells with a non-zero ionizing flux." %(len(binned_nion[binned_nion != 0]))

            binned_nion = binned_nion.reshape((GridSize,GridSize,GridSize))
            cut_slice = 0
            cut_width = 512
            nion_slice = binned_nion[:,:, cut_slice:cut_slice+cut_width].mean(axis=-1)*1.0e50

            ax1 = plt.subplot(211)
            
            im = ax1.imshow(np.log10(nion_slice), interpolation='bilinear', origin='low', extent =[0,AllVars.BoxSize,0,AllVars.BoxSize], cmap = 'Purples', vmin = 48, vmax = 53)

            cbar = plt.colorbar(im, ax = ax1)
            cbar.set_label(r'$\mathrm{log}_{10}N_{\gamma} [\mathrm{s}^{-1}]$')

            ax1.set_xlabel(r'$\mathrm{x}  (h^{-1}Mpc)$')
            ax1.set_ylabel(r'$\mathrm{y}  (h^{-1}Mpc)$')

            ax1.set_xlim([0.0, AllVars.BoxSize])
            ax1.set_ylim([0.0, AllVars.BoxSize])

            title = r"$z = %.3f$" %(RedshiftList[halo_z_idx])
            ax1.set_title(title)
    
            ax2 = plt.subplot(212)

            w = np.where((Halos['Halo_z'][:] > cut_slice) & (Halos['Halo_z'][:] <= cut_slice + cut_width))[0]

            x_plot = Halos['Halo_x'] * float(AllVars.BoxSize)/float(GridSize)
            y_plot = Halos['Halo_y'] * float(AllVars.BoxSize)/float(GridSize)
            z_plot = Halos['Halo_z'][w] * float(AllVars.BoxSize)/float(GridSize)
            
            ax2.scatter(x_plot[w], y_plot[w], s = 2, alpha = 0.5)

            ax2.set_xlabel(r'$\mathrm{x}  (h^{-1}Mpc)$')
            ax2.set_ylabel(r'$\mathrm{y}  (h^{-1}Mpc)$')

            ax2.set_xlim([0.0, AllVars.BoxSize])
            ax2.set_ylim([0.0, AllVars.BoxSize])

            tmp = "z%.3f" %(RedshiftList[halo_z_idx])
            
            plt.tight_layout()
            outputFile = './' + output_tag + tmp + output_format
            plt.savefig(outputFile)  # Save the figure
            print('Saved file to {0}'.format(outputFile))
            plt.close()
            
def plot_photoncount(SnapList, sum_nion, simulation_norm, FirstFile, LastFile, NumFiles, model_tags, output_tag): 
    '''
    Plots the ionizing emissivity as a function of redshift. 
    We normalize the emissivity to Mpc^-3 and this function allows the read-in of only a subset of the volume.
    Parallel compatible.

    Parameters
    ---------
    SnapList : Nested array, SnapList[model_number0] = [snapshot0_model0, ..., snapshotN_model0], with length equal to the number of models.
        Snapshots for each model, defines the x-axis we plot against.
    sum_nion : Nested 1-dimensional array, sum_nion[z0, z1, ..., zn], with length equal to the number of redshifts. 
        Number of escape ionizing photons (i.e., photon rate times the local escape fraction) at each redshift.
        In units of 1.0e50 s^-1.
    simulation_norm : array  of ints with length equal to the number of models.
        Denotes which simulation each model uses.  
        0 : MySim
        1 : Mini-Millennium
        2 : Tiamat (down to z = 5)
        3 : Extended Tiamat (down to z = 1.6ish).
        4 : Britton's Simulation
    FirstFile, LastFile, NumFile : array of integers with length equal to the number of models.
        The file numbers for each model that were read in (defined by the range between [FirstFile, LastFile] inclusive) and the TOTAL number of files for this model (we may only be plotting a subset of the volume). 
    model_tags : array of strings with length equal to the number of models.
        Strings that contain the tag for each model.  Will be placed on the plot.
    output_tag : string
        Name of the file that will be generated.

    Returns
    -------
    No returns.
    Generates and saves the plot (named via output_tag).

    Units
    -----
    sum_nion is in units of 1.0e50 s^-1.
    '''
    print("Plotting the ionizing emissivity.")

    sum_array = []

    for model_number in range(0, len(SnapList)):
        if(simulation_norm[model_number] == 0):
            AllVars.Set_Params_Mysim()
        elif(simulation_norm[model_number] == 3):
            AllVars.Set_Params_Tiamat_extended()
        elif(simulation_norm[model_number] == 4):
            AllVars.Set_Params_Britton()
        elif(simulation_norm[model_number] == 5):
            AllVars.Set_Params_Kali()
        else: 
            print("Simulation norm was set to {0}.".format(simulation_norm[model_number]))
            raise ValueError("This option has been implemented yet.  Get your head in the game Jacob!")
 
        sum_array.append([])
        for snapshot_idx in range(0, len(SnapList[model_number])):
       
            nion_sum_snapshot = comm.reduce(sum_nion[model_number][snapshot_idx], op = MPI.SUM, root = 0)
            if rank == 0:
                sum_array[model_number].append(nion_sum_snapshot * 1.0e50 / (pow(AllVars.BoxSize / AllVars.Hubble_h,3) * (float(LastFile[model_number] - FirstFile[model_number] + 1) / float(NumFiles[model_number]))))
            
    if (rank == 0):
        ax1 = plt.subplot(111)

        for model_number in range(0, len(SnapList)):

            if(simulation_norm[model_number] == 0):
                AllVars.Set_Params_Mysim()
            elif(simulation_norm[model_number] == 3):
                AllVars.Set_Params_Tiamat_extended()
            elif(simulation_norm[model_number] == 4):
                AllVars.Set_Params_Britton()
            elif(simulation_norm[model_number] == 5):
                AllVars.Set_Params_Kali()
            else: 
                print("Simulation norm was set to {0}.".format(simulation_norm[model_number]))
                raise ValueError("This option has been implemented yet.  Get your head in the game Jacob!")


            t = np.empty(len(SnapList[model_number]))
            for snapshot_idx in range(0, len(SnapList[model_number])):
                t[snapshot_idx] = (t_BigBang - cosmo.lookback_time(AllVars.SnapZ[SnapList[model_number][snapshot_idx]]).value) * 1.0e3     
            print("Model {0} has photon count {1}".format(model_number, sum_array[model_number]))           
           
            t = [t for t, N in zip(t, sum_array[model_number]) if N > 1.0]
            sum_array[model_number] = [x for x in sum_array[model_number] if x > 1.0]
 
            print(sum_array[model_number]) 
            ax1.plot(t, np.log10(sum_array[model_number]), color = PlotScripts.colors[model_number], linestyle = PlotScripts.linestyles[model_number], label = model_tags[model_number], linewidth = PlotScripts.global_linewidth)  
            #ax1.fill_between(t, np.subtract(mean,std), np.add(mean,std), color = colors[model_number], alpha = 0.25)

        ax1.xaxis.set_minor_locator(mtick.MultipleLocator(PlotScripts.time_tickinterval))
        #ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.025))
        ax1.set_xlim(PlotScripts.time_xlim)
        ax1.set_ylim([48.5, 51.5])

        ax2 = ax1.twiny()

        t_plot = (t_BigBang - cosmo.lookback_time(PlotScripts.z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
        z_labels = ["$%d$" % x for x in PlotScripts.z_plot] # Properly Latex-ize the labels.

        ax2.set_xlabel(r"$z$", size = PlotScripts.global_labelsize)
        ax2.set_xlim(PlotScripts.time_xlim)
        ax2.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
        ax2.set_xticklabels(z_labels) # But label them as redshifts.

        ax1.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = PlotScripts.global_fontsize)
        ax1.set_ylabel(r'$\sum f_\mathrm{esc}\dot{N}_\gamma \: [\mathrm{s}^{-1}\mathrm{Mpc}^{-3}]$', fontsize = PlotScripts.global_fontsize) 

        plot_time = 1
        bouwens_z = np.arange(6,16) # Redshift range for the observations.
        bouwens_t = (AllVars.t_BigBang - cosmo.lookback_time(bouwens_z).value) * 1.0e3 # Corresponding values for what we will plot on the x-axis.

        bouwens_1sigma_lower = [50.81, 50.73, 50.60, 50.41, 50.21, 50.00, 49.80, 49.60, 49.39, 49.18] # 68% Confidence Intervals for the ionizing emissitivity from Bouwens 2015.
        bouwens_1sigma_upper = [51.04, 50.85, 50.71, 50.62, 50.56, 50.49, 50.43, 50.36, 50.29, 50.23]

        bouwens_2sigma_lower = [50.72, 50.69, 50.52, 50.27, 50.01, 49.75, 49.51, 49.24, 48.99, 48.74] # 95% CI. 
        bouwens_2sigma_upper = [51.11, 50.90, 50.74, 50.69, 50.66, 50.64, 50.61, 50.59, 50.57, 50.55]

        if plot_time == 1:
            ax1.fill_between(bouwens_t, bouwens_1sigma_lower, bouwens_1sigma_upper, color = 'k', alpha = 0.2)
            ax1.fill_between(bouwens_t, bouwens_2sigma_lower, bouwens_2sigma_upper, color = 'k', alpha = 0.4, label = r"$\mathrm{Bouwens \: et \: al. \: (2015)}$")
        else:
            ax1.fill_between(bouwens_z, bouwens_1sigma_lower, bouwens_1sigma_upper, color = 'k', alpha = 0.2)
            ax1.fill_between(bouwens_z, bouwens_2sigma_lower, bouwens_2sigma_upper, color = 'k', alpha = 0.4, label = r"$\mathrm{Bouwens \: et \: al. \: (2015)}$")

    #   ax1.text(0.075, 0.965, '(a)', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
        ax1.text(350, 50.0, r"$68\%$", horizontalalignment='center', verticalalignment = 'center', fontsize = PlotScripts.global_labelsize) 
        ax1.text(350, 50.8, r"$95\%$", horizontalalignment='center', verticalalignment = 'center', fontsize = PlotScripts.global_labelsize)

        leg = ax1.legend(loc='lower right', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize(PlotScripts.global_legendsize)


        plt.tight_layout()
        outputFile = './{0}{1}'.format(output_tag, output_format)
        plt.savefig(outputFile)  # Save the figure
        print('Saved file to {0}'.format(outputFile))
        plt.close()

##

def plot_singleSFR(galaxies_filepath_array, merged_galaxies_filepath_array, number_snapshots, simulation_norm, model_tags, output_tag):
 
    SFR_gal = []
    SFR_ensemble = []

    ejected_gal = []
    ejected_ensemble = []

    infall_gal = []
    infall_ensemble = [] 

    ejectedmass_gal = []
    ejectedmass_ensemble = []

    N_random = 1

    ax1 = plt.subplot(111)
#    ax3 = plt.subplot(122)
    #ax5 = plt.subplot(133)

    look_for_alive = 1
    #idx_array = [20004, 20005, 20016]
    #halonr_array = [7381]
    halonr_array = [389106]
    #halonr_array = [36885]
    for model_number in range(0, len(model_tags)):
        if(simulation_norm[model_number] == 0):
            AllVars.Set_Params_Mysim()
        elif(simulation_norm[model_number] == 3):
            AllVars.Set_Params_Tiamat_extended()
        else:
            print("Simulation norm was set to {0}.".format(simulation_norm[model_number]))
            raise ValueError("This option has been implemented yet.  Get your head in the game Jacob!")

        SFR_gal.append([])
        SFR_ensemble.append([])

        ejected_gal.append([])
        ejected_ensemble.append([])

        infall_gal.append([])
        infall_ensemble.append([])

        ejectedmass_gal.append([])
        ejectedmass_ensemble.append([])

        GG, Gal_Desc = ReadScripts.ReadGals_SAGE_DelayedSN(galaxies_filepath_array[model_number], 0, number_snapshots[model_number], comm) # Read in the correct galaxy file. 
        G_Merged, Merged_Desc = ReadScripts.ReadGals_SAGE_DelayedSN(merged_galaxies_filepath_array[model_number], 0, number_snapshots[model_number], comm) # Also need the merged galaxies.
        G = ReadScripts.Join_Arrays(GG, G_Merged, Gal_Desc) # Then join them together for all galaxies that existed at this Redshift. 


        if look_for_alive == 1:
            G.GridHistory[G.GridHistory >= 0] = 1
            G.GridHistory[G.GridHistory < 0] = 0
            alive = np.sum(G.GridHistory, axis = 1)
#            print "The galaxy that was present in the most snapshots is %d which was in %d snaps" %(np.argmax(alive), np.amax(alive)) 
            most_alive = alive.argsort()[-10:][::-1] # Finds the 3 galaxies alive for the most snapshots.  Taken from https://stackoverflow.com/questions/6910641/how-to-get-indices-of-n-maximum-values-in-a-numpy-array
#            print G.HaloNr[most_alive]
         

        t = np.empty((number_snapshots[model_number])) 
 
       
        for snapshot_idx in range(0, number_snapshots[model_number]): 
            w = np.where((G.GridHistory[:, snapshot_idx] != -1) & (G.GridStellarMass[:, snapshot_idx] > 0.0) & (G.GridStellarMass[:, snapshot_idx] < 1e5) & (G.GridCentralGalaxyMass[:, snapshot_idx] >= m_low_SAGE) & (G.GridCentralGalaxyMass[:, snapshot_idx] <=  m_high_SAGE))[0] # Only include those galaxies that existed at the current snapshot, had positive (but not infinite) stellar/Halo mass and Star formation rate.
            
            SFR_ensemble[model_number].append(np.mean(G.GridSFR[w,snapshot_idx]))
            ejected_ensemble[model_number].append(np.mean(G.GridOutflowRate[w, snapshot_idx]))
            infall_ensemble[model_number].append(np.mean(G.GridInfallRate[w, snapshot_idx]))

            t[snapshot_idx] = (t_BigBang - cosmo.lookback_time(AllVars.SnapZ[snapshot_idx]).value) * 1.0e3 
            
        for p in range(0, N_random):
            random_idx = (np.where((G.HaloNr == halonr_array[p]))[0])[0] 
            SFR_gal[model_number].append(G.GridSFR[random_idx]) # Remember the star formation rate history of the galaxy.
            ejected_gal[model_number].append(G.GridOutflowRate[random_idx])
            infall_gal[model_number].append(G.GridInfallRate[random_idx])
            ejectedmass_gal[model_number].append(G.GridEjectedMass[random_idx])
        
            #SFR_gal[model_number][p][SFR_gal[model_number][p] < 1.0e-15] = 1 
            for snapshot_idx in range(0, number_snapshots[model_number]):  
                if snapshot_idx == 0:
                    pass 
                elif(G.GridHistory[random_idx, snapshot_idx] == -1):
                    SFR_gal[model_number][p][snapshot_idx] = SFR_gal[model_number][p][snapshot_idx - 1]

#        SFR_ensemble[model_number] = np.nan_to_num(SFR_ensemble[model_number])        
#        SFR_ensemble[model_number][SFR_ensemble[model_number] < 1.0e-15] = 1    

         
#        ejected_ensemble[model_number][ejected_ensemble[model_number] < 1.0e-15] = 1     
       
        
        ax1.plot(t, SFR_ensemble[model_number], color = PlotScripts.colors[0], linestyle = PlotScripts.linestyles[model_number], label = model_tags[model_number], linewidth = PlotScripts.global_linewidth)
        ax1.plot(t, ejected_ensemble[model_number], color = PlotScripts.colors[1], linestyle = PlotScripts.linestyles[model_number], linewidth = PlotScripts.global_linewidth, alpha = 1.0)
        #ax5.plot(t, infall_ensemble[model_number], color = PlotScripts.colors[2], linestyle = PlotScripts.linestyles[model_number], linewidth = PlotScripts.global_linewidth, alpha = 1.0)
        #ax5.plot(t, ejectedmass_ensemble[model_number], color = PlotScripts.colors[2], linestyle = PlotScripts.linestyles[model_number], linewidth = PlotScripts.global_linewidth, alpha = 1.0)
        
        for p in range(0, N_random):
            ax1.plot(t, SFR_gal[model_number][p], color = PlotScripts.colors[0], linestyle = PlotScripts.linestyles[model_number], alpha = 0.5, linewidth = 1)
            ax1.plot(t, ejected_gal[model_number][p], color = PlotScripts.colors[1], linestyle = PlotScripts.linestyles[model_number], alpha = 0.5, linewidth = 1)
            #ax5.plot(t, infall_gal[model_number][p], color = PlotScripts.colors[2], linestyle = PlotScripts.linestyles[model_number], alpha = 0.5, linewidth = 1)
            #ax5.plot(t, ejectedmass_gal[model_number][p], color = PlotScripts.colors[2], linestyle = PlotScripts.linestyles[model_number], alpha = 0.5, linewidth = 1)

            #ax1.plot(t, SFR_gal[model_number][p], color = PlotScripts.colors[0], linestyle = PlotScripts.linestyles[model_number], alpha = 1.0, linewidth = 1, label = model_tags[model_number])
            #ax1.plot(t, ejected_gal[model_number][p], color = PlotScripts.colors[1], linestyle = PlotScripts.linestyles[model_number], alpha = 1.0, linewidth = 1, label = model_tags[model_number])

    ax1.plot(np.nan, np.nan, color = 'r', linestyle = '-', label = "SFR")
    ax1.plot(np.nan, np.nan, color = 'b', linestyle = '-', label = "Outflow")

#    exit() 
    #ax1.plot(np.nan, np.nan, color = PlotScripts.colors[0], label = 'SFR')
    #ax1.plot(np.nan, np.nan, color = PlotScripts.colors[1], label = 'Outflow')

    ax1.set_yscale('log', nonposy='clip')
    ax1.set_ylabel(r"$\mathrm{Mass \: Flow} \: [\mathrm{M}_\odot \mathrm{yr}^{-1}]$")
    ax1.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = PlotScripts.global_fontsize)
    ax1.set_xlim(PlotScripts.time_xlim)
    ax1.set_ylim([1e-6, 1e3])

    '''
    ax3.set_yscale('log', nonposy='clip')
    ax3.set_ylabel(r"$\mathrm{Outflow \: Rate} \: [\mathrm{M}_\odot \mathrm{yr}^{-1}]$")
    ax3.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = PlotScripts.global_fontsize)
    ax3.set_xlim(PlotScripts.time_xlim)
    ax3.set_ylim([1e-8, 1e3])

    ax5.set_yscale('log', nonposy='clip')
    #ax5.set_ylabel(r"$\mathrm{Infall \: Rate} \: [\mathrm{M}_\odot \mathrm{yr}^{-1}]$")
    ax5.set_ylabel(r"$\mathrm{Ejected Mass} [\mathrm{M}_\odot]$")
    ax5.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = PlotScripts.global_fontsize)
    ax5.set_xlim(PlotScripts.time_xlim)
    #ax5.set_ylim([1e-8, 1e3])
    ax5.set_ylim([1e6, 1e10])
    '''
    ax2 = ax1.twiny()
    #ax4 = ax3.twiny()
    #ax6 = ax5.twiny()

    t_plot = (t_BigBang - cosmo.lookback_time(PlotScripts.z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
    z_labels = ["$%d$" % x for x in PlotScripts.z_plot] # Properly Latex-ize the labels.

    ax2.set_xlabel(r"$z$", size = PlotScripts.global_labelsize)
    ax2.set_xlim(PlotScripts.time_xlim)
    ax2.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
    ax2.set_xticklabels(z_labels) # But label them as redshifts.

    '''
    ax4.set_xlabel(r"$z$", size = PlotScripts.global_labelsize)
    ax4.set_xlim(PlotScripts.time_xlim)
    ax4.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
    ax4.set_xticklabels(z_labels) # But label them as redshifts.

    ax6.set_xlabel(r"$z$", size = PlotScripts.global_labelsize)
    ax6.set_xlim(PlotScripts.time_xlim)
    ax6.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
    ax6.set_xticklabels(z_labels) # But label them as redshifts.
    '''

    plt.tight_layout()
    leg = ax1.legend(loc='lower right', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize)


    outputFile = './Halo%d_mlow%.2f_%s%s' %(halonr_array[0], m_low_SAGE, output_tag, output_format) 
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print('Saved file to {0}'.format(outputFile))
    plt.close()

##

def plot_quasars_count(SnapList, N_quasars_z, N_quasars_boost_z, N_gal, fesc_prescription, simulation_norm, FirstFile, LastFile, NumFile, model_tags, output_tag):
    '''
    Parameters
    ---------
    SnapList : Nested 'array-like` of ints, SnapList[model_number0] = [snapshot0_model0, ..., snapshotN_model0], with length equal to the number of models.
        Snapshots that we plot the quasar density at for each model.
    N_quasars_z : Nested array of floats, N_quasars_z[model_number0] = [N_quasars_z0, N_quasars_z1, ..., N_quasars_zN]. Outer array has length equal to the number of models, inner array has length equal to length of the model's SnapList. 
        Number of quasars, THAT WENT OFF, during the given redshift. 
    N_quasars_boost_z : Nested array of floats, N_quasars_boost_z[model_number0] = [N_quasars_boost_z0, N_quasars_boost_z1, ..., N_quasars_boost_zN]. Outer array has length equal to the number of models, inner array has length equal to length of the model's SnapList. 
        Number of galaxies that had their escape fraction boosted by quasar activity. 
    N_gal : Nested array of floats, N_gal_z[model_number0] = [N_gal_z0, N_gal_z1, ..., N_gal_zN]. Outer array has length equal to the number of models, inner array has length equal to length of the model's SnapList. 
        Number of galaxies at each redshift.
    fesc_prescription : Array with length equal to the number of models.
        Denotes what escape fraction prescription each model used.  Quasars are only tracked when fesc_prescription == 3.
    simulation_norm : array with length equal to the number of models.
        Denotes which simulation each model uses.  
        0 : MySim
        1 : Mini-Millennium
        2 : Tiamat (down to z = 5)
        3 : Extended Tiamat (down to z = 1.6ish).
        4 : Britton's Simulation
        5 : Kali
    FirstFile, LastFile, NumFile : array of integers with length equal to the number of models.
        The file numbers for each model that were read in (defined by the range between [FirstFile, LastFile] inclusive) and the TOTAL number of files for this model (we may only be plotting a subset of the volume). 
    model_tags : array of strings with length equal to the number of models.
        Strings that contain the tag for each model.  Will be placed on the plot.
    output_tag : string
        Name of the file that will be generated. File will be saved in the current directory with the output format defined by the 'output_format' variable at the beggining of the file.

    Returns
    -------
    No returns.
    Generates and saves the plot (named via output_tag).

    Units
    -----
    No relevant units. 
    '''

    if rank == 0:       
        fig = plt.figure() 
        ax1 = fig.add_subplot(111)  
        ax6 = ax1.twinx()

        fig2 = plt.figure() 
        ax3 = fig2.add_subplot(111)  
        ax5 = ax3.twinx()       
 
    for model_number in range(0, len(SnapList)): # Does this for each of the models. 
        if (fesc_prescription[model_number] != 3): # Want to skip the models that didn't count quasars.
            continue
        ## Normalization for each model. ##
        if (simulation_norm[model_number] == 0):
            AllVars.Set_Params_Mysim()
        elif (simulation_norm[model_number] == 1):
            AllVars.Set_Params_MiniMill()
        elif (simulation_norm[model_number] == 2):
            AllVars.Set_Params_Tiamat()
        elif (simulation_norm[model_number] == 3):
            AllVars.Set_Params_Tiamat_extended()
        elif (simulation_norm[model_number] == 4):
            AllVars.Set_Params_Britton()       
        elif (simulation_norm[model_number] == 5):
            AllVars.Set_Params_Kali()
 
        box_factor = (LastFile[model_number] - FirstFile[model_number] + 1.0)/(NumFile[model_number]) # This factor allows us to take a sub-volume of the box and scale the results to represent the entire box.
        print("We are plotting the quasar density using {0:.4f} of the box's volume.".format(box_factor))
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * box_factor 

        ####    

        ## We perform the plotting on Rank 0 so only this rank requires the final counts array. ##
        if rank == 0:
            quasars_total = np.zeros_like((N_quasars_z[model_number]))
            boost_total = np.zeros_like(N_quasars_boost_z[model_number])
            gal_count_total = np.zeros_like(N_gal[model_number])
        else:
            quasars_total = None
            boost_total = None
            gal_count_total = None
       
        N_quasars_tmp = np.array((N_quasars_z[model_number])) # So we can use MPI.Reduce()
        comm.Reduce([N_quasars_tmp, MPI.DOUBLE], [quasars_total, MPI.DOUBLE], op = MPI.SUM, root = 0) # Sum the number of quasars and passes back to rank 0. 

        N_quasars_boost_tmp = np.array(N_quasars_boost_z[model_number]) # So we can use MPI.Reduce()
        comm.Reduce([N_quasars_boost_tmp, MPI.DOUBLE], [boost_total, MPI.DOUBLE], op = MPI.SUM, root = 0) # Sum the number of galaxies that had their fesc boosted. 

        N_gal_tmp = np.array(N_gal[model_number]) # So we can use MPI.Reduce()
        comm.Reduce([N_gal_tmp, MPI.DOUBLE], [gal_count_total, MPI.DOUBLE], op = MPI.SUM, root = 0) # Sum the number of total galaxies.

        if rank == 0:
            print("Quasars_total {0}".format(quasars_total))
      
        if rank == 0: 
                
            title = model_tags[model_number] 
            t = np.empty(len(SnapList[model_number]))
            for snapshot_idx in range(0, len(SnapList[model_number])):  
                t[snapshot_idx] = (t_BigBang - cosmo.lookback_time(AllVars.SnapZ[SnapList[model_number][snapshot_idx]]).value) * 1.0e3   
                    
            ax1.plot(t, quasars_total / norm, color = PlotScripts.colors[model_number], linestyle = PlotScripts.linestyles[0], rasterized = True, label = title, linewidth = PlotScripts.global_linewidth)

            ax3.plot(t, boost_total, color = PlotScripts.colors[model_number], linestyle = PlotScripts.linestyles[0], rasterized = True, label = title, linewidth = PlotScripts.global_linewidth) 
            w = np.where((gal_count_total > 0.0))[0] # Since we're doing a division, need to only plot those redshifts that actually have galaxies.

            ax5.plot(t[w], np.divide(boost_total[w], gal_count_total[w]), color = PlotScripts.colors[model_number], linestyle = PlotScripts.linestyles[1], rasterized = True, linewidth = PlotScripts.global_linewidth) 

            ax6.plot(t[w], gal_count_total[w] / norm, color = PlotScripts.colors[model_number], linestyle = PlotScripts.linestyles[1], rasterized = True, linewidth = PlotScripts.global_linewidth)
            
            ax1.plot(np.nan, np.nan, color = 'k', linestyle = PlotScripts.linestyles[0], label = "Quasar Density")
            ax1.plot(np.nan, np.nan, color = 'k', linestyle = PlotScripts.linestyles[1], label = "Galaxy Density")

            ax3.plot(np.nan, np.nan, color = 'k', linestyle = PlotScripts.linestyles[0], label = "Count")
            ax3.plot(np.nan, np.nan, color = 'k', linestyle = PlotScripts.linestyles[1], label = "Fraction of Galaxies")

            ax1.xaxis.set_minor_locator(mtick.MultipleLocator(PlotScripts.time_tickinterval))
            ax1.set_xlim(PlotScripts.time_xlim)
            ax1.set_yscale('log', nonposy='clip')

            ax3.xaxis.set_minor_locator(mtick.MultipleLocator(PlotScripts.time_tickinterval))
            ax3.set_xlim(PlotScripts.time_xlim)
            ax3.set_yscale('log', nonposy='clip')

            ## Create a second axis at the top that contains the corresponding redshifts. ##
            ## The redshift defined in the variable 'z_plot' will be displayed. ##
            ax2 = ax1.twiny()
            ax4 = ax3.twiny()

            t_plot = (t_BigBang - cosmo.lookback_time(PlotScripts.z_plot).value) * 1.0e3 # Corresponding time values on the bottom.
            z_labels = ["$%d$" % x for x in PlotScripts.z_plot] # Properly Latex-ize the labels.

            ax2.set_xlabel(r"$z$", size = PlotScripts.global_labelsize) 
            ax2.set_xlim(PlotScripts.time_xlim)
            ax2.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
            ax2.set_xticklabels(z_labels) # But label them as redshifts.

            ax4.set_xlabel(r"$z$", size = PlotScripts.global_labelsize) 
            ax4.set_xlim(PlotScripts.time_xlim)
            ax4.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
            ax4.set_xticklabels(z_labels) # But label them as redshifts.

            ax1.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = PlotScripts.global_labelsize) 
            ax1.set_ylabel(r'$N_\mathrm{Quasars} \: [\mathrm{Mpc}^{-3}]$', fontsize = PlotScripts.global_fontsize) 
            ax6.set_ylabel(r'$N_\mathrm{Gal} \: [\mathrm{Mpc}^{-3}]$', fontsize = PlotScripts.global_fontsize) 

            ax3.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = PlotScripts.global_labelsize) 
            ax3.set_ylabel(r'$N_\mathrm{Boosted}$', fontsize = PlotScripts.global_fontsize) 
            ax5.set_ylabel(r'$\mathrm{Fraction \: Boosted}$', fontsize = PlotScripts.global_fontsize) 

            leg = ax1.legend(loc='lower left', numpoints=1, labelspacing=0.1)
            leg.draw_frame(False)  # Don't want a box frame
            for t in leg.get_texts():  # Reduce the size of the text
                t.set_fontsize(PlotScripts.global_legendsize)

            leg = ax3.legend(loc='lower left', numpoints=1, labelspacing=0.1)
            leg.draw_frame(False)  # Don't want a box frame
            for t in leg.get_texts():  # Reduce the size of the text
                t.set_fontsize(PlotScripts.global_legendsize)

            fig.tight_layout()            
            fig2.tight_layout()           
 
            outputFile1 = './{0}_quasardensity_{1}'.format(output_tag, output_format)
            outputFile2 = './{0}_boostedcount_{1}'.format(output_tag, output_format)

            fig.savefig(outputFile1)  # Save the figure
            fig2.savefig(outputFile2)  # Save the figure

            print("Saved to {0}".format(outputFile1))
            print("Saved to {0}".format(outputFile2))

            plt.close(fig)
            plt.close(fig2)


##

def plot_photon_quasar_fraction(snapshot, filenr, output_tag, QuasarFractionalPhoton, QuasarActivityToggle, NumSubsteps):

    ax1 = plt.subplot(111)

    counts, bin_edges, bin_middle = AllVars.Calculate_Histogram(QuasarFractionalPhoton, 0.05, 0, 0, 1) 
            
    ax1.plot(bin_middle, counts, lw = PlotScripts.global_linewidth, color = 'r')
   
    ax1.axvline(np.mean(QuasarFractionalPhoton), lw = 0.5, ls = '-')
 
    ax1.set_yscale('log', nonposy='clip')
    ax1.set_xlabel(r"$\mathrm{Fractional \: Photon \: Boost}$")
    ax1.set_ylabel(r"$\mathrm{Count}$")

    ax1.set_ylim([1e1, 1e5])

    outputFile1 = './photonfraction/file{0}_snap{1}_{2}{3}'.format(filenr, snapshot, output_tag, output_format)

    plt.tight_layout()
    plt.savefig(outputFile1)

    print("Saved to {0}".format(outputFile1))
    
    plt.close()

### Here ends the plotting functions. ###
### Here begins the functions that calculate various properties for the galaxies (fesc, Magnitude etc). ###

def Calculate_HaloPartStellarMass(halo_part, stellar_mass, bound_low, bound_high):
    '''
    Calculates the stellar mass for galaxies whose host halos contain a specified number of particles.

    Parameters
    ----------
    halo_part : array
        Array containing the number of particles inside each halo.
    stellar_mass : array
        Array containing the Stellar Mass for each galaxy (entries align with HaloPart). Units of log10(Msun).
    bound_low, bound_high : int
        We calculate the Stellar Mass of galaxies whose host halo has, bound_low <= halo_part <= bound_high.    

    Return
    -----
    mass, mass_std : float
        Mean and standard deviation stellar mass of galaxies whose host halo has number of particles between the specified bounds.  Units of log10(Msun)

    Units
    -----
    Input Stellar Mass is in units of log10(Msun).
    Output mean/std Stellar Mass is in units of log10(Msun).
    '''

    w = np.where((halo_part >= bound_low) & (halo_part <= bound_high))[0] # Find the halos with particle number between the bounds.

    mass = np.mean(10**(stellar_mass[w]))
    mass_std = np.std(10**(stellar_mass[w]))

    return np.log10(mass), np.log10(mass_std)

##


def calculate_fesc(fesc_prescription, fesc_normalization, halo_mass, ejected_fraction, quasar_activity_toggle, quasar_fractional_boost, old):
    '''
    Calculate the escape fraction for a given prescription.

    Parameters
    ----------
    fesc_prescription : int
        Number that controls what escape fraction prescription we are using.
        0 : Constant, fesc = Constant.
        1 : Scaling with Halo Mass, fesc = A*Mh^B.
        2 : Scaling with ejected fraction, fesc = fej*A + B.
        3 : Depending on quasar activity, fesc = A (if no quasar within 1 dynamical time) or fesc = B (if quasar event within 1 dynamical time).
    fesc_normalization : float (if fesc_prescription == 0) or `numpy.darray' with length 2 (if fesc_prescription == 1 or == 2).
        If fesc_prescription == 0, gives the constant value for the escape fraction.
        If fesc_prescription == 1 or == 2 or ==3, gives A and B with the form [A, B].
    halo_mass : `numpy.darray'
        Array that contains the halo masses for this snapshot. Units are log10(Msun).
    ejected_fraction : `numpy.darray'
        Array that contains the ejected fraction of the galaxies for this snapshot.
    quasar_activity_toggle : `numpy.darray'
        Array that contains whether there has been a quasar event within 1 dynamical time for each galaxy.
    quasar_fractional_boost : Array of floats with length equal to the number of galaxies.
        Array that contains the fraction of photons that should be boosted with the quasar activity.
        If the escape fraction should not be boosted by quasar activity, the index value with be 0.0. 
    
    Returns
    -------
    fesc : array, length is same as input arrays.
    Array containing the escape fraction for each galaxy for this snapshot.

    Units
    -----
    Input Halo Mass is in units of log10(Msun).
    '''

    if(len(halo_mass) != len(ejected_fraction)):
        print("The length of the halo_mass array is {0} and the length of the ejected_fraction array is {1}".format(len(halo_mass), len(ejected_fraction)))
        raise ValueError("These two should be equal.")

    if(fesc_prescription > 3):
        print("The prescription value you've chosen is {0}".format(fesc_prescription))
        raise ValueError("Currently we only have prescriptions 0 (constant), 1 (scaling with halo mass), 2 (scaling with ejected fraction) and 3 (boosted by  quasar activity).")

    if((fesc_prescription == 0 and isinstance(fesc_normalization, list) == True)):
        print("The escape fraction prescription is 0 but fesc_normalization was a list with value of {0}".format(fesc_normalization))
        raise ValueError("For a constant escape fraction, fesc_noramlization should be a single float.")

    if((fesc_prescription == 1 or fesc_prescription == 2) and (isinstance(fesc_normalization, list) == False)):
        print("The escape fraction prescription is {0} but fesc_normalization was not a list; it instead had a value of {1}".format(fesc_normalization))
        raise ValueError("For a scaling escape fraction fesc_normalization should be a list of the form [A, B]") 
 
    if (fesc_prescription == 0):
        fesc = np.full((len(halo_mass)), fesc_normalization)
    elif (fesc_prescription == 1):
        fesc = fesc_normalization[0]*pow(10,halo_mass*fesc_normalization[1])
    elif (fesc_prescription == 2):
        fesc = fesc_normalization[0]*ejected_fraction + fesc_normalization[1]
    elif (fesc_prescription == 3):
        if (old == 1):
            fesc = np.full((len(halo_mass)), fesc_normalization[0]) # Initialize to the base value.
            fesc[quasar_activity_toggle == 1] = fesc_normalization[1] # Set the recent quasar galaxies to the boosted value.
        else:
            fesc = fesc_normalization[0] * (1 - quasar_fractional_boost)  + fesc_normalization[1] * quasar_fractional_boost # This is a tad subtle so I'll explain. 
        # quasar_fractional_boost is 0.0 if the escape fraction should not be boosted.  In this case, the escape fraction will be the base fraction.
        # If the quasar should boost the escape fraction for the entire snapshot, then quasar_fractional_boost = 1 and the escape fraction becomes the boosted fraction.
        # However if the quasar should boost the escape fraction for some fraction of the snapshot, then we perform the weighting and assign the escape fraction appropriately.
        # For example, consider a galaxy that emits 10^50 photons/s with a base escape fraction of 0.2 and boosted value of 1.0. 
        ## If we want three quarters of these photons to have the boosted escape fraction, then what should happen is (10^50 * (0.2 * 1/4) + 10^50 * (1.0 * 3/4)).
        ## In essence this is the exact same as assigning a single escape fraction of 0.2 * 1/4 + 1.0 * 3/4 = 0.8.
 
    if len(fesc) != 0:
    ## Adjust bad values, provided there isn't a riduculous amount of them. ##  
        w = np.where((halo_mass >= m_low) & (halo_mass <= m_high))[0] # Only care about the escape fraction values between these bounds. 
        nan_values = len(fesc[w][np.isnan(fesc[w])]) 
        aboveone_values = len(fesc[w][fesc[w] > 1.0])
        belowzero_values = len(fesc[w][fesc[w] < 0.0])
        bad_ratio =  (float(nan_values) + float(aboveone_values) + float(belowzero_values)) / float(len(fesc))
        
        if (bad_ratio > 0.10):
            print("The ratio of bad escape fractions to good is {0:.4f}".format(bad_ratio))
            print("Values above 1: {0}".format(fesc[fesc > 1.0]))
            w = np.where(fesc > 1.0)[0]
            print("Halo mass is: {0}".format(halo_mass[w]))
            print("Ejected fraction is: {0}".format(ejected_fraction[w]))
            
            print("Values below 0: {0}".format(fesc[fesc < 0.0]))
            w = np.where(fesc < 0.0)[0]
            print("Halo mass is: {0}".format(halo_mass[w]))
            print("Ejected fraction is: {0}".format(ejected_fraction[w]))
            if (len(fesc) > 100): # Only throw an error if we had an appreciable number of bad fesc.
                raise ValueError("This was above the tolerance level of 10%.")

        fesc[np.isnan(fesc)] = 0 # Get rid of any lingering Nans.
        fesc[fesc > 1] = 1.0  # Get rid of any values above 1. 
        fesc[fesc < 0] = 0.0  # Get rid of any values below 0.
    else:
        fesc = np.nan   
     
    return fesc

##
def calculate_photons(SFR, Z):
    '''
    This calculates the number of HI ionizing photons emitted by a galaxy of a given star formation rate and metallicity.
    Fit is based on the Stellar Synthesis of STARBURST99 (Leither et al., 1999).

    Parameters
    ----------
    SFR : `numpy.darray' with length equal to the number of galaxies.
    Star formation rate for each galaxy. Units of log10(Msun yr^-1).
    Z : `numpy.darray' with length equal to the number of galaxies.
    Metallicity (pure, not solar) for each galaxy.

    Returns
    -------
    ngamma_HI : `numpy.darray' with length equal to the number of galaxies.
    Number of HI ionizing photons emitted for each galaxy. Units of log10(s^-1).

    Units
    -----
    Star Formation Rate is in units of log10(Msun yr^-1).
    Metallicity is in proper metallicity (not solar).
    ngamma_HI is in units of log10(s^-1).
    '''

    ngamma_HI = []

    ## Fits are based on the blah tracks of STARBURST99. ##

    for i in range(0, len(SFR)):
        ngamma_HI_tmp = 0.0
        if (SFR[i] == 0.0):
            n_gamma_HI_tmp = 0
        elif (Z[i] < 0.0025):
            ngamma_HI_tmp = SFR[i] + 53.354
        elif (Z[i] >= 0.0025 and Z[i] < 0.006):
            ngamma_HI_tmp = SFR[i] + 53.290
        elif (Z[i] >= 0.006 and Z[i] < 0.014):
            ngamma_HI_tmp = SFR[i] + 53.240
        elif (Z[i] >= 0.014 and Z[i] < 0.30):
            ngamma_HI_tmp = SFR[i] + 53.166
        else:
            ngamma_HI_tmp = SFR[i] + 53.041  
        
        ngamma_HI.append(ngamma_HI_tmp)
    return ngamma_HI

##
def calculate_UV_extinction(z, L, M):
    '''
    Calculates the observed UV magnitude after dust extinction is accounted for.

    Parameters
    ----------
    z : float
    Redshift we are calculating the extinction at.
    L, M : array, length equal to the number of galaxies at this snapshot.
    Array containing the UV luminosities and magnitudes.

    Returns
    -------
    M_UV_obs : array, length equal to the number of galaxies at this snapshot.
    Array containing the observed UV magnitudes.

    Units
    -----
    Luminosities are in units of log10(erg s^-1 A^-1).
    Magnitudes are in the AB system.
    '''

    M_UV_bins = np.arange(-24, -16, 0.1)
    A_mean = np.zeros((len(MUV_bins))) # A_mean is the average UV extinction for a given UV bin.    

    for j in range(0, len(M_UV_bins)):
        beta = calculate_beta(M_UV_bins[j], AllVars.SnapZ[current_snap]) # Fits the beta parameter for the current redshift/UV bin. 
        dist = np.random.normal(beta, 0.34, 10000) # Generates a normal distribution with mean beta and standard deviation of 0.34.
        A = 4.43 + 1.99*dist 
        A[A < 0] = 0 # Negative extinctions don't make sense.
            
        A_Mean[j] = np.mean(A)

    indices = np.digitize(M, M_UV_bins) # Bins the simulation magnitude into the MUV bins. Note that digitize defines an index i if bin[i-1] <= x < bin[i] whereas I prefer bin[i] <= x < bin[i+1]
    dust = A_Mean[indices]
    flux = AllVars.Luminosity_to_Flux(L, 10.0) # Calculate the flux from a distance of 10 parsec, units of log10(erg s^-1 A^-1 cm^-2). 
    flux_observed = flux - 0.4*dust
    
    f_nu = ALlVars.spectralflux_wavelength_to_frequency(10**flux_observed, 1600) # Spectral flux desnity in Janksy.
    M_UV_obs(-2.5 * np.log10(f_nu) + 8.90) # AB Magnitude from http://www.astro.ljmu.ac.uk/~ikb/convert-units/node2.html

    return M_UV_obs

##

def update_cumulative_stats(mean_pool, std_pool, N_pool, mean_local, std_local, N_local):
    '''
    Update the cumulative statistics (such as Stellar Mass Function, Mvir-Ngamma, fesc-z) that are saved across files.
    Pooled mean formulae taken : from https://www.ncbi.nlm.nih.gov/books/NBK56512/
    Pooled variance formulae taken from : https://en.wikipedia.org/wiki/Pooled_variance

    Parameters
    ----------
    mean_pool, std_pool, N_pool : array of floats with length equal to the number of bins (e.g. the mass bins for the Stellar Mass Function).
        The current mean, standard deviation and number of data points within in each bin.  This is the array that will be updated in this function.
    mean_local, std_local, N_local : array of floats with length equal to the number of bins.
        The mean, standard deviation and number of data points within in each bin that will be added to the pool.

    Returns
    -------
    mean_pool, std_pool, N_pool : (See above)
    The updated arrays with the local values added and accounted for within the pools.

    Units
    -----
    All units are kept the same as the input units.
    Values are in real-space (not log-space).
    '''
   
    N_times_mean_local = np.multiply(N_local, mean_local)
    N_times_var_local = np.multiply(N_local - 1, np.multiply(std_local, std_local)) # Actually N - 1 because of Bessel's Correction 
                                        # https://en.wikipedia.org/wiki/Bessel%27s_correction).  #
    N_times_mean_pool = np.add(N_times_mean_local, np.multiply(N_pool, mean_pool))
    N_times_var_pool = np.add(N_times_var_local, np.multiply(N_pool - 1, np.multiply(std_pool, std_pool)))
    N_pool = np.add(N_local, N_pool)

    '''
    print(mean_local)
    print(type(mean_local))
    print((type(mean_local).__module__ == np.__name__))
    print(isinstance(mean_local, list))
    print(isinstance(mean_local,float64))
    print(isinstance(mean_local,float32))
    '''
    if (((type(mean_local).__module__ == np.__name__) == True or (isinstance(mean_local, list) == True)) and isinstance(mean_local, float) == False and isinstance(mean_local, int) == False and isinstance(mean_local,float32) == False and isinstance(mean_local, float64) == False): # Checks to see if we are dealing with arrays. 
        for i in range(0, len(N_pool)):
            if(N_pool[i] == 0): # This case is when we have no data points in the bin. 
                mean_pool[i] = 0.0
            else:
                mean_pool[i] = N_times_mean_pool[i]/N_pool[i]
            if(N_pool[i] < 3): # In this instance we don't have enough data points to properly calculate the standard deviation.
                std_pool[i] = 0.0
            else:
                std_pool[i] = np.sqrt(N_times_var_pool[i]/ (N_pool[i] - 2)) # We have -2 because there is two instances of N_pool contains two 'N - 1' terms. 
        
    else:
        mean_pool = N_times_mean_pool / N_pool

        if(N_pool < 3):
            std_pool = 0.0
        else:
            std_pool = np.sqrt(N_times_var_pool / (N_pool - 2))
 
    return mean_pool, std_pool

    ### Here ends the functions that deal with galaxy data manipulation. ###

#################################
### Here starts the main body of the code. ###

np.seterr(divide='ignore')
number_models = 1

galaxies_model1 = '/lustre/projects/p004_swin/jseiler/january/galaxies/tiamat_IRA_reion_quasarsubstep_10step_z1.827'
galaxies_model6 = '/lustre/projects/p004_swin/jseiler/kali/galaxies/kali_fiducial_quasarsubstep_z5.782'
galaxies_model7 = '/lustre/projects/p004_swin/jseiler/kali/galaxies/kali_fiducial_quasarsubstep_15step_z5.782'
galaxies_model8 = '/lustre/projects/p004_swin/jseiler/kali/galaxies/kali_fiducial_quasarsubstep_1step_z5.782'


merged_galaxies_model1 = '/lustre/projects/p004_swin/jseiler/january/galaxies/tiamat_IRA_reion_quasarsubstep_10step_MergedGalaxies'
merged_galaxies_model6 = '/lustre/projects/p004_swin/jseiler/kali/galaxies/kali_fiducial_quasarsubstep_MergedGalaxies'
merged_galaxies_model7 = '/lustre/projects/p004_swin/jseiler/kali/galaxies/kali_fiducial_quasarsubstep_15step_MergedGalaxies'
merged_galaxies_model8 = '/lustre/projects/p004_swin/jseiler/kali/galaxies/kali_fiducial_quasarsubstep_1step_MergedGalaxies'

galaxies_filepath_array = [galaxies_model6]
#galaxies_filepath_array = [galaxies_model1, galaxies_model2, galaxies_model3, galaxies_model4]
#galaxies_filepath_array = [galaxies_model3, galaxies_model4]
#galaxies_filepath_array = [galaxies_model8, galaxies_model6, galaxies_model7]

merged_galaxies_filepath_array = [merged_galaxies_model6] 
#merged_galaxies_filepath_array = [merged_galaxies_model1, merged_galaxies_model2, merged_galaxies_model3, merged_galaxies_model4]
#merged_galaxies_filepath_array = [merged_galaxies_model3, merged_galaxies_model4]
#merged_galaxies_filepath_array = [merged_galaxies_model8, merged_galaxies_model6, merged_galaxies_model7]

number_substeps = [10]
#number_snapshots = [92, 164] # Number of snapshots in the simulation (we don't have to do calculations for ALL snapshots).
#number_snapshots = [164, 164, 164, 164] # Number of snapshots in the simulation (we don't have to do calculations for ALL snapshots).
number_snapshots = [99] # Number of snapshots in the simulation (we don't have to do calculations for ALL snapshots).
# Tiamat extended has 164 snapshots.
#FirstFile = [0, 0, 0, 0] # The first file number THAT WE ARE PLOTTING.
#FirstFile = [0, 0] # The first file number THAT WE ARE PLOTTING.
FirstFile = [0] # The first file number THAT WE ARE PLOTTING.
#FirstFile = [0, 0] # The first file number THAT WE ARE PLOTTING.
#LastFile = [26, 26, 26, 26] # The last file number THAT WE ARE PLOTTING.
LastFile = [0] # The last file number THAT WE ARE PLOTTING.
#LastFile = [5, 1] # The last file number THAT WE ARE PLOTTING.
#NumFile = [27, 27, 27, 27] # The number of files for this simulation (plotting a subset of these files is allowed). 
NumFile = [64] # The number of files for this simulation (plotting a subset of these files is allowed). 
#NumFile = [64, 27] # The number of files for this simulation (plotting a subset of these files is allowed). 

same_files = [0] # In the case that model 1 and model 2 (index 0 and 1) have the same files, we don't want to read them in a second time.
# This array will tell us if we should keep the files for the next model or otherwise throw them away and keep them.
# The files will be kept until same_files[current_model_number] = 0.
# For example if we had 5 models we were plotting and model 1, 2, 3 shared the same files and models 4, 5 shared different files,
# Then same_files = [1, 1, 0, 1, 0] would be the correct values.
done_model = np.zeros((number_models)) # We use this to keep track of if we have done a model already.

#model_tags = [r"$f_\mathrm{esc} = 0.2$", r"$f_\mathrm{esc} \: \propto \: \mathrm{Quasar} \: (\mathrm{Boost} = 1.0, N_\mathrm{Dynamical} = 0.1)$", r"$f_\mathrm{esc} \: \propto \: \mathrm{Quasar} \: (\mathrm{Boost} = 1.0, N_\mathrm{Dynamical} = 1)$", r"$f_\mathrm{esc} \: \propto \: \mathrm{Quasar} \: (\mathrm{Boost} = 1.0, N_\mathrm{Dynamical} = 3)$"]
#model_tags = [r"$f_\mathrm{esc} \: \propto \: \mathrm{quasar} (\mathrm{boost} = 1.0, N_\mathrm{dynamical} = 1)$", r"$f_\mathrm{esc} \: \propto \: f_\mathrm{ej}$"]
model_tags = [r"$\mathrm{Kali}$"]
#model_tags = [r"$\mathrm{Tiamat}$"]
#model_tags = [r"$\mathrm{Kali \: 1 \: Step}$", r"$\mathrm{Kali \: 10 \: Step}$"]
#model_tags = [r"$\mathrm{ReionOn}$", r"$\mathrm{NoReion}$"]
#model_tags = [r"Pip Delayed", r"Tiamat Delayed"]

save_tags = ["Kali_Quasar_Newest"]
for model_number in range(0,number_models):
    assert(LastFile[model_number] - FirstFile[model_number] + 1 >= size)

## Constants used for each model. ##
# Need to add an entry for EACH model. #

sSFR_min = [1.0e100, 1.0e100, 1.0e100, 1.0e100]
sSFR_max = [-1.0e100, -1.0e100, 1.0e100, 1.0e100]
halo_cut = [1, 1, 1, 1] # Only calculate galaxy properties whose host halo has particle number greater than this.
source_efficiency = [1, 1, 1, 1] # Used for the halo based prescription for ionizing emissivity.

fesc_lyman_alpha = [0.3, 0.3, 0.3, 0.3] # Escape fraction of Lyman Alpha photons.
#fesc_prescription = [0, 3, 3, 3] # 0 is constant, 1 is scaling with halo mass, 2 is scaling with ejected fraction, 3 is a boosted escape fraction depending upon quasar activity.
fesc_prescription = [3] # 0 is constant, 1 is scaling with halo mass, 2 is scaling with ejected fraction, 3 is a boosted escape fraction depending upon quasar activity.
#fesc_prescription = [0, 0] # 0 is constant, 1 is scaling with halo mass, 2 is scaling with ejected fraction, 3 is a boosted escape fraction depending upon quasar activity.

## Normalizations for the escape fractions. ##
# For prescription 0, requires a number that defines the constant fesc.
# For prescription 1, fesc = A*M^B. Requires an array with 2 numbers the first being A and the second B.
# For prescription 2, fesc = A*fej + B.  Requires an array with 2 numbers the first being A and the second B.
# For prescription 3, fesc = A (if no quasar within 1 dynamical time) or fesc = B (if quasar event within 1 dynamical time). Requires an array with 2 numbers the first being A and the second B.

#fesc_normalization = [[0.4, 0.2], [0.4, 0.2]]
#fesc_normalization = [[0.2, 1.0], [0.4, 0.2]]
fesc_normalization = [[0.2, 1.0]]
#fesc_normalization = [0.2, 0.2]
#fesc_normalization = [0.2, [0.2, 1.0], [0.2, 1.0], [0.2, 1.0]]
N_dynamicaltime = [1] # For escape fraction prescription == 3, this defines how many dynamical times we wait until we turn off the boosted escape fraction.


#N_dynamicaltime = [0, 0.1, 1, 3] # For escape fraction prescription == 3, this defines how many dynamical times we wait until we turn off the boosted escape fraction.
##

#SnapList =  [[53, 66, 80], [51, 64, 78]]
#SnapList = [[64, 76, 93], [51, 64, 79]]
#SnapList = [[51, 64, 78]]
#SnapList = [np.arange(0, 101), np.arange(0, 101), np.arange(0, 101), np.arange(0, 101)]
# z = [6, 7, 8] are snapshots [78, 64, 51]
#PlotSnapList = [[51, 64, 78]]
PlotSnapList = [[64, 76, 93]]
SnapList = [np.arange(0, 99)]
#PlotSnapList = [[45, 48, 52, 55, 57], [45, 48, 52, 55, 57]]
#PlotSnapList = [[42, 43, 44], [42, 43, 44]]
#PlotSnapList = [[90], [90]]
simulation_norm = [5] # 0 for MySim, 1 for Mini-Millennium, 2 for Tiamat (up to z =5), 3 for extended Tiamat (down to z = 1.6ish), 4 for Britton's Sim, 5 for Kali.

stellar_mass_halolen_lower = [32, 95, 95, 95] # These limits are for the number of particles in a halo.  
stellar_mass_halolen_upper = [50, 105, 105, 105] # We calculate the average stellar mass for galaxies whose host halos have particle count between these limits.
calculate_observed_LF = 0
##############################################################################################################

files = [galaxies_filepath_array, merged_galaxies_filepath_array, number_snapshots, model_tags, SnapList, save_tags]
goodness_check = [len(x) == len(y) for i,x in enumerate(files) for j,y in enumerate(files) if i != j] # This goes through the files array and checks to see if all the files have the same lengths.

if False in goodness_check: # If any of the arrays had different lengths, throw an error.
    print(galaxies_filepath_array)
    print(merged_galaxies_filepath_array)
    print(number_snapshots)
    print(model_tags)
    print(SnapList)
    raise ValueError("All the printed arrays (galaxies_filepath_array, merged_galaxies_filepath_array, number_snapshots, SnapList and model_tags respectively) should have the same length.")

if number_models != len(galaxies_filepath_array):
    print("The number of models was given as {0} whereas we had {1} filepaths given.".format(number_models, len(galaxies_filepath_array)))

#bin_Simfast_halos(np.arange(13.500, 5.9, -0.250), AllVars.SnapZ, "/lustre/projects/p004_swin/jseiler/Simfast21/Halos/", "/lustre/projects/p004_swin/jseiler/tiamat/halo_ngamma", fesc_prescription, fesc_normalization, 512, "Simfast_nion")
#exit()

## Arrays for functions of stellar mass. ##
SMF = []
mean_fesc_galaxy_array = []
std_fesc_galaxy_array = []
N_galaxy_array = []

## Arrays for functions of halo mass. ##
mean_ejected_halo_array = []
std_ejected_halo_array = []
mean_fesc_halo_array = []
std_fesc_halo_array = []
mean_Ngamma_halo_array = []
std_Ngamma_halo_array = []
N_halo_array = []

## Arrays for functions of redshift. ##
sum_Ngamma_z_array = []
mean_fesc_z_array = []
std_fesc_z_array = []
N_z = []
galaxy_halo_mass_mean = []
number_galaxies = []
N_quasars_z = [] # This tracks how many quasars went off during a specified snapshot.
N_quasars_boost_z = [] # This tracks how many galaxies are having their escape fraction boosted by quasar activity.
dynamicaltime_quasars_mean_z = []
dynamicaltime_quasars_std_z = []
dynamicaltime_all_mean_z = []
dynamicaltime_all_std_z = []

AllVars.Set_Params_Tiamat_extended()
#for i in range(0, len(AllVars.SnapZ)-1):
#    print "Snapshot ", i, "and Snapshot", i + 1, "have a time difference of", (AllVars.Lookback_Time[i] - AllVars.Lookback_Time[i+1]) * 1.0e3, "Myr"

#plot_singleSFR(galaxies_filepath_array, merged_galaxies_filepath_array, number_snapshots, simulation_norm, model_tags, "singleSFR_Croatia")



######################################################################   
##################### SETTING UP ARRAYS ##############################
######################################################################   
for model_number in range(number_models):

    if(simulation_norm[model_number] == 1):
        AllVars.Set_Params_Mysim()
    elif(simulation_norm[model_number] == 3):
        AllVars.Set_Params_Tiamat_extended()
    elif(simulation_norm[model_number] == 4):
        AllVars.Set_Params_Britton()
    elif(simulation_norm[model_number] == 5):
        AllVars.Set_Params_Kali()
    else: 
        print("Simulation norm was set to {0}.".format(simulation_norm[model_number]))
        raise ValueError("This option has been implemented yet.  Get your head in the game Jacob!")

    if (number_snapshots[model_number] != len(AllVars.SnapZ)): # Here we do a check to ensure that the simulation we've defined correctly matches the number of snapshots we have also defined. 
        print("The number_snapshots array is {0}".format(number_snapshots))
        print("The simulation_norm array is {0}".format(simulation_norm))
        print("The number of snapshots for model_number {0} has {1} but you've said there is only {2}".format(model_number, len(AllVars.SnapZ), number_snapshots[model_number]))
        raise ValueError("Check either that the number of snapshots has been defined properly and that the normalization option is correct.")

    ## The statistics arrays will be 2D arrays.  E.g., Calling SMF[model_number] will yield an array with the count of the number of galaxies within each of the stellar mass bins for each snapshot. ##

    ## Galaxy Arrays ##
    SMF.append([])
    mean_fesc_galaxy_array.append([])
    std_fesc_galaxy_array.append([])
    N_galaxy_array.append([])

    ## Halo arrays. ##
    mean_ejected_halo_array.append([])
    std_ejected_halo_array.append([])
    mean_fesc_halo_array.append([])
    std_fesc_halo_array.append([])
    mean_Ngamma_halo_array.append([])
    std_Ngamma_halo_array.append([])
    N_halo_array.append([])

    ## Redshift arrays. ##
    sum_Ngamma_z_array.append([])
    mean_fesc_z_array.append([])
    std_fesc_z_array.append([])
    N_z.append([])
    galaxy_halo_mass_mean.append([])
    number_galaxies.append([])

    N_quasars_z.append([])
    N_quasars_boost_z.append([])
    dynamicaltime_quasars_mean_z.append([])
    dynamicaltime_quasars_std_z.append([])
    dynamicaltime_all_mean_z.append([])
    dynamicaltime_all_std_z.append([])

    for snapshot_idx in range(len(SnapList[model_number])): # These arrays are used for cumulative values across all files.
        ## Functions of stellar mass arrays. ##
        SMF[model_number].append(np.zeros((NB_gal), dtype = np.float32)) # Stellar Mass Function for each snapshot.
        mean_fesc_galaxy_array[model_number].append(np.zeros((NB_gal), dtype = np.float32)) # Escape fraction as a function of stellar mass. 
        std_fesc_galaxy_array[model_number].append(np.zeros((NB_gal), dtype = np.float32))
        N_galaxy_array[model_number].append(np.zeros((NB_gal), dtype = np.float32)) # How many galaxies have been binned (kept as float to ensure correct division).

        ## Function of halo mass arrays. ##
        mean_ejected_halo_array[model_number].append(np.zeros((NB), dtype = np.float32)) # Ejected fraction as a function of halo mass.
        std_ejected_halo_array[model_number].append(np.zeros((NB), dtype = np.float32)) 
        mean_fesc_halo_array[model_number].append(np.zeros((NB), dtype = np.float32)) # Escape fraction as a function of halo mass.
        std_fesc_halo_array[model_number].append(np.zeros((NB), dtype = np.float32))
        mean_Ngamma_halo_array[model_number].append(np.zeros((NB), dtype = np.float32)) # Number of escaping ionizing photons as a function of halo mass.
        std_Ngamma_halo_array[model_number].append(np.zeros((NB), dtype = np.float32))
        N_halo_array[model_number].append(np.zeros((NB), dtype = np.float32)) # How many halos have been binned to calculate the statistics.
 
        ## Function of Redshift arrays. ##
        sum_Ngamma_z_array[model_number].append(0.0) # This is the sum of ionizing photons emitted at each redshift.
        mean_fesc_z_array[model_number].append(0.0) # Escape fraction (of all galaxies) at each redshift.
        std_fesc_z_array[model_number].append(0.0)
        N_z[model_number].append(0.0) # Total number of galaxies at each redshift.
        galaxy_halo_mass_mean[model_number].append(0.0)
        number_galaxies[model_number].append(0.0)   
        N_quasars_z[model_number].append(0.0) # Total number of quasars that went off at each redshift.
        N_quasars_boost_z[model_number].append(0.0)
        dynamicaltime_quasars_mean_z[model_number].append(0.0)
        dynamicaltime_quasars_std_z[model_number].append(0.0)
        dynamicaltime_all_mean_z[model_number].append(0.0)
        dynamicaltime_all_std_z[model_number].append(0.0)


######################################################################   
#################### ALL ARRAYS SETUP ################################
######################################################################   

for model_number in range(number_models):

    if(simulation_norm[model_number] == 1):
        AllVars.Set_Params_Mysim()
    elif(simulation_norm[model_number] == 3):
        AllVars.Set_Params_Tiamat_extended()
    elif(simulation_norm[model_number] == 4):
        AllVars.Set_Params_Britton()
    elif(simulation_norm[model_number] == 5):
        AllVars.Set_Params_Kali()
    else: 
        print("Simulation norm was set to {0}.".format(simulation_norm[model_number]))
        raise ValueError("This option has been implemented yet.  Get your head in the game Jacob!")

    print("Done_model = {0}".format(done_model))
    if (done_model[model_number] == 1):
        assert(FirstFile[model_number] == FirstFile[model_number - 1]) 
        assert(LastFile[model_number] == LastFile[model_number - 1]) 
        continue
    
    for fnr in range(FirstFile[model_number] + rank, LastFile[model_number]+1, size): # Divide up the input files across the processors.

        ## These are instantaneous values for a single file. ##
        w_gal = []
        mass_gal = []
        photons_HI_gal = []
        photons_HI_tot_gal = [] 
        SFR_gal = []
        sSFR_gal = []
        metallicity_gal = []
        metallicity_tremonti_gal = []
        halo_part_count = []

        mass_central = []
        photons_HI_central = []
        photons_HI_tot_central = []

        Mfilt_gnedin = []
        Mfilt_sobacchi = []

        lyman_alpha = []

        L_UV = []
        M_UV = []

        M_UV_Obs = []
        mean_A = []

        ejected_fraction = []

        halo_count = []
        fesc_local = []

        ##  

        GG, Gal_Desc = ReadScripts.ReadGals_SAGE_DelayedSN(galaxies_filepath_array[model_number], fnr, number_snapshots[model_number], comm) # Read in the correct galaxy file. 

        G_Merged, _ = ReadScripts.ReadGals_SAGE_DelayedSN(merged_galaxies_filepath_array[model_number], fnr, number_snapshots[model_number], comm) # Also need the merged galaxies.
        G = ReadScripts.Join_Arrays(GG, G_Merged, Gal_Desc) # Then join them together for all galaxies that existed at this Redshift. 

        ## These three arrays are used to control the escape fraction that depends upon quasar activity. ##
        ## A galaxy has the base escape fraction (specified by fesc_normalization[model_number][0]) but after a quasar event the escape fraction changes to fesc_normalization[model_number][1]. ##
        ## This boosted escape fraction lasts a dynamical time. ##
 
        QuasarActivityToggle = np.zeros((len(G))) # Array to specify whether the galaxy should have the base or boosted escape fraction.
        QuasarActivitySubstep = np.full((len(G)), -1) # Array to specify in which substep the quasar boosting begins. 
        QuasarSnapshot = np.full((len(G)), -1) # Array to keep track of when the quasar went off. Initialize to -1. 
        TargetQuasarTime = np.zeros((len(G))) # Array to keep track of the amount of time we want to keep the boosted escape fraction for. This will be the dynamical time WHEN THE QUASAR GOES OFF.
        QuasarBoostActiveTime = np.zeros((len(G))) # Array that tracks how long it has been since the quasar boosting was turned on.
        QuasarFractionalPhoton = np.full((len(G)), 0.0) # Array to keep track of what fraction of photons we should boost with the prescription.             

        keep_files = 1 # Flips to 0 when we are done with this file.
        current_model_number = model_number # Used to differentiate between outer model_number and the inner model_number because we can keep files across model_numbers.

        while(keep_files == 1):

            current_fesc_lyman_alpha = fesc_lyman_alpha[current_model_number] # Just reduces a bit of extra clutter down the road.
            current_source_efficiency = source_efficiency[current_model_number]
            current_halo_cut = halo_cut[current_model_number]
            N_dyntime = N_dynamicaltime[current_model_number]
            NumSubsteps = number_substeps[current_model_number]

            for snapshot_idx in range(0, len(SnapList[current_model_number])):
            
                current_snap = SnapList[current_model_number][snapshot_idx]

                w_gal = np.where((G.GridHistory[:, current_snap] != -1) & (G.GridStellarMass[:, current_snap] > 0.0) & (G.GridStellarMass[:, current_snap] < 1e5) & (G.GridCentralGalaxyMass[:, current_snap] >= m_low_SAGE) & (G.GridCentralGalaxyMass[:, current_snap] <=  m_high_SAGE) & (G.LenHistory[:, current_snap] > current_halo_cut) & (G.GridSFR[:, current_snap] >= 0.0))[0] # Only include those galaxies that existed at the current snapshot, had positive (but not infinite) stellar/Halo mass and Star formation rate.
          
                print("There were {0} galaxies for snapshot {1} (Redshift {2:.3f}) model {3}.".format(len(w_gal), current_snap, AllVars.SnapZ[current_snap], current_model_number))
               
                if (len(w_gal) == 0): 
                    continue

                number_galaxies[current_model_number][snapshot_idx] += len(w_gal)
                halo_count = G.LenHistory[w_gal, current_snap]
                mass_gal = np.log10(G.GridStellarMass[w_gal, current_snap] * 1.0e10 / AllVars.Hubble_h)
                SFR_gal = np.log10(G.GridSFR[w_gal, current_snap]) # Msun yr^-1.  Log Units.        
           
                halo_part_count = G.LenHistory[w_gal, current_snap]
                metallicity_gal = G.GridZ[w_gal, current_snap]  
                metallicity_tremonti_gal = np.log10(G.GridZ[w_gal, current_snap] / 0.02) + 9.0
                mass_central = np.log10(G.GridCentralGalaxyMass[w_gal, current_snap] * 1.0e10 / AllVars.Hubble_h)  
                Mfilt_gnedin = G.MfiltGnedin[w_gal, current_snap]
                Mfilt_sobacchi = G.MfiltSobacchi[w_gal, current_snap]
                ejected_fraction = G.EjectedFraction[w_gal, current_snap]
            
                L_UV = SFR_gal + 39.927 # Using relationship from STARBURST99, units of erg s^-1 A^-1. Log Units.
                M_UV = AllVars.Luminosity_to_ABMag(L_UV, 1600)

                if (calculate_observed_LF == 1): # Calculate the UV extinction if requested. 
                    M_UV_obs = calculate_UV_extinction(AllVars.SnapZ[current_snap], L_UV, M_UV[snap_idx])
            
                #print("fesc_prescription[current_model_number] = {0} \t fesc_normalization[current_model_number] = {1}".format(fesc_prescription[current_model_number], fesc_normalization[current_model_number]))
                old = 0

                ## Here we do calculations involving the quasars temporarily boosting the escape fraction. ##            
                if (fesc_prescription[current_model_number] == 3): # Only do it if this model is using this escape fraction prescription.

                    if (old == 1):
                        N_quasars_tmp = 0.0
                        N_quasars_boost_tmp = 0.0
                        dynamicaltime_quasars_tmp = [] 
                        dynamicaltime_all_tmp = [] 
                        for gal_idx in w_gal:
                            dynamicaltime_all_tmp.append(G.DynamicalTime[gal_idx, current_snap])

                            ## Check to see if it's been N dynamical times since the quasar went off. ##
                            if (QuasarActivityToggle[gal_idx] == 1):
                                N_quasars_boost_tmp += 1 
                                quasar_snap = QuasarSnapshot[gal_idx]
                                assert(quasar_snap != -1)  
                                dt = (float(AllVars.Lookback_Time[quasar_snap]) - float(AllVars.Lookback_Time[SnapList[model_number][snapshot_idx]])) * 1.0e3 # Time since the last quasar went off.                                              
                                assert(dt >= 0.0) # Just make sure we're doing the lookback time difference correctly.
                                if (dt >= TargetQuasarTime[gal_idx] * N_dyntime): # If it has been longer than N dynamical times, turn the boosted escape fraction off.
                                    QuasarActivityToggle[gal_idx] = 0
                                    QuasarSnapshot[gal_idx] = -1
                                    TargetQuasarTime[gal_idx] = 0.0

                            ## Check to see if any quasars went off during this snapshot. ##       
                            if (G.QuasarActivity[gal_idx, current_snap] == 1): # If a quasar went off...
                                QuasarActivityToggle[gal_idx] = 1 # Turn on the boosted escape fraction.
                                QuasarSnapshot[gal_idx] = current_snap # Remember the snapshot that this quasar went off.
                                TargetQuasarTime[gal_idx] = G.DynamicalTime[gal_idx, current_snap] # Keep track of how long we want the boosted escape fraction to last for.
                                
                                N_quasars_tmp += 1
                                dynamicaltime_quasars_tmp.append(G.DynamicalTime[gal_idx, current_snap]) 

                    else: 
                        N_quasars_tmp = 0.0
                        N_quasars_boost_tmp = 0.0
                        dynamicaltime_quasars_tmp = [] 
                        dynamicaltime_all_tmp = []
                        dynamicaltime_all_tmp.extend(G.DynamicalTime[w_gal, current_snap])

                        ## Check to see if any quasars went off during this snapshot. ## 
                        w_quasar_activity = w_gal[np.where((G.QuasarActivity[w_gal, current_snap] == 1))[0]] 
                        if (len(w_quasar_activity) > 0): # If a quasar went off...
                            N_quasars_tmp += len(w_quasar_activity)                   
                            QuasarActivityToggle[w_quasar_activity] = 1 # Turn on the boosted escape fraction.
                            QuasarSnapshot[w_quasar_activity] = current_snap # Remember the snapshot that this quasar went off.
                            TargetQuasarTime[w_quasar_activity] = G.DynamicalTime[w_quasar_activity, current_snap] # Keep track of how long we want the boosted escape fraction to last for.
                            #print(G.QuasarSubstep[w_quasar_activity, current_snap])
                            QuasarActivitySubstep[w_quasar_activity] = G.QuasarSubstep[w_quasar_activity, current_snap] # The substep of the snapshot the quasar went off. 
                            dynamicaltime_quasars_tmp.extend(G.DynamicalTime[w_quasar_activity, current_snap])

                        ## Update the boost times and see if we need to turn anything off ## 
                        w_active_boost = w_gal[np.where((QuasarActivityToggle[w_gal] > 0.0))[0]] # Use greater than because comparing a float.
                        if (len(w_active_boost) > 0):
                            quasar_snap = QuasarSnapshot[w_active_boost]
                            N_quasars_boost_tmp += len(w_active_boost)
                            
                            assert(quasar_snap.any() != -1) # Just a sanity check.

                            dt = ((AllVars.Lookback_Time[SnapList[current_model_number][snapshot_idx] - 1]) - (AllVars.Lookback_Time[SnapList[current_model_number][snapshot_idx]])) * 1.0e3 # Time between the previous snapshot and now (in Myr). 
                          
                            w_just_turned_on = w_active_boost[np.where((quasar_snap == current_snap))[0]] # In this case, the boosted escape fraction turned on during this snapshot.  As a result, we need to account for it turning on during a substep and only providing a fractional boost.
                            w_not_just_turned_on = w_active_boost[np.where((quasar_snap != current_snap))[0]] 

                            QuasarBoostActiveTime[w_just_turned_on] += dt * (NumSubsteps - QuasarActivitySubstep[w_just_turned_on]) / NumSubsteps # This adds the time spanned by this snapshot weighted by when substep in which the quasar went off. E.g., if there are 10 substeps and the quasar went off during substep 3, then the boosted escape fraction occurs for (10 - 3) / 10 = 0.7 of the time.  
                            QuasarFractionalPhoton[w_just_turned_on] = (NumSubsteps - QuasarActivitySubstep[w_just_turned_on]) / NumSubsteps # Only boost the photons for a fraction of the time during the snapshot step. 

                            QuasarBoostActiveTime[w_not_just_turned_on] += dt # The boosted escape fraction lasts for the entire snapshot.
                            QuasarFractionalPhoton[w_not_just_turned_on] = 1.0
                    
                            w_shutoff = w_active_boost[np.greater(QuasarBoostActiveTime[w_active_boost], TargetQuasarTime[w_active_boost] * N_dyntime).nonzero()] # If it's been longer than N dynamical times, turn the boosted escape fraction off.
# np.greater returns True if it's time to turn off, .nonzero() returns the indices of those, and then we want the index within the GLOBAL galaxy array so wrap it all within w_active_boost.

                            # Now for those galaxies for which we wish to shut off their boosted escape fraction, we need to determine on which substep we are turning it off.
                            # This will then give us a fraction of the photons emitted that we should boost for this final time.

                            TimeIntoSnapshot = TargetQuasarTime[w_shutoff] - (QuasarBoostActiveTime[w_shutoff] - dt) # This tells us how much time into this snapshot this quasar will be turned on for.
                            FractionIntoSnapshot = TimeIntoSnapshot / dt # Then this is the fraction of time the quasar will be on for.
                            QuasarFractionalPhoton[w_shutoff] = FractionIntoSnapshot # Weight the boosted escape fraction to be correct. E.g., if a galaxy is turned off 40% into the snapshot, then the escape fraction should only be boosted for 40% of the time.
                            #TimePerSubstep = dt / NumSubsteps
                            #TimeOverTarget = np.subtract(QuasarBoostActiveTime[w_active_boost], TargetQuasarTime[w_active_boost]) # For galaxies which have not reached their target time, this number will be negative.
                            #StepsOverTarget = int(TimeOverTarget / TimePerSubstep) # Cast as an integer to match the definition of the initial substep boosting.
                           

#                            w_gal_over_target = w_active_boost[np.where((StepsOverTarget > 0.0))[0]] # Gets those galaxies who have gone over target time.
#                            print("w_gal_over_target {0}".format(w_gal_over_target))
#                            print("NumSubsteps {0}".format(NumSubsteps))
#                            if (len(w_gal_over_target) > 0):
#                                QuasarFractionalPhoton[w_gal_over_target] = StepsOverTarget[StepsOverTarget > 0.0] / NumSubsteps # Sets the fraction of time boosted to be correct. E.g., if a galaxy is 4 substeps over the target 
                            
                            QuasarActivityToggle[w_shutoff] = 0 # Then turn off all the toggles.
                            QuasarSnapshot[w_shutoff] = -1
                            TargetQuasarTime[w_shutoff] = 0.0
                            QuasarBoostActiveTime[w_shutoff] = 0.0
                            QuasarActivitySubstep[w_shutoff] = -1 
                     
                    '''
                    if (N_quasars_tmp > 0):
                        print("QuasarSnapshot {0}".format(QuasarSnapshot[QuasarSnapshot != -1]))
                        print("QuasarActivityToggle {0}".format(QuasarActivityToggle[QuasarActivityToggle != 0]))
                        print("TargetQuasarTime {0}".format(TargetQuasarTime[TargetQuasarTime != 0]))
                        exit()
                    '''
                    if (N_quasars_tmp > 0):           
                        (dynamicaltime_quasars_mean_z[current_model_number][snapshot_idx], dynamicaltime_quasars_std_z[current_model_number][snapshot_idx]) = update_cumulative_stats(dynamicaltime_quasars_mean_z[current_model_number][snapshot_idx], dynamicaltime_quasars_std_z[current_model_number][snapshot_idx], N_quasars_z[current_model_number][snapshot_idx], np.mean(dynamicaltime_quasars_tmp), np.std(dynamicaltime_quasars_tmp), N_quasars_tmp)               
                    (dynamicaltime_all_mean_z[current_model_number][snapshot_idx], dynamicaltime_all_std_z[current_model_number][snapshot_idx]) = update_cumulative_stats(dynamicaltime_all_mean_z[current_model_number][snapshot_idx], dynamicaltime_all_std_z[current_model_number][snapshot_idx], N_z[current_model_number][snapshot_idx], np.mean(dynamicaltime_all_tmp), np.std(dynamicaltime_all_tmp), len(w_gal))
                                
                    N_quasars_z[current_model_number][snapshot_idx] += N_quasars_tmp
                    N_quasars_boost_z[current_model_number][snapshot_idx] += N_quasars_boost_tmp
  
                galaxy_halo_mass_mean_local, galaxy_halo_mass_std_local = Calculate_HaloPartStellarMass(halo_part_count, mass_gal, stellar_mass_halolen_lower[current_model_number], stellar_mass_halolen_upper[current_model_number])
                galaxy_halo_mass_mean[current_model_number][snapshot_idx] += pow(10, galaxy_halo_mass_mean_local) / (LastFile[current_model_number] + 1) # Adds to the average of the mean. 

                #plot_photon_quasar_fraction(SnapList[current_model_number][snapshot_idx], fnr, save_tags[current_model_number], QuasarFractionalPhoton[w_gal], QuasarActivityToggle[w_gal], NumSubsteps)

                fesc_local = calculate_fesc(fesc_prescription[current_model_number], fesc_normalization[current_model_number], mass_central, ejected_fraction, QuasarActivityToggle[w_gal], QuasarFractionalPhoton[w_gal], old) 
                     
                photons_HI_gal = calculate_photons(SFR_gal, metallicity_gal)    
                photons_HI_gal_nonlog = [10**x for x in photons_HI_gal]
                ionizing_photons = np.multiply(photons_HI_gal_nonlog, fesc_local)

                ## We have now calculated all the base properties for galaxies within this snapshot.  Calculate the relevant statistics and put them into their arrays. ##

                ### Functions of Galaxies/Stellar Mass ###
                (counts_local, bin_edges, bin_middle) = AllVars.Calculate_Histogram(mass_gal, bin_width, 0, m_gal_low, m_gal_high) # Bin the Stellar Mass 
                SMF[current_model_number][snapshot_idx] += counts_local 

                (mean_fesc_galaxy_local, std_fesc_galaxy_local, N_local, bin_middle) = AllVars.Calculate_2D_Mean(mass_gal, fesc_local, bin_width, m_gal_low, m_gal_high)

                (mean_fesc_galaxy_array[current_model_number][snapshot_idx], std_fesc_galaxy_array[current_model_number][snapshot_idx]) = update_cumulative_stats(mean_fesc_galaxy_array[current_model_number][snapshot_idx], std_fesc_galaxy_array[current_model_number][snapshot_idx], N_galaxy_array[current_model_number][snapshot_idx], mean_fesc_galaxy_local, std_fesc_galaxy_local, N_local)
#
                N_galaxy_array[current_model_number][snapshot_idx] += N_local 

                ### Functions of Halos/Halo Mass ###
                (mean_ejected_halo_local, std_ejected_halo_local, N_local, bin_middle) = AllVars.Calculate_2D_Mean(mass_central, ejected_fraction, bin_width, m_low, m_high) 
                (mean_ejected_halo_array[current_model_number][snapshot_idx], std_ejected_halo_array[current_model_number][snapshot_idx]) = update_cumulative_stats(mean_ejected_halo_array[current_model_number][snapshot_idx], std_ejected_halo_array[current_model_number][snapshot_idx], N_halo_array[current_model_number][snapshot_idx], mean_ejected_halo_local, std_ejected_halo_local, N_local) 
                
                (mean_fesc_halo_local, std_fesc_halo_local, N_local, bin_middle) = AllVars.Calculate_2D_Mean(mass_central, fesc_local, bin_width, m_low, m_high) 
                (mean_fesc_halo_array[current_model_number][snapshot_idx], std_fesc_halo_array[current_model_number][snapshot_idx]) = update_cumulative_stats(mean_fesc_halo_array[current_model_number][snapshot_idx], std_fesc_halo_array[current_model_number][snapshot_idx], N_halo_array[current_model_number][snapshot_idx], mean_fesc_halo_local, std_fesc_halo_local, N_local) 

                (mean_Ngamma_halo_local, std_Ngamma_halo_local, N_local, bin_middle) = AllVars.Calculate_2D_Mean(mass_central, ionizing_photons, bin_width, m_low, m_high) 

                mean_Ngamma_halo_local = np.divide(mean_Ngamma_halo_local, 1.0e50)
                std_Ngamma_halo_local = np.divide(std_Ngamma_halo_local, 1.0e50)

                (mean_Ngamma_halo_array[current_model_number][snapshot_idx], std_Ngamma_halo_array[current_model_number][snapshot_idx]) = update_cumulative_stats(mean_Ngamma_halo_array[current_model_number][snapshot_idx], std_Ngamma_halo_array[current_model_number][snapshot_idx], N_halo_array[current_model_number][snapshot_idx], mean_Ngamma_halo_local, std_Ngamma_halo_local, N_local) 
                
                N_halo_array[current_model_number][snapshot_idx] += N_local 

                ### Functions of redshift ###
                sum_Ngamma_z_array[current_model_number][snapshot_idx] += np.sum(np.divide(ionizing_photons, 1.0e50))   

                (mean_fesc_z_array[current_model_number][snapshot_idx], std_fesc_z_array[current_model_number][snapshot_idx]) = update_cumulative_stats(mean_fesc_z_array[current_model_number][snapshot_idx], std_fesc_z_array[current_model_number][snapshot_idx], N_z[current_model_number][snapshot_idx], np.mean(fesc_local), np.std(fesc_local), len(w_gal))
     
                N_z[current_model_number][snapshot_idx] += len(w_gal)


            done_model[current_model_number] = 1
            if (current_model_number < number_models):                
                keep_files =  same_files[current_model_number] # If we want to keep the files and bin the next model, this will do it.
                current_model_number += 1 # Update the inner loop model number.
for snap in SnapList[0]:
    print("Mean Dynamical Time for quasars at redshift {0} is {1} with std {2}".format(AllVars.SnapZ[snap], dynamicaltime_quasars_mean_z[0][snap], dynamicaltime_quasars_std_z[0][snap]))

exit()
#StellarMassFunction(SnapList, SMF, simulation_norm, FirstFile, LastFile, NumFile, galaxy_halo_mass_mean, model_tags, 1, "Kali_Tiamat_SMF") ## PARALLEL COMPATIBLE
#plot_ejectedfraction(SnapList, mean_ejected_halo_array, std_ejected_halo_array, N_halo_array, model_tags, "tiamat_newDelayedComp_ejectedfract_highz") ## PARALELL COMPATIBLE # Ejected fraction as a function of Halo Mass 
#plot_fesc(SnapList, mean_fesc_z_array, std_fesc_z_array, N_z, model_tags, "Quasarfesc_z_DynamicalTimes") ## PARALELL COMPATIBLE 
plot_quasars_count(SnapList, N_quasars_z, N_quasars_boost_z, N_z, fesc_prescription, simulation_norm, FirstFile, LastFile, NumFile, model_tags, "Kali_Quasar_Newest")
plot_fesc_galaxy(SnapList, PlotSnapList, simulation_norm, mean_fesc_galaxy_array, std_fesc_galaxy_array, N_galaxy_array, mean_fesc_halo_array, std_fesc_halo_array,  N_halo_array, galaxy_halo_mass_mean, model_tags, "fesc_Kali_Quasar_Newest", 1, save_tags)
#plot_photoncount(SnapList, sum_Ngamma_z_array, simulation_norm, FirstFile, LastFile, NumFile, model_tags, "Ngamma_KaliTiamat_Quasar") ## PARALELL COMPATIBLE
#plot_mvir_Ngamma(SnapList, mean_Ngamma_halo_array, std_Ngamma_halo_array, N_halo_array, model_tags, "Mvir_Ngamma_test", fesc_prescription, fesc_normalization, "/lustre/projects/p004_swin/jseiler/tiamat/halo_ngamma/") ## PARALELL COMPATIBLE 

