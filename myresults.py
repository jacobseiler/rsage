#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')

import os
import heapq
import h5py as h5
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

comm= MPI.COMM_WORLD
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

m_low = 8.5 # We only sum the photons coming from halos within the mass range m_low < Halo Mass < m_high
m_high = 14

m_gal_low = 4
m_gal_high = 14

m_low_SAGE = pow(10, m_low)/1.0e10 * AllVars.Hubble_h
m_high_SAGE = pow(10, m_high)/1.0e10 * AllVars.Hubble_h

bin_width = 0.1
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
    for i in xrange(0, len(array)):
        total *= array[i]
    return total

##

def Calculate_2D_Mean(data_x, data_y, bin_width, min_hist_x = None, max_hist_x = None):     

    '''
    Calculates the mean of the y-data that lies within binned x-data.   

    Parameters
    ----------
    data_x : `numpy.darray'
    Data that will be binned and provide the bins for the y-mean.
    data_y : `numpy.darray'
    Data that will be averaged in each of the bins defined by the x-data.
    bin_width : float
    Width of each x-bin.
    min_hist_x, max_hist_x: float (optional)
    Defines the x-bounds that we will be binning over.
    If no values defined, range will be the minimum/maximum data point +/- 10 times the bin_width.

    Returns
    -------
    mean_data_y, std_data_y : `numpy.darray'    
    Arrays that contain the mean and standard deviation for the y-data as binned by the x-data.
    N_data_y : `numpy.darray'
    Array that contains the number of data points in each x-bin.    
    bins_mid : `numpy.darray'
    The mid-point coordinate for the x-bins. 

    Units
    -----
    All units are kept the same as the inputs.
    '''

    if not np.isfinite(min_hist_x):
        raise ValueError("xmin should be finite")

    if not np.isfinite(max_hist_x):
        raise ValueError("xmax should be finite")

    if (min_hist_x == None): 
        range_low = np.floor(min(data_x)) - 10*bin_width
        range_high = np.floor(max(data_x)) + 10*bin_width
    else:
        range_low = min_hist_x 
        range_high = max_hist_x 

    if range_high <= range_low: 
        raise ValueError("The upper bin range should be less than the lower bin range")
 
    NB = round((range_high - range_low) / bin_width) 

    bins = np.arange(range_low, range_high + bin_width, bin_width)
    bins_mid = bins + bin_width/2.0
    bins_mid = bins_mid[:-1] # len(bins_mid) should be 1 less than len(bins) as the last bin doesn't have a midpoint.   

    bins_x = np.digitize(data_x, bins) # Provides the indices for which bin each x-data point belongs to.
    N_data_y = np.zeros((len(bins)))  
 
    data_y_binned = []
    for i in xrange(0, len(bins)): 
        data_y_binned.append([])
    for i in xrange(0, len(data_x)):
        idx = bins_x[i]
        if idx == len(data_y_binned): # Fixes up binning index edge case.
            idx -= 1

        data_y_binned[idx].append(data_y[i])
        N_data_y[idx] += 1 

    mean_data_y = np.zeros((len(bins))) 
    std_data_y = np.zeros((len(bins))) 

    for i in xrange(0, len(bins)):
        if len(data_y_binned[i]) != 0: # If there was any y-data placed into this bin. 
            mean_data_y[i] = np.mean(data_y_binned[i])
            std_data_y[i] = np.std(data_y_binned[i])
            
        else: # Otherwise if there were no data points, simply fill it with a nan.
            mean_data_y[i] = 0.0 
            std_data_y[i] = 0.0

    return mean_data_y[:-1], std_data_y[:-1], N_data_y[:-1], bins_mid

##

def Sum_Log(array):
    '''
    Performs an element wise sum of an array who's elements are in log-space.

    Parameters
    ----------
    array : `array-like'
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
    for i in xrange(0, len(array)):
        sum_total += 10**array[i]

    return sum_total

##

def Std_Log(array, mean):
    '''
    Calculates the standard deviation of an array with elements in log-space. 

    Parameters
    ----------
    array : `array-like'
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
    for i in xrange(0, len(array)):
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

    Parameters
    ----------
    mean_pool, std_pool : `array-like' of floats.
    Arrays that contain the current pooled means/standard deviation (for rank 0) or just a junk input (for other ranks).
    mean_local, mean_std : float.
    The non-pooled mean and standard deviation unique for each process.
    N_local : int 
    Number of data points used to calculate the mean/standard deviation that is going to be added to the pool.

    Returns
    -------
    mean_pool, std_pool : `array-like' of floats.
    Original array with the new pooled mean/standard deviation appended (for rank 0) or the new pooled mean/standard deviation only (for other ranks).

    Units
    -----
    All units are the same as the input.
    All inputs are in real-space (not log-space).
    '''

    if isinstance(mean_local, list) == True:    
        if len(mean_local) != len(std_local):
            print "len(mean_local) = %d \t len(std_local) = %d" %(len(mean_local), len(std_local))
            raise ValueError("Lengths of mean_local and std_local should be equal")
   
    if ((type(mean_local).__module__ == np.__name__) == True or (isinstance(mean_local, list) == True)): # Checks to see if we are dealing with arrays. 
    
        N_times_mean_local = np.multiply(N_local, mean_local)
        N_times_var_local = np.multiply(N_local, np.multiply(std_local, std_local))
        
        N_local = np.array(N_local).astype(int)
        N_times_mean_local = np.array(N_times_mean_local).astype(np.float32)

        if rank == 0: # Only rank 0 holds the final arrays so only it requires proper definitions.
            N_times_mean_pool = np.zeros_like(N_times_mean_local) 
            N_pool = np.zeros_like(N_local)
            N_times_var_pool = np.zeros_like(N_times_var_local)
        else:
            N_times_mean_pool = None
            N_pool = None
            N_times_var_pool = None

        comm.Barrier()
        comm.Reduce([N_times_mean_local, MPI.FLOAT], [N_times_mean_pool, MPI.FLOAT], op = MPI.SUM, root = 0) # Sum the arrays across processors.
        comm.Reduce([N_local, MPI.INT],[N_pool, MPI.INT], op = MPI.SUM, root = 0)   
        comm.Reduce([N_times_var_local, MPI.FLOAT], [N_times_var_pool, MPI.FLOAT], op = MPI.SUM, root = 0)

    else:
    
        N_times_mean_local = N_local * mean_local
        N_times_var_local = N_local * std_local * std_local

        N_times_mean_pool = comm.reduce(N_times_mean_local, op = MPI.SUM, root = 0)
        N_pool = comm.reduce(N_local, op = MPI.SUM, root = 0)
        N_times_var_pool = comm.reduce(N_times_var_local, op = MPI.SUM, root = 0)
    
    if rank == 0:

        mean_pool_function = np.zeros((len(N_pool)))
        std_pool_function = np.zeros((len(N_pool)))

        for i in xrange(0, len(N_pool)):
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


def StellarMassFunction(SnapList, SMF, simulation_norm, FirstFile, LastFile, NumFile, model_tags, observations, output_tag):
    '''
    Calculates the stellar mass function for given galaxies with the option to overplot observations by Song et al. (2013) at z = 6, 7, 8 and/or Baldry et al. (2008) at z = 0.1. 
    Parallel compatible.
    NOTE: The plotting assumes the redshifts we are plotting at are the same for each model. 

    Parameters
    ---------
    SnapList : Nested 'array-like`, SnapList[model_number0] = [snapshot0_model0, ..., snapshotN_model0], with length equal to the number of models.
        Snapshots that we plot the stellar mass function at for each model.
    SMF : Nested 2-dimensional `array-like', SMF[model_number0][snapshot0]  = [bin0galaxies, ..., binNgalaxies], with length equal to the number of bins (NB_gal). 
        The count of galaxies within each stellar mass bin.  Bounds are given by 'm_gal_low' and 'm_gal_high' in bins given by 'bin_width'. 
    simulation_norm : `array-like' with length equal to the number of models.
        Denotes which simulation each model uses.  
        0 : MySim
        1 : Mini-Millennium
        2 : Tiamat (down to z = 5)
        3 : Extended Tiamat (down to z = 1.6ish).
    FirstFile, LastFile, NumFile : `array-like' of integers with length equal to the number of models.
        The file numbers for each model that were read in (defined by the range between [FirstFile, LastFile] inclusive) and the TOTAL number of files for this model (we may only be plotting a subset of the volume). 
    model_tags : `array-like' of strings with length equal to the number of models.
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

    for model_number in xrange(0, len(SnapList)):
        counts_array.append([])
        bin_middle_array.append([])
        redshift_labels.append([])

    ####
    for model_number in xrange(0, len(SnapList)): # Does this for each of the models. 

        ## Normalization for each model. ##
        if (simulation_norm[model_number] == 0):
            AllVars.Set_Params_Mysim()
        elif (simulation_norm[model_number] == 1):
            AllVars.Set_Params_MiniMill()
        elif (simulation_norm[model_number] == 2):
            AllVars.Set_Params_Tiamat()
        elif (simulation_norm[model_number] == 3):
            AllVars.Set_Params_Tiamat_extended()

        
        box_factor = (LastFile[model_number] - FirstFile[model_number] + 1.0)/(NumFile[model_number]) # This factor allows us to take a sub-volume of the box and scale the results to represent the entire box.
        print "We are creating the stellar mass function using %.4f of the box's volume." %(box_factor)
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * bin_width * box_factor 
        normalization_array.append(norm)
        print "norm = %.4f" %(norm)
        ####
    
        for snapshot_idx in xrange(0, len(SnapList[model_number])): # Loops for each snapshot in each model.

            tmp = 'z = %.2f' %(AllVars.SnapZ[SnapList[model_number][snapshot_idx]]) # Assigns a redshift label.
            redshift_labels[model_number].append(tmp)

            ## We perform the plotting on Rank 0 so only this rank requires the final counts array. ##
            if rank == 0:
                counts_total = np.zeros_like(SMF[model_number][snapshot_idx])
            else:
                counts_total = None

            comm.Reduce([SMF[model_number][snapshot_idx], MPI.DOUBLE], [counts_total, MPI.DOUBLE], op = MPI.SUM, root = 0) # Sum all the stellar mass and pass to Rank 0.

            if rank == 0:
                counts_array[model_number].append(counts_total)
    #           print counts_array[model_number][snapshot_idx]
                bin_middle_array[model_number].append(np.arange(m_gal_low, m_gal_high+bin_width, bin_width)[:-1] + bin_width * 0.5)
            ####
    

    ## Plotting ##

    if rank == 0: # Plot only on rank 0.
        f = plt.figure()  
        ax = plt.subplot(111)  

        for model_number in xrange(0, len(SnapList)):
            print model_tags[model_number]
            for snapshot_idx in xrange(0, len(SnapList[model_number])):
                if model_number == 0: # We assume the redshifts for each model are the same, we only want to put a legend label for each redshift once.
                    title = redshift_labels[model_number][snapshot_idx]
                else:
                    title = ''
                plt.plot(bin_middle_array[model_number][snapshot_idx], counts_array[model_number][snapshot_idx] / normalization_array[model_number], color = PlotScripts.colors[snapshot_idx], linestyle = PlotScripts.linestyles[model_number], rasterized = True, label = title, linewidth = PlotScripts.global_linewidth) 

                #print counts_array[model_number][snapshot_idx] / normalization_array[model_number]
 
        for model_number in xrange(0, len(SnapList)): # Place legend labels for each of the models. NOTE: Placed after previous loop for proper formatting of labels. 
            plt.plot(1e100, 1e100, color = 'k', linestyle = PlotScripts.linestyles[model_number], label = model_tags[model_number], rasterized=True, linewidth = PlotScripts.global_linewidth)
    
        ## Adjusting axis labels/limits. ##

        plt.yscale('log', nonposy='clip')
        plt.axis([6, 11.5, 1e-6, 1e-1])

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
            plt.errorbar(Song_z6[:,0], 10**Song_z6[:,1], yerr= (10**Song_z6[:,1] - 10**Song_z6[:,3], 10**Song_z6[:,2] - 10**Song_z6[:,1]), xerr = 0.25, capsize = caps, elinewidth = PlotScripts.global_errorwidth, alpha = 1.0, lw=2.0, marker='o', ls='none', label = 'Song 2015, z = 6', color = '#bd0026', rasterized=True)
            plt.errorbar(Song_z7[:,0], 10**Song_z7[:,1], yerr= (10**Song_z7[:,1] - 10**Song_z7[:,3], 10**Song_z7[:,2] - 10**Song_z7[:,1]), xerr = 0.25, capsize = caps, alpha=0.75, elinewidth = PlotScripts.global_errorwidth, lw=1.0, marker='o', ls='none', label = 'Song 2015, z = 7', color = 'blue', rasterized=True)
            plt.errorbar(Song_z8[:,0], 10**Song_z8[:,1], yerr= (10**Song_z8[:,1] - 10**Song_z8[:,3], 10**Song_z8[:,2] - 10**Song_z8[:,1]), xerr = 0.25, capsize = caps, alpha=0.75, elinewidth = PlotScripts.global_errorwidth, lw=1.0, marker='o', ls='none', label = 'Song 2015, z = 8', color = 'green', rasterized=True)
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

        leg = plt.legend(loc='upper right', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize(PlotScripts.global_legendsize)

        outputFile = './%s%s' %(output_tag, output_format) 
        plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

##

def plot_fesc(SnapList, mean_z_fesc, std_z_fesc, N_fesc, model_tags, output_tag):
    '''
    Plots the escape fraction as a function of redshift for the given galaxies. 
    Parallel compatible.
    Accepts 3D arrays of the escape fraction binned into Halo Mass bins to plot the escape fraction for multiple models. 

    Parameters
    ---------
    SnapList : Nested `array-like', SnapList[model_number0] = [snapshot0_model0, ..., snapshotN_model0], with length equal to the number of models.
    Snapshots for each model. 
    mean_z_fesc, std_z_fesc, N_fesc : Nested 2-dimensional `array-like', mean_z_fesc[model_number0][snapshot0]  = [z0_meanfesc, ..., zN_meanfesc], with length equal to the number of models. 
    Mean/Standard deviation for fesc at each redshift. N_ejected is the number of data points in each bin. 
    model_tags : `array-like' of strings with length equal to the number of models.
    Strings that contain the tag for each model.  Will be placed on the plot.
    output_tag : string
    Name of the file that will be generated.

    Returns
    -------
    No returns.
    Generates and saves the plot (named via output_tag).   
    '''

    print "Plotting fesc as a function of redshift."

    ## Array initialization ##
    pooled_mean_fesc = []
    pooled_std_fesc = []

    for model_number in xrange(0, len(SnapList)): # Loop for each model. 

        pooled_mean_fesc, pooled_std_fesc = calculate_pooled_stats(pooled_mean_fesc, pooled_std_fesc, mean_z_fesc[model_number], std_z_fesc[model_number], N_fesc[model_number]) # Calculates the pooled mean/standard deviation for this snapshot.  Only rank 0 receives a proper value here; the other ranks don't need this information. 
    
    if (rank == 0):
        ax1 = plt.subplot(111)

    for model_number in xrange(0, len(SnapList)):
    
        ## Calculate lookback time for each snapshot ##
        t = np.empty(len(SnapList[model_number]))
        for snapshot_idx in xrange(0, len(SnapList[model_number])):  
            t[snapshot_idx] = (t_BigBang - cosmo.lookback_time(AllVars.SnapZ[SnapList[model_number][snapshot_idx]]).value) * 1.0e3   
                
        mean = pooled_mean_fesc[model_number]
        std = pooled_std_fesc[model_number]   

        ax1.plot(t, mean, color = PlotScripts.colors[model_number], linestyle = PlotScripts.linestyles[model_number], label = model_tags[model_number], linewidth = PlotScripts.global_linewidth)  
        #ax1.fill_between(t, np.subtract(mean,std), np.add(mean,std), color = PlotScripts.colors[model_number], alpha = 0.25)

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

        ax1.set_ylim([0.0, 0.20])
        ax1.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = PlotScripts.global_labelsize) 
        ax1.set_ylabel(r'$f_\mathrm{esc}$', fontsize = PlotScripts.global_fontsize) 

        leg = ax1.legend(loc=1, numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize(PlotScripts.global_legendsize)

        plt.tight_layout()
        outputFile = './%s%s' %(output_tag, output_format)
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()

##

def plot_ejectedfraction(SnapList, mean_mvir_ejected, std_mvir_ejected, N_ejected, model_tags, output_tag): 
    '''
    Plots the ejected fraction as a function of the halo mass. 
    Parallel compatible.
    Accepts a 3D array of the ejected fraction so we can plot for multiple models and redshifts. 

    Parameters
    ---------
    SnapList : Nested `array-like', SnapList[model_number0] = [snapshot0_model0, ..., snapshotN_model0], with length equal to the number of models.
    Snapshots for each model. 
    mean_mvir_ejected, std_mvir_ejected, N_ejected : Nested 2-dimensional `array-like', mean_mvir_ejected[model_number0][snapshot0]  = [bin0_meanejected, ..., binN_meanejected], with length equal to the number of models. 
    Mean/Standard deviation for the escape fraction binned into Halo Mass bins. N_ejected is the number of data points in each bin. Bounds are given by 'm_low' and 'm_high' in bins given by 'bin_width'.   
    model_tags : `array-like' of strings with length equal to the number of models.
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

    print "Plotting the Ejected Fraction as a function of halo mass."

    ## Array initialization. ##
    title = []
    redshift_labels = []

    mean_ejected_array = []
    std_ejected_array = []

    mean_halomass_array = []
    std_halomass_array = []

    bin_middle_array = []

    for model_number in xrange(0, len(SnapList)):
        redshift_labels.append([])

        mean_ejected_array.append([])
        std_ejected_array.append([])

        mean_halomass_array.append([])
        std_halomass_array.append([])

        bin_middle_array.append([])
    
    bin_width = 0.1
 
    for model_number in xrange(0, len(SnapList)): 
        for snapshot_idx in xrange(0, len(SnapList[model_number])):
            print "Doing Snapshot %d" %(SnapList[model_number][snapshot_idx])
            tmp = 'z = %.2f' %(AllVars.SnapZ[SnapList[model_number][snapshot_idx]])
            redshift_labels[model_number].append(tmp)
            
            mean_ejected_array[model_number], std_ejected_array[model_number] = calculate_pooled_stats(mean_ejected_array[model_number], std_ejected_array[model_number], mean_mvir_ejected[model_number][snapshot_idx], std_mvir_ejected[model_number][snapshot_idx], N_ejected[model_number][snapshot_idx]) # Calculates the pooled mean/standard deviation for this snapshot.  Only rank 0 receives a proper value here; the other ranks don't need this information. 

            bin_middle_array[model_number].append(np.arange(m_low, m_high+bin_width, bin_width)[:-1] + bin_width * 0.5)
    
    if rank == 0:
        f = plt.figure()  
        ax1 = plt.subplot(111)  

        for model_number in xrange(0, len(SnapList)):
            for snapshot_idx in xrange(0, len(SnapList[model_number])): 
                ax1.plot(bin_middle_array[model_number][snapshot_idx], mean_ejected_array[model_number][snapshot_idx], color = PlotScripts.colors[snapshot_idx], linestyle = PlotScripts.linestyles[model_number], rasterized = True, label = redshift_labels[model_number][snapshot_idx], linewidth = PlotScripts.global_linewidth) 
                

        for model_number in xrange(0, len(SnapList)): # Just plot some garbage to get the legend labels correct.
            ax1.plot(np.nan, np.nan, color = 'k', linestyle = PlotScripts.linestyles[model_number], rasterized = True, label = model_tags[model_number], linewidth = PlotScripts.global_linewidth)

        ax1.set_xlabel(r'$\log_{10}\ M_{\mathrm{vir}}\ [M_{\odot}]$', size = PlotScripts.global_fontsize) 
        ax1.set_ylabel(r'$\mathrm{Ejected \: Fraction}$', size = PlotScripts.global_fontsize)
        ax1.set_xlim([9.5, 12])
        ax1.set_ylim([-0.05, 1.0])   

        ax1.xaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.025))
        
        leg = ax1.legend(loc=1, numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')

        outputFile = './%s%s' %(output_tag, output_format)
        plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
        print 'Saved file to', outputFile
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

    for model_number in xrange(0, len(SnapList)):
        redshift_labels.append([])

    mean_fesc_array.append([])
    std_fesc_array.append([])

    mean_halomass_array.append([])
    std_halomass_array.append([])

    bin_middle_array.append([])
    print "Plotting fesc against Mvir" 
    
    binwidth = 0.1
    Frequency = 1  
 
    for model_number in xrange(0, len(SnapList)): 
        for snapshot_idx in xrange(0, len(SnapList[model_number])):
            print "Doing Snapshot %d" %(SnapList[model_number][snapshot_idx])
            tmp = 'z = %.2f' %(AllVars.SnapZ[SnapList[model_number][snapshot_idx]])
            redshift_labels[model_number].append(tmp)

            minimum_mass = np.floor(min(mass_central[model_number][snapshot_idx])) - 10*binwidth
            maximum_mass = np.floor(max(mass_central[model_number][snapshot_idx])) + 10*binwidth

            minimum_mass = 6.0
            maximum_mass = 12.0

            binning_minimum = comm.allreduce(minimum_mass, op = MPI.MIN)
            binning_maximum = comm.allreduce(maximum_mass, op = MPI.MAX)
            
            halomass_nonlog = [10**x for x in mass_central[model_number][snapshot_idx]]
            (mean_fesc, std_fesc, N, bin_middle) = Calculate_2D_Mean(mass_central[model_number][snapshot_idx], fesc[model_number][snapshot_idx], binwidth, binning_minimum, binning_maximum)

            mean_fesc_array[model_number], std_fesc_array[model_number] = calculate_pooled_stats(mean_fesc_array[model_number], std_fesc_array[model_number], mean_fesc, std_fesc, N)
            mean_halomass_array[model_number], std_halomass_array[model_number] = calculate_pooled_stats(mean_halomass_array[model_number], std_halomass_array[model_number], np.mean(halomass_nonlog), np.std(halomass_nonlog), len(mass_central[model_number][snapshot_idx]))

            ## If want to do mean/etc of halo mass need to update script. ##
            bin_middle_array[model_number].append(bin_middle)
        
        mean_halomass_array[model_number] = np.log10(mean_halomass_array[model_number]) 
        
    if rank == 0:
        f = plt.figure()  
        ax1 = plt.subplot(111)  

        for model_number in xrange(0, len(SnapList)):
            for snapshot_idx in xrange(0, len(SnapList[model_number])):
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
#       for model_number in xrange(0, len(SnapList)):
#       ax1.plot(1e100, 1e100, color = 'k', ls = linestyles[model_number], label = model_tags[model_number], rasterized=True)
    
    
        leg = ax1.legend(loc='upper left', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')

        outputFile = './' + output_tag + output_format
        plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
        print 'Saved file to', outputFile
        plt.close()
##

def plot_mvir_Ngamma(SnapList, mean_mvir_Ngamma, std_mvir_Ngamma, N_Ngamma, model_tags, output_tag,fesc_prescription=None, fesc_normalization=None, fitpath=None): 
    '''
    Plots the number of ionizing photons (pure ngamma times fesc) as a function of halo mass. 
    Parallel compatible.
    The input data has been binned as a function of halo virial mass (Mvir), with the bins defined at the top of the file (m_low, m_high, bin_width). 
    Accepts 3D arrays to plot ngamma for multiple models. 

    Parameters
    ---------
    SnapList : Nested `array-like', SnapList[model_number0] = [snapshot0_model0, ..., snapshotN_model0], with length equal to the number of models.
    Snapshots for each model. 
    mean_mvir_Ngamma, std_mvir_Ngamma, N_Ngamma : Nested 2-dimensional `array-like', mean_mvir_Ngamma[model_number0][snapshot0]  = [bin0_meanNgamma, ..., binN_meanNgamma], with length equal to the number of bins. 
    Mean/Standard deviation/number of data points in each halo mass (Mvir) bin.
    The number of photons is in units of 1.0e50 s^-1.   
    model_tags : `array-like' of strings with length equal to the number of models.
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

    print "Plotting ngamma*fesc against the halo mass" 

    ## Array initialization. ##
    title = []
    redshift_labels = []

    mean_ngammafesc_array = []
    std_ngammafesc_array = []

    mean_halomass_array = []
    std_halomass_array = []

    bin_middle_array = []

    for model_number in xrange(0, len(SnapList)):
        redshift_labels.append([])

    mean_ngammafesc_array.append([])
    std_ngammafesc_array.append([])

    mean_halomass_array.append([])
    std_halomass_array.append([])

    bin_middle_array.append([])
    
    for model_number in xrange(0, len(SnapList)): 
        for snapshot_idx in xrange(0, len(SnapList[model_number])):
            print "Doing Snapshot %d" %(SnapList[model_number][snapshot_idx])
            tmp = 'z = %.2f' %(AllVars.SnapZ[SnapList[model_number][snapshot_idx]])
            redshift_labels[model_number].append(tmp)

            N = N_Ngamma[model_number][snapshot_idx]
            
            mean_ngammafesc_array[model_number], std_ngammafesc_array[model_number] = calculate_pooled_stats(mean_ngammafesc_array[model_number], std_ngammafesc_array[model_number], mean_mvir_Ngamma[model_number][snapshot_idx], std_mvir_Ngamma[model_number][snapshot_idx], N) # Collate the values from all processors.   
            bin_middle_array[model_number].append(np.arange(m_low, m_high+bin_width, bin_width)[:-1] + bin_width * 0.5) 
    
    if rank == 0:
        f = plt.figure()  
        ax1 = plt.subplot(111)  

        for model_number in xrange(0, len(SnapList)):
            count = 0
            for snapshot_idx in xrange(0, len(SnapList[model_number])):
                if model_number == 0:
                    title = redshift_labels[model_number][snapshot_idx]
                else:
                    title = ''

                mean = np.zeros((len(mean_ngammafesc_array[model_number][snapshot_idx])), dtype = np.float32)
                std = np.zeros((len(mean_ngammafesc_array[model_number][snapshot_idx])), dtype=np.float32)

                for i in xrange(0, len(mean)):
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
                    print "The filename is %s" %(fname)
                    raise ValueError("Can't write to this file.")
        
                for i in xrange(0, len(bin_middle)):
                    f.write("%.4f %.4f %.4f %d\n" %(bin_middle[i], mean[i], std[i], N_Ngamma[model_number][snapshot_idx][i]))
                f.close() 
                print "Wrote successfully to file %s" %(fname)
            ##

        for model_number in xrange(0, len(SnapList)): # Just plot some garbage to get the legend labels correct.
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
        print 'Saved file to', outputFile
    
        plt.close()


def bin_Simfast_halos(RedshiftList, SnapList, halopath, fitpath, fesc_prescription, fesc_normalization, GridSize, output_tag):
   
    for model_number in xrange(0, len(fesc_prescription)):
        for halo_z_idx in xrange(0, len(RedshiftList)):
            snapshot_idx = min(range(len(SnapList)), key=lambda i: abs(SnapList[i]-RedshiftList[halo_z_idx])) # This finds the index of the simulation redshift that most closely matches the Halo redshift.
            print "Binning Halo redshift %.4f" %(RedshiftList[halo_z_idx])
            print "For the Halo redshift %.3f the nearest simulation redshift is %.3f" %(RedshiftList[halo_z_idx], SnapList[snapshot_idx])  
            if (fesc_prescription[model_number] == 0):
                fname = "%s/fesc%d_%.3f_z%.3f.txt" %(fitpath, fesc_prescription[model_number], fesc_normalization[model_number], AllVars.SnapZ[snapshot_idx]) 
            elif (fesc_prescription[model_number] == 1 or fesc_prescription[model_number] == 2):
                fname = "%s/fesc%d_A%.3eB%.3f_z%.3f.txt" %(fitpath, fesc_prescription[model_number], fesc_normalization[model_number][0], fesc_normalization[model_number][1], AllVars.SnapZ[snapshot_idx])

            print "Reading in file %s" %(fname)
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

            names = [Halodesc_full[i][0] for i in xrange(len(Halodesc_full))]
            formats = [Halodesc_full[i][1] for i in xrange(len(Halodesc_full))]
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
            for i in xrange(0, N_Halos):
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
            print "We had %d halos (out of %d, so %.4f fraction) that had halo mass that was not covered by the Mvir-Ngamma results." %(fit_nan, N_Halos, float(fit_nan)/float(N_Halos))
            print "There were %d cells with a non-zero ionizing flux." %(len(binned_nion[binned_nion != 0]))

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
            print 'Saved file to', outputFile
            plt.close()
            
def plot_photoncount(SnapList, sum_nion, num_files, max_files, model_tags, output_tag): 
    '''
    Plots the ionizing emissivity as a function of redshift. 
    We normalize the emissivity to Mpc^-3 and this function allows the read-in of only a subset of the volume.
    Parallel compatible.
    Accepts 3D arrays of ngamma/fesc to plot Ngamma for multiple models. 

    Parameters
    ---------
    SnapList : Nested `array-like', SnapList[model_number0] = [snapshot0_model0, ..., snapshotN_model0], with length equal to the number of models.
    Snapshots for each model, defines the x-axis we plot against.
    sum_nion : Nested 1-dimensional `array-like', sum_nion[z0, z1, ..., zn], with length equal to the number of redshifts. 
    Number of escape ionizing photons (i.e., photon rate times the local escape fraction) at each redshift.
    In units of 1.0e50 s^-1.
    num_files : int
    Number of files that were read in to create the plot.
    max_files : int
    Number of files that are required to span the entire volume.
    model_tags : `array-like' of strings with length equal to the number of models.
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
    print "Plotting the ionizing emissivity." 

    sum_array = []

    for model_number in xrange(0, len(SnapList)): 
        sum_array.append([])
        for snapshot_idx in xrange(0, len(SnapList[model_number])):
       
            nion_sum_snapshot = comm.reduce(sum_nion[model_number][snapshot_idx], op = MPI.SUM, root = 0)
            if rank == 0:
                sum_array[model_number].append(nion_sum_snapshot * 1.0e50 / (pow(AllVars.BoxSize / AllVars.Hubble_h,3) * (float(num_files) / float(max_files))))
            
    if (rank == 0):
        ax1 = plt.subplot(111)

        for model_number in xrange(0, len(SnapList)):
            t = np.empty(len(SnapList[model_number]))
            for snapshot_idx in xrange(0, len(SnapList[model_number])):
                t[snapshot_idx] = (t_BigBang - cosmo.lookback_time(AllVars.SnapZ[SnapList[model_number][snapshot_idx]]).value) * 1.0e3      
    
            ax1.plot(t, np.log10(sum_array[model_number]), color = PlotScripts.colors[model_number], linestyle = PlotScripts.linestyles[model_number], label = model_tags[model_number], linewidth = PlotScripts.global_linewidth)  
            #ax1.fill_between(t, np.subtract(mean,std), np.add(mean,std), color = colors[model_number], alpha = 0.25)

        #ax1.xaxis.set_minor_locator(mtick.MultipleLocator(time_tick_interval))
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

        leg = ax1.legend(loc='lower right', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize(PlotScripts.global_legendsize)

        plot_time = 1
        bouwens_z = np.arange(6,16) # Redshift range for the observations.
        bouwens_t = (AllVars.t_BigBang - cosmo.lookback_time(bouwens_z).value) * 1.0e3 # Corresponding values for what we will plot on the x-axis.

        bouwens_1sigma_lower = [50.81, 50.73, 50.60, 50.41, 50.21, 50.00, 49.80, 49.60, 49.39, 49.18] # 68% Confidence Intervals for the ionizing emissitivity from Bouwens 2015.
        bouwens_1sigma_upper = [51.04, 50.85, 50.71, 50.62, 50.56, 50.49, 50.43, 50.36, 50.29, 50.23]

        bouwens_2sigma_lower = [50.72, 50.69, 50.52, 50.27, 50.01, 49.75, 49.51, 49.24, 48.99, 48.74] # 95% CI. 
        bouwens_2sigma_upper = [51.11, 50.90, 50.74, 50.69, 50.66, 50.64, 50.61, 50.59, 50.57, 50.55]

        if plot_time == 1:
            ax1.fill_between(bouwens_t, bouwens_1sigma_lower, bouwens_1sigma_upper, color = 'k', alpha = 0.2, label = r"$\mathrm{Bouwens \: et \: al. \: (2015)}$")
            ax1.fill_between(bouwens_t, bouwens_2sigma_lower, bouwens_2sigma_upper, color = 'k', alpha = 0.4, label = r"$\mathrm{Bouwens \: et \: al. \: (2015)}$")
        else:
            ax1.fill_between(bouwens_z, bouwens_1sigma_lower, bouwens_1sigma_upper, color = 'k', alpha = 0.2, label = r"$\mathrm{Bouwens \: et \: al. \: (2015)}$")
            ax1.fill_between(bouwens_z, bouwens_2sigma_lower, bouwens_2sigma_upper, color = 'k', alpha = 0.4, label = r"$\mathrm{Bouwens \: et \: al. \: (2015)}$")

    #   ax1.text(0.075, 0.965, '(a)', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
        ax1.text(350, 50.0, r"$68\%$", horizontalalignment='center', verticalalignment = 'center', fontsize = PlotScripts.global_labelsize) 
        ax1.text(350, 50.8, r"$95\%$", horizontalalignment='center', verticalalignment = 'center', fontsize = PlotScripts.global_labelsize)

        plt.tight_layout()
        outputFile = './' + output_tag + output_format
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
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

    ax1 = plt.subplot(121)
    ax3 = plt.subplot(122)
    #ax5 = plt.subplot(133)

    look_for_alive = 0
    #idx_array = [20004, 20005, 20016]
    #halonr_array = [7381]
    #halonr_array = [389106]
    halonr_array = [36885]
    for model_number in xrange(0, len(model_tags)):
        if(simulation_norm[model_number] == 0):
            AllVars.Set_Params_Mysim()
        elif(simulation_norm[model_number] == 3):
            AllVars.Set_Params_Tiamat_extended()
        else:
            print "Simulation norm was set to %d." %(simulation_norm[model_number])
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
            print "The galaxy that was present in the most snapshots is %d which was in %d snaps" %(np.argmax(alive), np.amax(alive)) 
            most_alive = alive.argsort()[-10:][::-1] # Finds the 3 galaxies alive for the most snapshots.  Taken from https://stackoverflow.com/questions/6910641/how-to-get-indices-of-n-maximum-values-in-a-numpy-array
            print G.HaloNr[most_alive]
            exit() 

        t = np.empty((number_snapshots[model_number])) 
 
       
        for snapshot_idx in xrange(0, number_snapshots[model_number]): 
            SFR_ensemble[model_number].append(np.mean(G.GridSFR[:,snapshot_idx]))
            ejected_ensemble[model_number].append(np.mean(G.GridOutflowRate[:, snapshot_idx]))
            infall_ensemble[model_number].append(np.mean(G.GridInfallRate[:, snapshot_idx]))
            ejectedmass_ensemble[model_number].append(np.mean(G.GridEjectedMass[:, snapshot_idx]) * 1.0e10 / AllVars.Hubble_h)
            
            t[snapshot_idx] = (t_BigBang - cosmo.lookback_time(AllVars.SnapZ[snapshot_idx]).value) * 1.0e3      
    
            
        for p in xrange(0, N_random):
            random_idx = (np.where((G.HaloNr == halonr_array[p]))[0])[0] 
            SFR_gal[model_number].append(G.GridSFR[random_idx]) # Remember the star formation rate history of the galaxy.
            ejected_gal[model_number].append(G.GridOutflowRate[random_idx])
            infall_gal[model_number].append(G.GridInfallRate[random_idx])
            ejectedmass_gal[model_number].append(G.GridEjectedMass[random_idx] * 1.0e10 / AllVars.Hubble_h)
            for snapshot_idx in xrange(0, number_snapshots[model_number]):  
                if snapshot_idx == 0:
                    pass 
                elif(G.GridHistory[random_idx, snapshot_idx] == -1):
                    print snapshot_idx
                    print len(SFR_gal[model_number])
                    print np.shape(SFR_gal[model_number])
                    print SFR_gal[model_number][p][snapshot_idx]
                    SFR_gal[model_number][p][snapshot_idx] = SFR_gal[model_number][p][snapshot_idx - 1]
        
       
        


        ax1.plot(t, SFR_ensemble[model_number], color = PlotScripts.colors[0], linestyle = PlotScripts.linestyles[model_number], label = model_tags[model_number], linewidth = PlotScripts.global_linewidth)
        ax3.plot(t, ejected_ensemble[model_number], color = PlotScripts.colors[1], linestyle = PlotScripts.linestyles[model_number], linewidth = PlotScripts.global_linewidth, alpha = 1.0)
        #ax5.plot(t, infall_ensemble[model_number], color = PlotScripts.colors[2], linestyle = PlotScripts.linestyles[model_number], linewidth = PlotScripts.global_linewidth, alpha = 1.0)
        #ax5.plot(t, ejectedmass_ensemble[model_number], color = PlotScripts.colors[2], linestyle = PlotScripts.linestyles[model_number], linewidth = PlotScripts.global_linewidth, alpha = 1.0)
    
        print ejected_ensemble[model_number]
    
        for p in xrange(0, N_random):
            ax1.plot(t, SFR_gal[model_number][p], color = PlotScripts.colors[0], linestyle = PlotScripts.linestyles[model_number], alpha = 0.5, linewidth = 1)
            ax3.plot(t, ejected_gal[model_number][p], color = PlotScripts.colors[1], linestyle = PlotScripts.linestyles[model_number], alpha = 0.5, linewidth = 1)
            #ax5.plot(t, infall_gal[model_number][p], color = PlotScripts.colors[2], linestyle = PlotScripts.linestyles[model_number], alpha = 0.5, linewidth = 1)
            #ax5.plot(t, ejectedmass_gal[model_number][p], color = PlotScripts.colors[2], linestyle = PlotScripts.linestyles[model_number], alpha = 0.5, linewidth = 1)

            #ax1.plot(t, SFR_gal[model_number][p], color = PlotScripts.colors[0], linestyle = PlotScripts.linestyles[model_number], alpha = 1.0, linewidth = 1, label = model_tags[model_number])
            #ax1.plot(t, ejected_gal[model_number][p], color = PlotScripts.colors[1], linestyle = PlotScripts.linestyles[model_number], alpha = 1.0, linewidth = 1, label = model_tags[model_number])

    #print abs(ejected_ensemble[1] - ejected_ensemble[0])/ejected_ensembles[0]

    #ax1.plot(np.nan, np.nan, color = PlotScripts.colors[0], label = 'SFR')
    #ax1.plot(np.nan, np.nan, color = PlotScripts.colors[1], label = 'Outflow')

    ax1.set_yscale('log', nonposy='clip')
    ax1.set_ylabel(r"$\mathrm{SFR} \: [\mathrm{M}_\odot \mathrm{yr}^{-1}]$")
    ax1.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = PlotScripts.global_fontsize)
    ax1.set_xlim(PlotScripts.time_xlim)
    ax1.set_ylim([1e-8, 1e3])

    
    ax3.set_yscale('log', nonposy='clip')
    ax3.set_ylabel(r"$\mathrm{Outflow \: Rate} \: [\mathrm{M}_\odot \mathrm{yr}^{-1}]$")
    ax3.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = PlotScripts.global_fontsize)
    ax3.set_xlim(PlotScripts.time_xlim)
    ax3.set_ylim([1e-8, 1e3])

    '''
    ax5.set_yscale('log', nonposy='clip')
    #ax5.set_ylabel(r"$\mathrm{Infall \: Rate} \: [\mathrm{M}_\odot \mathrm{yr}^{-1}]$")
    ax5.set_ylabel(r"$\mathrm{Ejected Mass} [\mathrm{M}_\odot]$")
    ax5.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = PlotScripts.global_fontsize)
    ax5.set_xlim(PlotScripts.time_xlim)
    #ax5.set_ylim([1e-8, 1e3])
    ax5.set_ylim([1e6, 1e10])
    '''
    ax2 = ax1.twiny()
    ax4 = ax3.twiny()
    #ax6 = ax5.twiny()

    t_plot = (t_BigBang - cosmo.lookback_time(PlotScripts.z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
    z_labels = ["$%d$" % x for x in PlotScripts.z_plot] # Properly Latex-ize the labels.

    ax2.set_xlabel(r"$z$", size = PlotScripts.global_labelsize)
    ax2.set_xlim(PlotScripts.time_xlim)
    ax2.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
    ax2.set_xticklabels(z_labels) # But label them as redshifts.

    ax4.set_xlabel(r"$z$", size = PlotScripts.global_labelsize)
    ax4.set_xlim(PlotScripts.time_xlim)
    ax4.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
    ax4.set_xticklabels(z_labels) # But label them as redshifts.

    '''
    ax6.set_xlabel(r"$z$", size = PlotScripts.global_labelsize)
    ax6.set_xlim(PlotScripts.time_xlim)
    ax6.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
    ax6.set_xticklabels(z_labels) # But label them as redshifts.
    '''

    plt.tight_layout()
    leg = ax1.legend(loc='upper left', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize)


    outputFile = './Halo%d_%s%s' %(halonr_array[0], output_tag, output_format) 
    plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    print 'Saved file to', outputFile
    plt.close()


### Here ends the plotting functions. ###
### Here begins the functions that calculate various properties for the galaxies (fesc, Magnitude etc). ###

def Calculate_HaloPartStellarMass(halo_part, stellar_mass, bound_low, bound_high):
    '''
    Calculates the stellar mass for galaxies whose host halos contain a specified number of particles.

    Parameters
    ----------
    halo_part : `array-like'
    Array containing the number of particles inside each halo.
    stella_rmass : `array-like'
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


def calculate_fesc(fesc_prescription, halo_mass, ejected_fraction, fesc_normalization):
    '''
    Calculate the escape fraction for a given prescription.

    Parameters
    ----------
    fesc_prescription : int
    Number that controls what escape fraction prescription we are using.
    0 : Constant, fesc = Constant.
    1 : Scaling with Halo Mass, fesc = A*Mh^B.
    2 : Scaling with ejected fraction, fesc = fej*A + B.
    halo_mass : `numpy.darray'
    Array that contains the halo masses for this snapshot. Units are log10(Msun).
    ejected_fraction : `numpy.darray'
    Array that contains the ejected fraction of the galaxies for this snapshot.
    fesc_normalization : float (if fesc_prescription == 0) or `numpy.darray' with length 2 (if fesc_prescription == 1 or == 2).
    If fesc_prescription == 0, gives the constant value for the escape fraction.
    If fesc_prescription == 1 or == 2, gives A and B with the form [A, B].

    Returns
    -------
    fesc : `array-like', length is same as input arrays.
    Array containing the escape fraction for each galaxy for this snapshot.

    Units
    -----
    Input Halo Mass is in units of log10(Msun).
    '''

    if(len(halo_mass) != len(ejected_fraction)):
        print "The length of the halo_mass array is %d and the length of the ejected_fraction array is %d" %(len(halo_mass), len(ejected_fraction))
        raise ValueError("These two should be equal.")

    if(fesc_prescription > 2):
        print "The prescription value you've chosen is %d" %(fesc_prescription)
        raise ValueError("Currently we only have prescriptions 0 (constant), 1 (scaling with halo mass) and 2 (scaling with ejected fraction).")

    if((fesc_prescription == 0 and isinstance(fesc_normalization, list) == True)):
        print "The escape fraction prescription is 0 but fesc_normalization was a list with value of ", fesc_normalization
        raise ValueError("For a constant escape fraction, fesc_noramlization should be a single float.")

    if((fesc_prescription == 1 or fesc_prescription == 2) and (isinstance(fesc_normalization, list) == False)):
        print "The escape fraction prescription is ", fesc_prescription, " but fesc_normalization was not a list; it instead had a value of ", fesc_normalization
        raise ValueError("For a scaling escape fraction fesc_normalization should be a list of the form [A, B]") 
 
    if (fesc_prescription == 0):
        fesc = np.full((len(halo_mass)), fesc_normalization)
    elif (fesc_prescription == 1):
        fesc = fesc_normalization[0]*pow(10,halo_mass*fesc_normalization[1])
    elif (fesc_prescription == 2):
        fesc = fesc_normalization[0]*ejected_fraction + fesc_normalization[1]

    if len(fesc) != 0:
    ## Adjust bad values, provided there isn't a riduculous amount of them. ##  
        w = np.where((halo_mass >= m_low) & (halo_mass <= m_high))[0] # Only care about the escape fraction values between these bounds. 
        nan_values = len(fesc[w][np.isnan(fesc[w])]) 
        aboveone_values = len(fesc[w][fesc[w] > 1.0])
        belowzero_values = len(fesc[w][fesc[w] < 0.0])
#   print "There was %d escape fraction values that were NaN, %d > 1.0 and %d < 0.0.  There was %d in total." %(nan_values, aboveone_values, belowzero_values, len(fesc))
        bad_ratio =  (float(nan_values) + float(aboveone_values) + float(belowzero_values)) / float(len(fesc))
        
        if (bad_ratio > 0.10):
            print "The ratio of bad escape fractions to good is %.4f" %(bad_ratio)
            print "Values above 1: ", fesc[fesc > 1.0]
            w = np.where(fesc > 1.0)[0]
            print "Halo mass is; ", halo_mass[w]
            print "Ejected fraction is: ", ejected_fraction[w]
            print
            print "Values below 0: ", fesc[fesc < 0.0]
            w = np.where(fesc < 0.0)[0]
            print "Halo mass is; ", halo_mass[w]
            print "Ejected fraction is: ", ejected_fraction[w]
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

    if len(SFR[SFR == 0.0] > 0):
        print SFR
        raise ValueError("When trying to calculate the number of photons I encountered a galaxy with SFR = 0.0.  This should've been taken care of in the np.where() condition.") 

    for i in xrange(0, len(SFR)):
        ngamma_HI_tmp = 0.0
        if (Z[i] < 0.0025):
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
    L, M : `array-like', length equal to the number of galaxies at this snapshot.
    Array containing the UV luminosities and magnitudes.

    Returns
    -------
    M_UV_obs : `array-like', length equal to the number of galaxies at this snapshot.
    Array containing the observed UV magnitudes.

    Units
    -----
    Luminosities are in units of log10(erg s^-1 A^-1).
    Magnitudes are in the AB system.
    '''

    M_UV_bins = np.arange(-24, -16, 0.1)
    A_mean = np.zeros((len(MUV_bins))) # A_mean is the average UV extinction for a given UV bin.    

    for j in xrange(0, len(M_UV_bins)):
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
    mean_pool, std_pool, N_pool : `array-like' of floats with length equal to the number of bins (e.g. the mass bins for the Stellar Mass Function).
    The current mean, standard deviation and number of data points within in each bin.  This is the array that will be updated in this function.
    mean_local, std_local, N_local : `array-like' of floats with length equal to the number of bins.
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
 
    for i in xrange(0, len(N_pool)):
        if(N_pool[i] == 0): # This case is when we have no data points in the bin. 
            mean_pool[i] = 0.0
        else:
            mean_pool[i] = N_times_mean_pool[i]/N_pool[i]
        if(N_pool[i] < 3): # In this instance we don't have enough data points to properly calculate the standard deviation.
            std_pool[i] = 0.0
        else:
            std_pool[i] = np.sqrt(N_times_var_pool[i]/ (N_pool[i] - 2)) # We have -2 because there is two instances of N_pool contains two 'N - 1' terms. 
         
    return mean_pool, std_pool

    ### Here ends the functions that deal with galaxy data manipulation. ###

#################################
### Here starts the main body of the code. ###

HaloPart_Low = 50 # Bounds for where we define the cutoff for a 'Dark Matter Halo'. Below this we can't be sure of the results.
HaloPart_High = 52

calculate_observed_LF = 0

### The arrays in this block control constants for each model ##

number_models = 2

galaxies_model1 = '/lustre/projects/p004_swin/jseiler/september/galaxies/tiamat_IRA_10step_z1.827'
#galaxies_model2 = '/lustre/projects/p004_swin/jseiler/september/galaxies/tiamat_Delayed5Myr_10step_z1.827'
galaxies_model3 = '/lustre/projects/p004_swin/jseiler/september/galaxies/tiamat_Delayed5Myr_10step_z1.827'
#galaxies_model2 = '/lustre/projects/p004_swin/jseiler/september/galaxies/tiamat_Delayed10Myr_10step_z1.827'
#galaxies_model1 = '/lustre/projects/p004_swin/jseiler/september/galaxies/tiamat_Delayed5Myr_SF0.0075_z1.827'
#galaxies_model2 = '/lustre/projects/p004_swin/jseiler/18month/galaxies/tiamat_test_z1.827'
#galaxies_model2 = '/lustre/projects/p004_swin/jseiler/september/galaxies/tiamat_Delayed5Myr_SF0.01_z1.827'
#galaxies_model1 = '/lustre/projects/p004_swin/jseiler/september/galaxies/mysim_IRA_z5.000'
#galaxies_model2 = '/lustre/projects/p004_swin/jseiler/september/galaxies/mysim_Delayed5Myr_z5.000'

#merged_galaxies_model1 = '/lustre/projects/p004_swin/jseiler/september/galaxies/mysim_IRA_MergedGalaxies'
#merged_galaxies_model2 = '/lustre/projects/p004_swin/jseiler/september/galaxies/mysim_Delayed5Myr_MergedGalaxies'

merged_galaxies_model1 = '/lustre/projects/p004_swin/jseiler/september/galaxies/tiamat_IRA_10step_MergedGalaxies'
#merged_galaxies_model2 = '/lustre/projects/p004_swin/jseiler/september/galaxies/tiamat_Delayed5Myr_10step_MergedGalaxies'
merged_galaxies_model3 = '/lustre/projects/p004_swin/jseiler/september/galaxies/tiamat_Delayed5Myr_10step_MergedGalaxies'
#merged_galaxies_model2 = '/lustre/projects/p004_swin/jseiler/september/galaxies/tiamat_Delayed10Myr_10step_MergedGalaxies'
#merged_galaxies_model1 = '/lustre/projects/p004_swin/jseiler/september/galaxies/tiamat_Delayed5Myr_SF0.0075_MergedGalaxies'
#merged_galaxies_model2 = '/lustre/projects/p004_swin/jseiler/18month/galaxies/tiamat_test_MergedGalaxies'
#merged_galaxies_model2 = '/lustre/projects/p004_swin/jseiler/september/galaxies/tiamat_Delayed5Myr_SF0.01_MergedGalaxies'

galaxies_filepath_array = [galaxies_model1, galaxies_model3]
merged_galaxies_filepath_array = [merged_galaxies_model1, merged_galaxies_model3]

number_snapshots = [164, 164] # Number of snapshots in the simulation (we don't have to do calculations for ALL snapshots).
# Tiamat extended has 164 snapshots.
FirstFile = [0, 0] # The first file number THAT WE ARE PLOTTING.
LastFile = [9, 9] # The last file number THAT WE ARE PLOTTING.
NumFile = [27, 27] # The number of files for this simulation (plotting a subset of these files is allowed). 

#model_tags = [r"$f_\mathrm{esc} = \mathrm{Constant}$", r"$f_\mathrm{esc} \: \propto \: M_\mathrm{H}^{-1}$", r"$f_\mathrm{esc} \: \propto \: M_\mathrm{H}$", r"$f_\mathrm{esc} \: \propto \: f_\mathrm{ej}$"]
#model_tags = [r"No SN", r"Delayed - 5Myr", r"Delayed - 10 Myr"]
#model_tags = [r"IRA", r"Delayed - 5Myr", r"Delayed - 10Myr"]
model_tags = [r"IRA", r"Delayed"]
#model_tags = [r"IRA - 1 Step", r"IRA - 10 Step"]
#model_tags = [r"Delayed - 1 Step", r"Delayed - 10 Step"]
#model_tags =  [r"Delayed5Myr - $\alpha = 0.0075$", r"Delayed5Myr - $\alpha = 0.01$"]

for model_number in xrange(0,number_models):
    assert(LastFile[model_number] - FirstFile[model_number] + 1 >= size)

## Constants used for each model. ##
# Need to add an entry for EACH model. #

sSFR_min = [1.0e100, 1.0e100, 1.0e100, 1.0e100]
sSFR_max = [-1.0e100, -1.0e100, 1.0e100, 1.0e100]
halo_cut = [50, 50, 50, 100] # Only calculate galaxy properties whose host halo has particle number greater than this.
source_efficiency = [1, 1, 1, 1] # Used for the halo based prescription for ionizing emissivity.

fesc_lyman_alpha = [0.3, 0.3, 0.3, 0.3] # Escape fraction of Lyman Alpha photons.
fesc_prescription = [0, 1, 2] # 0 is constant, 1 is scaling with halo mass, 2 is scaling with ejected fraction.

## Normalizations for the escape fractions. ##
# For prescription 0, requires a number that defines the constant fesc.
# For prescription 1, fesc = A*M^B. Requires an array with 2 numbers the first being A and the second B.
# For prescription 2, fesc = A*fej + B.  Requires an array with 2 numbers the first being A and the second B.
fesc_normalization = [0.10, [1000.0, -0.4], [0.2, 0.00]] 
#fesc_normalization = [0.50, [1000.0, -0.4], [1.0, 0.0]] 
##

SnapList =  [[78, 64, 51], [78, 64, 51]]
#SnapList =  [[163], [163], [163]]
# z = [6, 7, 8] are snapshots [78, 64, 51]
simulation_norm = [3, 3, 3] # 0 for MySim, 1 for Mini-Millennium, 2 for Tiamat (up to z =5), 3 for extended Tiamat (down to z = 1.6ish).

galaxy_halo_mass_lower = [95, 95, 95, 95] # These limits are for the number of particles in a halo.  
galaxy_halo_mass_upper = [105, 105, 105, 105] # We calculate the average stellar mass for galaxies whose host halos have particle count between these limits.

##############################################################################################################

files = [galaxies_filepath_array, merged_galaxies_filepath_array, number_snapshots, model_tags, SnapList]
goodness_check = [len(x) == len(y) for i,x in enumerate(files) for j,y in enumerate(files) if i != j] # This goes through the files array and checks to see if all the files have the same lengths.

if False in goodness_check: # If any of the arrays had different lengths, throw an error.
    print galaxies_filepath_array
    print merged_galaxies_filepath_array
    print number_snapshots
    print model_tags
    print SnapList
    raise ValueError("All the printed arrays (galaxies_filepath_array, merged_galaxies_filepath_array, number_snapshots, SnapList and model_tags respectively) should have the same length.")

if number_models != len(galaxies_filepath_array):
    print "The number of models was given as %d whereas we had %d filepaths given." %(number_models, len(galaxies_filepath_array))

#bin_Simfast_halos(np.arange(13.500, 5.9, -0.250), AllVars.SnapZ, "/lustre/projects/p004_swin/jseiler/Simfast21/Halos/", "/lustre/projects/p004_swin/jseiler/tiamat/halo_ngamma", fesc_prescription, fesc_normalization, 512, "Simfast_nion")
#exit()

SMF = []

## Arrays for functions of halo mass. ##
mean_ejected_halo_array = []
std_ejected_halo_array = []
mean_fesc_halo_array = []
std_fesc_halo_array = []
mean_Ngamma_halo_array = []
std_Ngamma_halo_array = []
N_array = []

## Arrays for functions of redshift. ##
sum_Ngamma_z_array = []
mean_fesc_z_array = []
std_fesc_z_array = []
N_z = []

AllVars.Set_Params_Tiamat_extended()
#for i in xrange(0, len(AllVars.SnapZ)-1):
#    print "Snapshot ", i, "and Snapshot", i + 1, "have a time difference of", (AllVars.Lookback_Time[i] - AllVars.Lookback_Time[i+1]) * 1.0e3, "Myr"

plot_singleSFR(galaxies_filepath_array, merged_galaxies_filepath_array, number_snapshots, simulation_norm, model_tags, "singleSFR_ejectedmass_IRAcomp_10step")
exit()
for model_number in xrange(0, number_models):

    if(simulation_norm[model_number] == 0):
        AllVars.Set_Params_Mysim()
    elif(simulation_norm[model_number] == 3):
        AllVars.Set_Params_Tiamat_extended()
    else:
        print "Simulation norm was set to %d." %(simulation_norm[model_number])
        raise ValueError("This option has been implemented yet.  Get your head in the game Jacob!")

    if (number_snapshots[model_number] != len(AllVars.SnapZ)): # Here we do a check to ensure that the simulation we've defined correctly matches the number of snapshots we have also defined. 
        print "The number_snapshots array is", number_snapshots
        print "The simulation_norm array is", simulation_norm
        print "The number of snapshots for model_number ", model_number, "has ", len(AllVars.SnapZ), "but you've said there is only", number_snapshots[model_number]
        raise ValueError("Check either that the number of snapshots has been defined properly and that the normalization option is correct.")

    ## The statistics arrays will be 2D arrays.  E.g., Calling SMF[model_number] will yield an array with the count of the number of galaxies within each of the stellar mass bins for each snapshot. ##
    SMF.append([])

    ## Halo arrays. ##
    mean_ejected_halo_array.append([])
    std_ejected_halo_array.append([])
    mean_fesc_halo_array.append([])
    std_fesc_halo_array.append([])
    mean_Ngamma_halo_array.append([])
    std_Ngamma_halo_array.append([])
    N_array.append([])
    ## Redshift arrays. ##
    sum_Ngamma_z_array.append([])
    mean_fesc_z_array.append([])
    std_fesc_z_array.append([])
    N_z.append([])

    for snapshot_idx in xrange(0, len(SnapList[model_number])): # These arrays are used for cumulative values across all files.
        SMF[model_number].append(np.zeros((NB_gal), dtype = np.int32)) # Stellar Mass Function for each snapshot.

        ## Function of halo mass arrays. ##
        mean_ejected_halo_array[model_number].append(np.zeros((NB), dtype = np.float32)) 
        std_ejected_halo_array[model_number].append(np.zeros((NB), dtype = np.float32)) 
        mean_fesc_halo_array[model_number].append(np.zeros((NB), dtype = np.float32)) 
        std_fesc_halo_array[model_number].append(np.zeros((NB), dtype = np.float32))
        mean_Ngamma_halo_array[model_number].append(np.zeros((NB), dtype = np.float32)) 
        std_Ngamma_halo_array[model_number].append(np.zeros((NB), dtype = np.float32))
        N_array[model_number].append(np.zeros((NB), dtype = np.float32)) # How many galaxies have been binned to calculate the statistics.
        ## Function of Redshift arrays. ##
        sum_Ngamma_z_array[model_number].append(0.0) # This is the sum of ionizing photons emitted at each redshift.
        mean_fesc_z_array[model_number].append(0.0)
        std_fesc_z_array[model_number].append(0.0)
        N_z[model_number].append(0.0)
    
    for fnr in xrange(FirstFile[model_number] + rank, LastFile[model_number]+1, size): # Divide up the input files across the processors.

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

        galaxy_halo_mass_mean = []
        galaxy_halo_mass_std = []
        ##  

        GG, Gal_Desc = ReadScripts.ReadGals_SAGE_DelayedSN(galaxies_filepath_array[model_number], fnr, number_snapshots[model_number], comm) # Read in the correct galaxy file. 
        G_Merged, Merged_Desc = ReadScripts.ReadGals_SAGE_DelayedSN(merged_galaxies_filepath_array[model_number], fnr, number_snapshots[model_number], comm) # Also need the merged galaxies.
        G = ReadScripts.Join_Arrays(GG, G_Merged, Gal_Desc) # Then join them together for all galaxies that existed at this Redshift. 

        current_fesc_lyman_alpha = fesc_lyman_alpha[model_number] # Just reduces a bit of extra clutter down the road.
        current_source_efficiency = source_efficiency[model_number]
        current_halo_cut = halo_cut[model_number]

        for snapshot_idx in xrange(0, len(SnapList[model_number])):
        
            current_snap = SnapList[model_number][snapshot_idx]

            w_gal = np.where((G.GridHistory[:, current_snap] != -1) & (G.GridStellarMass[:, current_snap] > 0.0) & (G.GridStellarMass[:, current_snap] < 1e5) & (G.GridCentralGalaxyMass[:, current_snap] >= m_low_SAGE) & (G.GridCentralGalaxyMass[:, current_snap] <=  m_high_SAGE) & (G.GridSFR[:, current_snap] > 0.0) & (G.LenHistory[:, current_snap] > current_halo_cut) & (G.GridSFR[:, current_snap] > 0.0))[0] # Only include those galaxies that existed at the current snapshot, had positive (but not infinite) stellar/Halo mass and Star formation rate.
        
#            print "There were %d galaxies for snapshot %d (Redshift %.4f) model %d." %(len(w_gal), current_snap, AllVars.SnapZ[current_snap], model_number)

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
        
            fesc_local = calculate_fesc(fesc_prescription[model_number], mass_central, ejected_fraction, fesc_normalization[model_number]) 
            
            galaxy_halo_mass_mean, galaxy_halo_mass_std = Calculate_HaloPartStellarMass(halo_part_count, mass_gal, galaxy_halo_mass_lower[model_number], galaxy_halo_mass_upper[model_number])

            print galaxy_halo_mass_mean

            
            photons_HI_gal = calculate_photons(SFR_gal, metallicity_gal)    
            photons_HI_gal_nonlog = [10**x for x in photons_HI_gal]
            ionizing_photons = np.multiply(photons_HI_gal_nonlog, fesc_local)

            ## We have now calculated all the base properties for galaxies within this snapshot.  Calculate the relevant statistics and put them into their arrays. ##

            (counts_local, bin_edges, bin_middle) = AllVars.Calculate_Histogram(mass_gal, bin_width, 0, m_gal_low, m_gal_high) # Bin the Stellar Mass 
            SMF[model_number][snapshot_idx] += counts_local 

            (mean_ejected_halo_local, std_ejected_halo_local, N_local, bin_middle) = Calculate_2D_Mean(mass_central, ejected_fraction, bin_width, m_low, m_high) # This bins the halo_mass (x-axis) and then calculates the mean ejected fraction (y-axis) within each of these bins.
            (mean_ejected_halo_array[model_number][snapshot_idx], std_ejected_halo_array[model_number][snapshot_idx]) = update_cumulative_stats(mean_ejected_halo_array[model_number][snapshot_idx], std_ejected_halo_array[model_number][snapshot_idx], N_array[model_number][snapshot_idx], mean_ejected_halo_local, std_ejected_halo_local, N_local) # Update the mean ejected fraction for this Snapshot.
            
            (mean_fesc_halo_local, std_fesc_halo_local, N_local, bin_middle) = Calculate_2D_Mean(mass_central, fesc_local, bin_width, m_low, m_high) # Do the same for the escape fraction as a function of halo mass.
            (mean_fesc_halo_array[model_number][snapshot_idx], std_fesc_halo_array[model_number][snapshot_idx]) = update_cumulative_stats(mean_fesc_halo_array[model_number][snapshot_idx], std_fesc_halo_array[model_number][snapshot_idx], N_array[model_number][snapshot_idx], mean_fesc_halo_local, std_fesc_halo_local, N_local) 

            mean_fesc_z_array[model_number][snapshot_idx] = np.mean(fesc_local)
            std_fesc_z_array[model_number][snapshot_idx] = np.std(fesc_local)

            (mean_Ngamma_halo_local, std_Ngamma_halo_local, N_local, bin_middle) = Calculate_2D_Mean(mass_central, ionizing_photons, bin_width, m_low, m_high) # Do the same for the escape fraction as a function of halo mass.

            mean_Ngamma_halo_local = np.divide(mean_Ngamma_halo_local, 1.0e50)
            std_Ngamma_halo_local = np.divide(std_Ngamma_halo_local, 1.0e50)

            (mean_Ngamma_halo_array[model_number][snapshot_idx], std_Ngamma_halo_array[model_number][snapshot_idx]) = update_cumulative_stats(mean_Ngamma_halo_array[model_number][snapshot_idx], std_Ngamma_halo_array[model_number][snapshot_idx], N_array[model_number][snapshot_idx], mean_Ngamma_halo_local, std_Ngamma_halo_local, N_local) 

            sum_Ngamma_z_array[model_number][snapshot_idx] += np.sum(np.divide(ionizing_photons, 1.0e50))   

            N_z[model_number][snapshot_idx] += len(w_gal)
            N_array[model_number][snapshot_idx] += N_local 



StellarMassFunction(SnapList, SMF, simulation_norm, FirstFile, LastFile, NumFile, model_tags, 0, "tiamat_IRAcomp_10step_lowz_SMF") ## PARALLEL COMPATIBLE
#plot_ejectedfraction(SnapList, mean_ejected_halo_array, std_ejected_halo_array, N_array, model_tags, "18month_Ejected_test") ## PARALELL COMPATIBLE # Ejected fraction as a function of Halo Mass 
#plot_fesc(SnapList, mean_fesc_z_array, std_fesc_z_array, N_z, model_tags, "18month_fesc_talk") ## PARALELL COMPATIBLE 
#plot_photoncount(SnapList, sum_Ngamma_z_array, 50, 125, model_tags, "18month_nion_talk") ## PARALELL COMPATIBLE
#plot_mvir_Ngamma(SnapList, mean_Ngamma_halo_array, std_Ngamma_halo_array, N_array, model_tags, "Mvir_Ngamma_test", fesc_prescription, fesc_normalization, "/lustre/projects/p004_swin/jseiler/tiamat/halo_ngamma/") ## PARALELL COMPATIBLE 

