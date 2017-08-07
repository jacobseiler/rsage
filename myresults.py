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

matplotlib.rcdefaults()
plt.rc('axes', color_cycle=['k', 'b', 'r', 'g', 'm', 'c','y', '0.5'], labelsize='x-large')
plt.rc('lines', linewidth='2.0')
# plt.rc('font', variant='monospace')
plt.rc('legend', numpoints=1, fontsize='x-large')
plt.rc('text', usetex=True)

label_size = 20
extra_size = 2
talk_fontsize = 20
talk_legendsize = 18
plt.rc('xtick', labelsize=talk_fontsize)
plt.rc('ytick', labelsize=talk_fontsize)
plt.rc('text', usetex=True)
tick_interval = 0.25
np.set_printoptions(formatter={'float': lambda x: "{0:0.10e}".format(x)})

colors = ['r', 'b', 'g', 'm', 'c', 'k']
markers = ['x', 'o', '^', 's', 'D']
linestyles = ['-', '--', '-.', ':']

AllVars.Set_Constants()
AllVars.Set_Params_Mysim()

cosmo = cosmology.FlatLambdaCDM(H0 = AllVars.Hubble_h*100, Om0 = AllVars.Omega_m) 
t_BigBang = cosmo.lookback_time(100000).value # Lookback time to the Big Bang in Gyr.
z_plot = np.arange(6, 14)  #Range of redshift we wish to plot.
time_xlim = [315, 930]
time_tick_interval = 25

output_format = ".png"

kink_low = 10.3
kink_high = 10.30000001
m_low = 8.5
m_high = 100
	
def calculate_beta(MUV, z):
	
## Calculation of the dust attenuation parameter Beta ##
## Fit values are from Bouwens (2015) ApJ 793, 115 ##
## For z = 5 and 6, Bouwens uses a piece-wise linear relationship and a linear relationship for higher redshift. ##

## INPUT ##
# MUV: A value of the absolute magnitude in the UV (generally M1600) in the AB magnitude system.
# z: Redshift value.

## OUTPUT ##
# beta: Value of the UV continuum paramaeter beta. 

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
		
def multiply(n):
	total = 1
	for i in xrange(0, len(n)):
		total *= n[i]
		print "%.4e" %(total)
	return total


def Calculate_2D_Histogram(Data_x, Data_y, Bin_Width, min_hist_x = None, max_hist_x = None):		

    if (min_hist_x == None): 
	mi = np.floor(min(Data)) - 10*Bin_Width
	ma = np.floor(max(Data)) + 10*Bin_Width
    else:
	mi = min_hist_x 
	ma = max_hist_x

    print "The minimum of the data being binned is %.4e" %(mi)
    print "The maximum of the data being binned is %.4e" %(ma) 
	    
    NB = (ma - mi) / Bin_Width 

    bins = np.arange(mi, ma + Bin_Width, Bin_Width)
    bins_mid = bins + Bin_Width/2.0
    bins_mid = bins_mid[:-1] # len(bins_mid) should be 1 less than len(bins) as the last bin doesn't have a midpoint.	

    bins_x = np.digitize(Data_x, bins)
    N_data_y = np.zeros((len(bins)))  
 
    data_y_binned = []
    for i in xrange(0, len(bins)): 
	data_y_binned.append([])
    for i in xrange(0, len(Data_x)):
	idx = bins_x[i]
	if idx == len(data_y_binned):
		idx -= 1

    	data_y_binned[idx].append(Data_y[i])
	N_data_y[idx] += 1 


    mean_data_y = np.zeros((len(bins))) 
    std_data_y = np.zeros((len(bins))) 


    for i in xrange(0, len(bins)):
	#print data_y_binned[i]
	if data_y_binned != None: 
    		mean_data_y[i] = np.mean(data_y_binned[i])
    		std_data_y[i] = np.std(data_y_binned[i])
	else:
		mean_data_y[i] = np.nan
		std_data_y[i] = np.nan

    return mean_data_y[:-1], std_data_y[:-1], N_data_y[:-1], bins_mid


def Calculate_Histogram(Data, Bin_Width, Weights, min_hist=None, max_hist=None):

# This calculates the counts and Bin_Edges for a given set of data.

## Input ##
# Data is an array containing the data to be binned.
# Bin_Width is the width of the bins.
# Weights is either 0 or 1.  0: Denotes that the histogram should be a frequency (count) histogram. 1: Denotes that the histogram should be a probability histogram.
# min_hist: Minimum value that will be binned.  OPTIONAL: If not specified will be given by the minimum of the data - 10 times the binwidth.
# max_hist: Minimum value that will be binned.  OPTIONAL: If not specified will be given by the maximum of the data - 10 times the binwidth.


## Output ##
# Counts: The count (either frequency or probability) in each bin.
# Bin_Edges: The location of the edges of the bins.
# Bin_Middle: The middle of the bins.

    if (min_hist == None): 
	mi = np.floor(min(Data)) - 10*Bin_Width
	ma = np.floor(max(Data)) + 10*Bin_Width
    else:
	mi = min_hist 
	ma = max_hist 

    print "The minimum of the data being binned is %.4e" %(mi)
    print "The maximum of the data being binned is %.4e" %(ma) 
	    
    NB = (ma - mi) / Bin_Width 

#    print "The total of the data being binned is %.4e" %(sum(Data[Data >= mi and Data <= ma]))

    if (Weights == 0):
        (counts, Bin_Edges) = np.histogram(Data, range=(mi, ma), bins=NB)
    else:
        weights = np.ones_like(Data)/len(Data)
        (counts, Bin_Edges) = np.histogram(Data, range=(mi, ma), bins=NB, weights = weights)

    Bin_Middle = Bin_Edges[:-1] + 0.5 * Bin_Width

    return (counts, Bin_Edges, Bin_Middle)

def Sum_Log(Array):

    Sum = 0
    for i in xrange(0, len(Array)):
    	Sum += 10**Array[i]

    return Sum


def Std_Log(Array, Mean):

    Sum = 0
    for i in xrange(0, len(Array)):
	Sum += (10**Array[i] - Mean)**2

    Sum *= 1.0/len(Array)

    Std = np.sqrt(Sum)
    return Std

def StellarMassFunction(SnapList, mass, simulation_norm, halo_part_stellar_mass, model_tags, output_tag):

    title = []
    normalization_array = []
    redshift_labels = []

    counts_array = []
    bin_middle_array = []

    for model_number in xrange(0, len(SnapList)):
	counts_array.append([])
	bin_middle_array.append([])
	redshift_labels.append([])

    ## Plot Parameters ##
    binwidth = 0.1
    Observations = 1

    Frequency = 0 # 0 for a frequency (count) histogram, 1 for a probability histogram.
    errorwidth = 2
    delta = 0.05
    caps = 5
    ##

    ## Normalization for each model. ##
    for model_number in xrange(0, len(SnapList)): 

	if (simulation_norm[model_number] == 0):
		AllVars.Set_Params_Mysim()
	elif (simulation_norm[model_number] == 1):
		AllVars.Set_Params_MiniMill()
	elif (simulation_norm[model_number] == 3):
		AllVars.Set_Params_Tiamat_extended()


        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
        normalization_array.append(norm)

	## Redshift labels for each redshift within each model. ##
	for snapshot_idx in xrange(0, len(SnapList[model_number])):
		print "Doing Snapshot %d" %(SnapList[model_number][snapshot_idx])
		tmp = 'z = %.2f' %(AllVars.SnapZ[SnapList[model_number][snapshot_idx]])
		redshift_labels[model_number].append(tmp)

		#minimum_mass = np.floor(min(mass[model_number][snapshot_idx])) - 10*binwidth
		#maximum_mass = np.floor(max(mass[model_number][snapshot_idx])) + 10*binwidth

		minimum_mass = -2 
		maximum_mass = 10 

		binning_minimum = comm.allreduce(minimum_mass, op = MPI.MIN)
		binning_maximum = comm.allreduce(maximum_mass, op = MPI.MAX)

		print "I am rank %d and my binning_minimum is %.3f and my binning_maximum is %.3f" %(rank, binning_minimum, binning_maximum)

		(counts_local, bin_edges, bin_middle) = Calculate_Histogram(mass[model_number][snapshot_idx], binwidth, Frequency, binning_minimum, binning_maximum)
		#counts_array[model_number].append(counts)
		#bin_middle_array[model_number].append(bin_middle) 

		if rank == 0:
			counts_total = np.zeros_like(counts_local)
		else:
			counts_total = None
		comm.Reduce([counts_local, MPI.DOUBLE], [counts_total, MPI.DOUBLE], op = MPI.SUM, root = 0)

		if rank == 0:
			counts_array[model_number].append(counts_total)
			bin_middle_array[model_number].append(bin_middle)

#			for p in xrange(1, size):
#				counts_array[model_number][snapshot_idx] += comm.recv(source = p, tag = 2)
#			comm.send(counts_array[model_number][snapshot_idx], dest = 0, tag = 2)	

	

### Plotting ###

    if rank == 0:
    	f = plt.figure()  
	ax = plt.subplot(111)  

	for model_number in xrange(0, len(SnapList)):
		for snapshot_idx in xrange(0, len(SnapList[model_number])):
			if model_number == 0:
				title = redshift_labels[model_number][snapshot_idx]
			else:
				title = ''
			print counts_array[model_number][snapshot_idx]	
			plt.plot(bin_middle_array[model_number][snapshot_idx], counts_array[model_number][snapshot_idx] / normalization_array[model_number], color = colors[snapshot_idx], linestyle = linestyles[model_number], rasterized = True, label = title)
	
##

	
    	for model_number in xrange(0, len(SnapList)):
		plt.plot(1e100, 1e100, color = 'k', ls = linestyles[model_number], label = model_tags[model_number], rasterized=True)
	
### Draws a vertical line to denote lower bounds for what is an 'acceptable' Stellar Mass ### 

#    plt.axvline(x = np.log10(HaloPartStellarMass), ymin = 0, ymax = 10, linestyle = '-.', rasterized=True)

## 

    	plt.yscale('log', nonposy='clip')

    	plt.axis([6, 11.5, 1e-6, 1e-1])

    	ax.set_xlabel(r'$\log_{10}\ m_{\mathrm{*}} \:[M_{\odot}]$', fontsize = talk_fontsize)
    	ax.set_ylabel(r'$\Phi\ [\mathrm{Mpc}^{-3}\: \mathrm{dex}^{-1}]$', fontsize = talk_fontsize)
    	ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

### If we want to put observations on the figure ###

    if (Observations == 1 and rank == 0):
    	Gonzalez_z6 = np.array([[7.77 + delta, -2.0956, -1.8596, -2.3539],
                                [8.27 + delta, -2.1742, -1.9494, -2.4101],
                                [8.77 + delta, -2.5674, -2.3876, -2.7921],
                                [9.27 + delta, -2.8483, -2.6573, -3.0843],
                                [9.77 + delta, -3.5787, -3.3764, -3.8258],
                                [10.27 + delta, -4.3202, -4.0281, -4.5674]], dtype = np.float32)

                                #plt.errorbar(Gonzalez_z6[:,0], 10**Gonzalez_z6[:,1], yerr= (10**Gonzalez_z6[:,3], 10**Gonzalez_z6[:,2]) , alpha=0.75, lw=1.0, marker='o', ls='none', label = 'Gonzalez 2011, z = 6', color = 'cyan')


        Gonzalez_z7 = np.array([[7.75, -2.1828, -1.7463, -2.5858],
                                [8.26, -2.25, -1.8694, -2.2631],
                                [8.77, -2.7425, -2.3731, -3.1231],
                                [9.27, -3.0672, -2.6753, -3.4142],
                                [9.76, -3.8731, -3.4831, -4.2537]], dtype = np.float32)

#plt.errorbar(Gonzalez_z7[:,0], 10**Gonzalez_z7[:,1], yerr= (10**Gonzalez_z7[:,3], 10**Gonzalez_z7[:,2]) , alpha=0.75, lw=1.0, marker='o', ls='none'    , label = 'Gonzalez 2011, z = 7', color = 'magenta')

        Song_z6 = np.array([[7.25 - delta, -1.47, -1.47 + 0.35, -1.47 - 0.23],
                            [7.75 - delta, -1.81, -1.81 + 0.23, -1.81 - 0.28],
                            [8.25 - delta, -2.26, -2.26 + 0.21, -2.26 - 0.16],
                            [8.75 - delta, -2.65, -2.65 + 0.15, -2.65 - 0.15],
                            [9.25 - delta, -3.14, -3.14 + 0.12, -3.14 - 0.11],
                            [9.75 - delta, -3.69, -3.69 + 0.12, -3.69 - 0.13],
                            [10.25 - delta, -4.27, -4.27 + 0.38, -4.27 - 0.86]], dtype = np.float32)

        plt.errorbar(Song_z6[:,0], 10**Song_z6[:,1], yerr= (10**Song_z6[:,1] - 10**Song_z6[:,3], 10**Song_z6[:,2] - 10**Song_z6[:,1]), xerr = 0.25, capsize = caps, elinewidth = errorwidth, alpha = 1.0, lw=2.0, marker='o', ls='none', label = 'Song 2015, z = 6', color = '#bd0026', rasterized=True)

        Song_z7 = np.array([[7.25, -1.63, -1.63 + 0.54, -1.63 - 0.54],
                            [7.75, -2.07, -2.07 + 0.45, -2.07 - 0.41],
                            [8.25, -2.49, -2.49 + 0.38, -2.49 - 0.32],
                            [8.75, -2.96, -2.96 + 0.32, -2.96 - 0.30],
                            [9.25, -3.47, -3.47 + 0.32, -3.47 - 0.35],
                            [9.75, -4.11, -4.11 + 0.41, -4.11 - 0.57],
                            [10.25, -4.61, -4.61 + 0.72, -4.61 - 0.82],
                            [10.75, -5.24, -5.24 + 0.90, -5.25 - 0.57]], dtype = np.float32)

        plt.errorbar(Song_z7[:,0], 10**Song_z7[:,1], yerr= (10**Song_z7[:,1] - 10**Song_z7[:,3], 10**Song_z7[:,2] - 10**Song_z7[:,1]), xerr = 0.25, capsize = caps, alpha=0.75, elinewidth = errorwidth, lw=1.0, marker='o', ls='none', label = 'Song 2015, z = 7', color = 'blue', rasterized=True)

        Song_z8 = np.array([[7.25, -1.73, -1.73 + 1.01, -1.73 - 0.84],
                            [7.75, -2.28, -2.28 + 0.84, -2.28 - 0.64],
                            [8.25, -2.88, -2.88 + 0.75, -2.88 - 0.57],
                            [8.75, -3.45, -3.45 + 0.57, -3.45 - 0.60],
                            [9.25, -4.21, -4.21 + 0.63, -4.21 - 0.78],
                            [9.75, -5.31, -5.31 + 1.01, -5.31 - 1.64]], dtype = np.float32)


        plt.errorbar(Song_z8[:,0], 10**Song_z8[:,1], yerr= (10**Song_z8[:,1] - 10**Song_z8[:,3], 10**Song_z8[:,2] - 10**Song_z8[:,1]), xerr = 0.25, capsize = caps, alpha=0.75, elinewidth = errorwidth, lw=1.0, marker='o', ls='none', label = 'Song 2015, z = 8', color = 'green', rasterized=True)

    if (Observations == 2 and rank == 0):
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



    if (rank == 0):
    	leg = plt.legend(loc='lower left', numpoints=1, labelspacing=0.1)
    	leg.draw_frame(False)  # Don't want a box frame
    	for t in leg.get_texts():  # Reduce the size of the text
        	t.set_fontsize(talk_legendsize)

    	#plt.tight_layout()


	outputFile = './%s%s' %(output_tag, output_format)
    	#f.savefig("foo.pdf", bbox_inches='tight')
    	plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    	print 'Saved file to', outputFile
    	plt.close()
 

#####

def calculate_pooled_stats(pooled_mean, pooled_std, mean_local, std_local, N_local):

    mean_times_N = mean_local*N_local
    var_times_N = pow(std_local,2)*N_local
	
    if (isinstance(mean_local, list) == True):   
	if rank == 0:
		pooled_mean_numerator = np.zeros_like(mean_times_N)
		pooled_std_numerator = np.zeros_like(var_times_N)
		denominator = np.zeros_like(N_local)
	else:
		pooled_mean_numerator = None
		pooled_std_numerator = None
		denominator = None

	comm.Reduce([mean_times_N, MPI.DOUBLE], [pooled_mean_numerator, MPI.DOUBLE], op = MPI.SUM , root = 0)
	comm.Reduce([var_times_N, MPI.DOUBLE], [pooled_std_numerator, MPI.DOUBLE], op = MPI.SUM, root = 0)
	comm.Reduce([N_local, MPI.DOUBLE], [denominator, MPI.DOUBLE], op = MPI.SUM, root = 0)
    else:
    	pooled_mean_numerator = comm.reduce(mean_times_N, op = MPI.SUM, root = 0)
    	pooled_std_numerator = comm.reduce(var_times_N, op = MPI.SUM, root = 0)
    	denominator = comm.reduce(N_local, op = MPI.SUM, root = 0)
		
    if rank == 0:
	if (isinstance(mean_local, list) == True):
		pooled_std_sum = np.zeros_like(mean_local)
	else:
    		pooled_std_sum = 0.0
    	for p in xrange(1, size):
		if (isinstance(mean_local, list) == True):	
			mean_other_process = np.zeros_like(mean_local)
			N_other_process = np.zeros_like(mean_local)
			
			comm.Recv(mean_other_process, source = p, tag = 3)
			comm.Recv(N_other_process, source = p, tag = 4)
	
			pooled_std_sum += N_other_process * N_local * pow(np.subtract(mean_other_process, mean_local),2)
			
		else:
			
			mean_other_process = comm.recv(source = p, tag = 3)
			N_other_process = comm.recv(source = p, tag = 4)

			pooled_std_sum += N_other_process*N_local * pow(mean_other_process - mean_local, 2)

    else:		
    		if (isinstance(mean_local, list) == True):
			comm.Send(mean_local, dest = 0, tag = 3)
			comm.Send(N_local, dest = 0, tag = 4)
		else:
			comm.send(mean_local, dest = 0, tag = 3)
			comm.send(N_local, dest = 0, tag = 4)
	
    if (rank == 0):	
    	if (isinstance(mean_local, list) == True):
		pooled_mean.append(np.divide(pooled_mean_numerator, denominator))
		pooled_std.append(np.add(np.sqrt(np.divide(pooled_std_numerator,denominator), np.divide(pooled_std_sum,pow(denominator,2)))))
	else:
		pooled_mean.append(pooled_mean_numerator / denominator)
		pooled_std.append(np.sqrt(pooled_std_numerator / denominator + pooled_std_sum/pow(denominator,2)))
	return pooled_mean, pooled_std 
    else:
	return pooled_mean, pooled_std 
##

def plot_fesc(SnapList, fesc, galaxy_mass, halo_mass, model_tags, output_tag):

    mean_fesc = []
    std_fesc = []

    mean_times_N = []
    std_times_N = []

    pooled_mean_fesc = []
    pooled_std_fesc = []

    print "Plotting fesc."
    for model_number in xrange(0, len(SnapList)): 
	mean_fesc.append([])
	std_fesc.append([])

	mean_times_N.append([])	
	std_times_N.append([])

	pooled_mean_fesc.append([])
	pooled_std_fesc.append([])

    	for snapshot_idx in xrange(0, len(SnapList[model_number])):
		w = np.where(((halo_mass[model_number][snapshot_idx] < kink_low) & (halo_mass[model_number][snapshot_idx] > m_low))| ((halo_mass[model_number][snapshot_idx] > kink_high) & (halo_mass[model_number][snapshot_idx] < m_high)))[0]
		
	    	mean_fesc[model_number].append(np.mean(fesc[model_number][snapshot_idx][w])) # Mean of this slice of the fesc.
		std_fesc[model_number].append(np.std(fesc[model_number][snapshot_idx][w])) # Standard deviation of this slice of the fesc.

		## Calculation of pooled mean and standard deviation (unecessary for this case because they should be the same as the base mean/std.##
		## Taken from https://en.wikipedia.org/wiki/Pooled_variance ## 

		N = len(fesc[model_number][snapshot_idx][w])
		comm.Barrier()
		pooled_mean_fesc[model_number], pooled_std_fesc[model_number] = calculate_pooled_stats(pooled_mean_fesc[model_number], pooled_std_fesc[model_number], mean_fesc[model_number][snapshot_idx], std_fesc[model_number][snapshot_idx], N)

    if (rank == 0):
    	ax1 = plt.subplot(111)

	for model_number in xrange(0, len(SnapList)):
    		t = np.empty(len(SnapList[model_number]))
		for snapshot_idx in xrange(0, len(SnapList[model_number])):
       			t[snapshot_idx] = (t_BigBang - cosmo.lookback_time(AllVars.SnapZ[SnapList[model_number][snapshot_idx]]).value) * 1.0e3   

		
		mean = pooled_mean_fesc[model_number]
		std = pooled_std_fesc[model_number]   

    		ax1.plot(t, mean, color = colors[model_number], ls = linestyles[model_number], label = model_tags[model_number], lw = 4)  
    		ax1.fill_between(t, np.subtract(mean,std), np.add(mean,std), color = colors[model_number], alpha = 0.25)

    	ax1.xaxis.set_minor_locator(mtick.MultipleLocator(time_tick_interval))
    	ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.025))
    	ax1.set_xlim(time_xlim)
    	#ax1.set_ylim([0.0, 0.3])

    	ax2 = ax1.twiny()

    	t_plot = (t_BigBang - cosmo.lookback_time(z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
    	z_labels = ["$%d$" % x for x in z_plot] # Properly Latex-ize the labels.

    	ax2.set_xlabel(r"$z$", size = label_size)
    	ax2.set_xlim(time_xlim)
    	ax2.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
    	ax2.set_xticklabels(z_labels) # But label them as redshifts.

    	ax1.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = label_size)
    	ax1.set_ylabel(r'$f_\mathrm{esc}$', fontsize = talk_fontsize)

	leg = ax1.legend(loc=1, numpoints=1, labelspacing=0.1)
    	leg.draw_frame(False)  # Don't want a box frame
    	for t in leg.get_texts():  # Reduce the size of the text
        	t.set_fontsize(talk_legendsize)

    	plt.tight_layout()
    	outputFile = './' + output_tag + output_format
    	plt.savefig(outputFile)  # Save the figure
    	print 'Saved file to', outputFile
    	plt.close()

##

def plot_ejectedfraction(SnapList, mass_central, ejected_fraction, model_tags, output_tag): 

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
    print "Plotting the Ejected Fraction"
    
    binwidth = 0.1
    Frequency = 1  
 
    for model_number in xrange(0, len(SnapList)): 
	for snapshot_idx in xrange(0, len(SnapList[model_number])):
		print "Doing Snapshot %d" %(SnapList[model_number][snapshot_idx])
		tmp = 'z = %.2f' %(AllVars.SnapZ[SnapList[model_number][snapshot_idx]])
		redshift_labels[model_number].append(tmp)

		minimum_mass = np.floor(min(mass_central[model_number][snapshot_idx])) - 10*binwidth
		maximum_mass = np.floor(max(mass_central[model_number][snapshot_idx])) + 10*binwidth

		binning_minimum = comm.allreduce(minimum_mass, op = MPI.MIN)
		binning_maximum = comm.allreduce(maximum_mass, op = MPI.MAX)

		(mean_ejected_fraction, std_ejected_fraction, N, bin_middle) = Calculate_2D_Histogram(mass_central[model_number][snapshot_idx], ejected_fraction[model_number][snapshot_idx], binwidth, binning_minimum, binning_maximum)


		mean_ejected_array[model_number], std_ejected_array[model_number] = calculate_pooled_stats(mean_ejected_array[model_number], std_ejected_array[model_number], mean_ejected_fraction, std_ejected_fraction, N)
		mean_halomass_array[model_number], std_halomass_array[model_number] = calculate_pooled_stats(mean_halomass_array[model_number], std_halomass_array[model_number], np.mean(mass_central[model_number][snapshot_idx]), np.std(mass_central[model_number][snapshot_idx]), len(mass_central[model_number][snapshot_idx]))

		bin_middle_array[model_number].append(bin_middle)
		

    if rank == 0:
	f = plt.figure()  
	ax1 = plt.subplot(111)  

	for model_number in xrange(0, len(SnapList)):
		for snapshot_idx in xrange(0, len(SnapList[model_number])):
			if model_number == 0:
				title = redshift_labels[model_number][snapshot_idx]
			else:
				title = ''
				
			mean = mean_ejected_array[model_number][snapshot_idx]
			std = std_ejected_array[model_number][snapshot_idx]
			bin_middle = bin_middle_array[model_number][snapshot_idx]

			ax1.plot(bin_middle, mean, color = colors[snapshot_idx], linestyle = linestyles[model_number], rasterized = True, label = title)
			ax1.scatter(mean_halomass_array[model_number][snapshot_idx], np.mean(~np.isnan(mean)), color = colors[snapshot_idx], marker = 'o', rasterized = True, s = 40, lw = 3)	
			#if (len(SnapList) == 1):
    			#	ax1.fill_between(bin_middle, np.subtract(mean,std), np.add(mean,std), color = colors[snapshot_idx], alpha = 0.25)

	ax1.set_xlabel(r'$\log_{10}\ M_{\mathrm{vir}}\ [M_{\odot}]$', size = talk_fontsize) 
    	ax1.set_ylabel(r'$\mathrm{Ejected \: Fraction}$', size = talk_fontsize)
    	ax1.set_xlim([8.5, 12])
    	ax1.set_ylim([-0.05, 1.0])   

    	ax1.xaxis.set_minor_locator(mtick.MultipleLocator(0.1))
    	ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.025))
    			#ax1.set_xticklabels([])
 


		
    	for model_number in xrange(0, len(SnapList)):
		ax1.plot(1e100, 1e100, color = 'k', ls = linestyles[model_number], label = model_tags[model_number], rasterized=True)
	
	
    	leg = ax1.legend(loc=1, numpoints=1, labelspacing=0.1)
    	leg.draw_frame(False)  # Don't want a box frame
    	for t in leg.get_texts():  # Reduce the size of the text
        	t.set_fontsize('medium')

	outputFile = './' + output_tag + output_format
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
		(mean_fesc, std_fesc, N, bin_middle) = Calculate_2D_Histogram(mass_central[model_number][snapshot_idx], fesc[model_number][snapshot_idx], binwidth, binning_minimum, binning_maximum)

		mean_fesc_array[model_number], std_fesc_array[model_number] = calculate_pooled_stats(mean_fesc_array[model_number], std_fesc_array[model_number], mean_fesc, std_fesc, N)
		mean_halomass_array[model_number], std_halomass_array[model_number] = calculate_pooled_stats(mean_halomass_array[model_number], std_halomass_array[model_number], np.mean(halomass_nonlog), np.std(halomass_nonlog), len(mass_central[model_number][snapshot_idx]))

		## If want to do mean/etc of halo mass need to update script. ##
		bin_middle_array[model_number].append(bin_middle)


		'''
		if rank == 0:
			mean_ngammafesc_array[model_number][snapshot_idx] = np.nan_to_num(mean_ngammafesc_array[model_number][snapshot_idx])
			std_ngammafesc_array[model_number][snapshot_idx] = np.nan_to_num(std_ngammafesc_array[model_number][snapshot_idx])
			outfile = '/lustre/projects/p004_swin/jseiler/tiamat/mean_mvir_ngammafesc_z%.3f.dat' %(AllVars.SnapZ[SnapList[model_number][snapshot_idx]])
			with open(outfile, 'w') as f:
				np.savetxt(outfile, mean_ngammafesc_array[model_number][snapshot_idx])
			print "Saved %s" %(outfile)


			outfile = '/lustre/projects/p004_swin/jseiler/tiamat/std_mvir_ngammafesc_z%.3f.dat' %(AllVars.SnapZ[SnapList[model_number][snapshot_idx]])
			with open(outfile, 'w') as f:
				np.savetxt(outfile, std_ngammafesc_array[model_number][snapshot_idx])
			print "Saved %s" %(outfile)
		'''
	#std_ngammafesc_array[model_number] = 0.434 * np.divide(std_ngammafesc_array[model_number], mean_ngammafesc_array[model_number])
	
	#mean_ngammafesc_array[model_number] = np.log10(~np.isnan(mean_ngammafesc_array[model_number]))
	
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

	ax1.set_xlabel(r'$\log_{10}\ M_{\mathrm{vir}}\ [M_{\odot}]$', size = talk_fontsize) 
    	ax1.set_ylabel(r'$f_\mathrm{esc}$', size = talk_fontsize)
    	#ax1.set_xlim([8.5, 12])
    	#ax1.set_ylim([0.0, 1.0])   

    	ax1.xaxis.set_minor_locator(mtick.MultipleLocator(0.1))
#    	ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
 	
#    	ax1.set_yscale('log', nonposy='clip')
#    	for model_number in xrange(0, len(SnapList)):
#		ax1.plot(1e100, 1e100, color = 'k', ls = linestyles[model_number], label = model_tags[model_number], rasterized=True)
	
	
    	leg = ax1.legend(loc='upper left', numpoints=1, labelspacing=0.1)
    	leg.draw_frame(False)  # Don't want a box frame
    	for t in leg.get_texts():  # Reduce the size of the text
        	t.set_fontsize('medium')

	outputFile = './' + output_tag + output_format
    	plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    	print 'Saved file to', outputFile
    	plt.close()



##

def plot_mvir_Ngamma(SnapList, mass_central, Ngamma, fesc, model_tags, output_tag): 

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
    print "Plotting Ngamma*fesc against Mvir" 
    
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


		Ngamma_nonlog = [10**x for x in Ngamma[model_number][snapshot_idx]]
		halomass_nonlog = [10**x for x in mass_central[model_number][snapshot_idx]]
		(mean_ngammafesc, std_ngammafesc, N, bin_middle) = Calculate_2D_Histogram(mass_central[model_number][snapshot_idx], np.multiply(Ngamma_nonlog, fesc[model_number][snapshot_idx]), binwidth, binning_minimum, binning_maximum)

		print "mean", mean_ngammafesc
		print "std", std_ngammafesc
		print "mean/std", std_ngammafesc/mean_ngammafesc
		print 0.434*std_ngammafesc/mean_ngammafesc


		mean_ngammafesc_array[model_number], std_ngammafesc_array[model_number] = calculate_pooled_stats(mean_ngammafesc_array[model_number], std_ngammafesc_array[model_number], mean_ngammafesc, std_ngammafesc, N)
		mean_halomass_array[model_number], std_halomass_array[model_number] = calculate_pooled_stats(mean_halomass_array[model_number], std_halomass_array[model_number], np.mean(halomass_nonlog), np.std(halomass_nonlog), len(mass_central[model_number][snapshot_idx]))

		## If want to do mean/etc of halo mass need to update script. ##
		bin_middle_array[model_number].append(bin_middle)

		'''
		if rank == 0:
			mean_ngammafesc_array[model_number][snapshot_idx] = np.nan_to_num(mean_ngammafesc_array[model_number][snapshot_idx])
			std_ngammafesc_array[model_number][snapshot_idx] = np.nan_to_num(std_ngammafesc_array[model_number][snapshot_idx])
			outfile = '/lustre/projects/p004_swin/jseiler/tiamat/mean_mvir_ngammafesc_z%.3f.dat' %(AllVars.SnapZ[SnapList[model_number][snapshot_idx]])
			with open(outfile, 'w') as f:
				np.savetxt(outfile, mean_ngammafesc_array[model_number][snapshot_idx])
			print "Saved %s" %(outfile)


			outfile = '/lustre/projects/p004_swin/jseiler/tiamat/std_mvir_ngammafesc_z%.3f.dat' %(AllVars.SnapZ[SnapList[model_number][snapshot_idx]])
			with open(outfile, 'w') as f:
				np.savetxt(outfile, std_ngammafesc_array[model_number][snapshot_idx])
			print "Saved %s" %(outfile)
		'''
	#std_ngammafesc_array[model_number] = 0.434 * np.divide(std_ngammafesc_array[model_number], mean_ngammafesc_array[model_number])
	
	#mean_ngammafesc_array[model_number] = np.log10(~np.isnan(mean_ngammafesc_array[model_number]))
	
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
			
			mean = np.log10(mean_ngammafesc_array[model_number][snapshot_idx])
			std = 0.434 * np.divide(std_ngammafesc_array[model_number][snapshot_idx], mean_ngammafesc_array[model_number][snapshot_idx])
			bin_middle = bin_middle_array[model_number][snapshot_idx]



			ax1.plot(bin_middle, mean, color = colors[snapshot_idx], linestyle = linestyles[model_number], rasterized = True, label = title)
			#ax1.scatter(mean_halomass_array[model_number][snapshot_idx], np.mean(~np.isnan(mean)), color = colors[snapshot_idx], marker = 'o', rasterized = True, s = 40, lw = 3)	
			if (len(SnapList) == 1):
    				ax1.fill_between(bin_middle, np.subtract(mean,std), np.add(mean,std), color = colors[snapshot_idx], alpha = 0.25)

	ax1.set_xlabel(r'$\log_{10}\ M_{\mathrm{vir}}\ [M_{\odot}]$', size = talk_fontsize) 
    	ax1.set_ylabel(r'$\dot{N}_\gamma \: f_\mathrm{esc} \: [\mathrm{s}^{-1}]$', size = talk_fontsize)
    	ax1.set_xlim([8.5, 12])
    	#ax1.set_ylim([1e50, 1e54])   

    	ax1.xaxis.set_minor_locator(mtick.MultipleLocator(0.1))
#    	ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
 	
#    	ax1.set_yscale('log', nonposy='clip')
#    	for model_number in xrange(0, len(SnapList)):
#		ax1.plot(1e100, 1e100, color = 'k', ls = linestyles[model_number], label = model_tags[model_number], rasterized=True)
	
	
    	leg = ax1.legend(loc='upper left', numpoints=1, labelspacing=0.1)
    	leg.draw_frame(False)  # Don't want a box frame
    	for t in leg.get_texts():  # Reduce the size of the text
        	t.set_fontsize('medium')

	outputFile = './' + output_tag + output_format
    	plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
    	print 'Saved file to', outputFile
    	plt.close()




'''
At this point each process has it's own slice of the data set.  
Each of these slices has inputs at each of the redshifts and for each model considered.
What we do here is calculate the mean and standard deviation for each processor and then pool it and pass back to the master task.
'''

def plot_photoncount(SnapList, ngamma, fesc, halo_mass, num_files, model_tags, output_tag): 

    print "Plotting a count of the photons."

    sum_array = []

    for model_number in xrange(0, len(SnapList)): 
	sum_array.append([])

    	for snapshot_idx in xrange(0, len(SnapList[model_number])):
		
		w = np.array(np.where(((halo_mass[model_number][snapshot_idx] < kink_low) & (halo_mass[model_number][snapshot_idx] > m_low))| ((halo_mass[model_number][snapshot_idx] > kink_high) & (halo_mass[model_number][snapshot_idx] < m_high)))[0])

		tmp = np.array(ngamma[model_number][snapshot_idx])
		tmp2 = np.array(fesc[model_number][snapshot_idx])
#		ngamma_cut = ngamma[model_number][snapshot_idx][w]
		ngamma_cut = tmp[w] 
		fesc_cut = tmp2[w] 
		ngammafesc = np.multiply(np.power(10, ngamma_cut), fesc_cut)
	    

		sum_local = np.sum(ngammafesc)


		if rank == 0:
			sum_total = np.zeros_like(sum_local)
		else:
			sum_total = None
		comm.Reduce([sum_local, MPI.DOUBLE], [sum_total, MPI.DOUBLE], op = MPI.SUM, root = 0)

		if rank == 0:
			sum_array[model_number].append(sum_total / (pow(AllVars.BoxSize / AllVars.Hubble_h,3) * (float(num_files) / float(125.0))))
	

    if (rank == 0):
    	ax1 = plt.subplot(111)

	for model_number in xrange(0, len(SnapList)):
    		t = np.empty(len(SnapList[model_number]))
		for snapshot_idx in xrange(0, len(SnapList[model_number])):
       			t[snapshot_idx] = (t_BigBang - cosmo.lookback_time(AllVars.SnapZ[SnapList[model_number][snapshot_idx]]).value) * 1.0e3   
	
    		ax1.plot(t, np.log10(sum_array[model_number]), color = colors[model_number], ls = linestyles[model_number], label = model_tags[model_number], lw = 4)  
    		#ax1.fill_between(t, np.subtract(mean,std), np.add(mean,std), color = colors[model_number], alpha = 0.25)

    	#ax1.xaxis.set_minor_locator(mtick.MultipleLocator(time_tick_interval))
    	#ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.025))
    	ax1.set_xlim(time_xlim)
    	ax1.set_ylim([48.5, 51.5])

    	ax2 = ax1.twiny()

    	t_plot = (t_BigBang - cosmo.lookback_time(z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
    	z_labels = ["$%d$" % x for x in z_plot] # Properly Latex-ize the labels.

    	ax2.set_xlabel(r"$z$", size = label_size)
    	ax2.set_xlim(time_xlim)
    	ax2.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
    	ax2.set_xticklabels(z_labels) # But label them as redshifts.

    	ax1.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = label_size)
    	ax1.set_ylabel(r'$\sum f_\mathrm{esc}\dot{N}_\gamma \: [\mathrm{s}^{-1}\mathrm{Mpc}^{-3}]$', fontsize = talk_fontsize)

	leg = ax1.legend(loc='lower right', numpoints=1, labelspacing=0.1)
    	leg.draw_frame(False)  # Don't want a box frame
    	for t in leg.get_texts():  # Reduce the size of the text
        	t.set_fontsize(talk_legendsize)

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

#	ax1.text(0.075, 0.965, '(a)', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
	ax1.text(350, 50.0, r"$68\%$", horizontalalignment='center', verticalalignment = 'center', fontsize = label_size - 2)
	ax1.text(350, 50.8, r"$95\%$", horizontalalignment='center', verticalalignment = 'center', fontsize = label_size - 2)

    	plt.tight_layout()
    	outputFile = './' + output_tag + output_format
    	plt.savefig(outputFile)  # Save the figure
    	print 'Saved file to', outputFile
    	plt.close()




## Calculates the stellar mass for galaxies with host halos that have a certain amount of particles.

def Calculate_HaloPartStellarMass(HaloPart, StellarMass, BoundLow, BoundHigh):

#	print HaloPart
	print 
	w = np.where((HaloPart > BoundLow) & (HaloPart < BoundHigh))[0]

	Mass = np.mean(10**(StellarMass[w]))
	Mass_std = np.std(10**(StellarMass[w]))
	count = len(StellarMass[w])
#	print "Mass %.4e \t Mass_std = %.4e \t Count = %d" %(Mass, Mass_std, count) 

	return Mass, Mass_std


#################################

 
'''
countSnap_low = 30
countSnap_high = 99

idx_low = [i for i in range(len(G_MySim.GridHistory)) if G_MySim.GridHistory[i,countSnap_low] != -1] # Indices for all galaxies that exist at snapshot 'countSnap_low'.
numgals_low = len(idx_low)

numgals_high = len(G_MySim.GridHistory[G_MySim.GridHistory[idx_low, countSnap_high] != -1])
 

print numgals_low
print numgals_high
'''

HaloPart_Low = 41 # Bounds for where we define the cutoff for a 'Dark Matter Halo'. Below this we can't be sure of the results.
HaloPart_High = 51

calculate_observed_LF = 0

### The arrays in this block control constants for each model ##

number_models = 4

'''
galaxies_model1 = '/lustre/projects/p004_swin/jseiler/tiamat/lowz_z4.929'
galaxies_model2 = '/lustre/projects/p004_swin/jseiler/tiamat/lowz_z4.929'
galaxies_model3 = '/lustre/projects/p004_swin/jseiler/tiamat/lowz_z4.929'
galaxies_model4 = '/lustre/projects/p004_swin/jseiler/tiamat/lowz_z4.929'
'''
galaxies_model1 = '/lustre/projects/p004_swin/jseiler/18month/IRA_z5.000'
galaxies_model2 = '/lustre/projects/p004_swin/jseiler/18month/IRA_z5.000'
galaxies_model3 = '/lustre/projects/p004_swin/jseiler/18month/IRA_z5.000'
galaxies_model4 = '/lustre/projects/p004_swin/jseiler/18month/IRA_z5.000'

'''
merged_galaxies_model1 = '/lustre/projects/p004_swin/jseiler/tiamat/lowz_MergedGalaxies'
merged_galaxies_model2 = '/lustre/projects/p004_swin/jseiler/tiamat/lowz_MergedGalaxies'
merged_galaxies_model3 = '/lustre/projects/p004_swin/jseiler/tiamat/lowz_MergedGalaxies'
merged_galaxies_model4 = '/lustre/projects/p004_swin/jseiler/tiamat/lowz_MergedGalaxies'
'''

merged_galaxies_model1 = '/lustre/projects/p004_swin/jseiler/18month/IRA_MergedGalaxies'
merged_galaxies_model2 = '/lustre/projects/p004_swin/jseiler/18month/IRA_MergedGalaxies'
merged_galaxies_model3 = '/lustre/projects/p004_swin/jseiler/18month/IRA_MergedGalaxies'
merged_galaxies_model4 = '/lustre/projects/p004_swin/jseiler/18month/IRA_MergedGalaxies'

galaxies_filepath_array = [galaxies_model1, galaxies_model2, galaxies_model3, galaxies_model4]
merged_galaxies_filepath_array = [merged_galaxies_model1, merged_galaxies_model2, merged_galaxies_model3, merged_galaxies_model4]

#number_snapshots = [164, 164, 164, 164] # Property of the simulation.
number_snapshots = [101, 101, 101, 101] # Property of the simulation.
#number_snapshots = [101] # Property of the simulation.
FirstFile = [0, 0, 0, 0]
LastFile = [124, 124, 124, 124]

#model_tags = [r"Tiamat"]
#model_tags = [r"$f_\mathrm{esc} = \mathrm{Constant}$", r"$f_\mathrm{esc} \: \propto \: M_\mathrm{H}^{-1}$", r"$f_\mathrm{esc} \: \propto \: M_\mathrm{H}$", r"$f_\mathrm{esc} \: \propto \: f_\mathrm{ej}$"]
model_tags = [r"$f_\mathrm{esc} = 0.50$", r"$f_\mathrm{esc} \: \propto \: M_\mathrm{H}^{-1}$", r"$f_\mathrm{esc} \: \propto \: M_\mathrm{H}$", r"$f_\mathrm{esc} \: \propto \: f_\mathrm{ej}$"]
#model_tags = [r"No Cut", r"10 Particle Cut", r"100 Particle Cut"]

for model_number in xrange(0,number_models):
	assert(LastFile[model_number] - FirstFile[model_number] + 1 >= size)

## Constants used for each model. ##
# Need to add an entry for EACH model. #

sSFR_min = [1.0e100, 1.0e100, 1.0e100, 1.0e100]
sSFR_max = [-1.0e100, -1.0e100, 1.0e100, 1.0e100]
halo_cut = [100, 100, 100, 100]
source_efficiency = [1, 1, 1, 1]

fesc_lyman_alpha = [0.3, 0.3, 0.3, 0.3]
fesc_prescription = [0, 1, 1, 2] # 0 is constant, 1 is scaling with halo mass, 2 is scaling with ejected fraction.

## Normalizations for the escape fractions. ##
# For prescription 0, requires a number that defines the constant fesc.
# For prescription 1, fesc = A*M^B. Requires an array with 2 numbers the first being A and the second B.
# For prescription 2, fesc = A*fej + B.  Requires an array with 2 numbers the first being A and the second B.
fesc_normalization = [0.50, [6000.0, -0.4], [5.0e-5, 0.375], [1.0, 0.0]] 
##

SnapList = [np.arange(100, 20, -1), np.arange(100, 20, -1), np.arange(100, 20, -1), np.arange(100, 20, -1)]
# z = [6, 7, 8] are snapshots [78, 64, 51]

simulation_norm = [0, 0, 0, 0] # 0 for MySim, 1 for Mini-Millennium, 2 for Tiamat (up to z =5), 3 for extended Tiamat (down to z = 1.6ish).

galaxy_halo_mass_lower = [95, 95, 95, 95] # These limits are for the number of particles in a halo.  
galaxy_halo_mass_upper = [105, 105, 105, 105] # We calculate the average stellar mass for galaxies whose host halos have particle count between these limits.

##############################################################################################################

## Halo Initialization ##

w_halo = []
mass_halo = []
photons_HI_halo = []
photons_HI_tot_halo = []
source_efficiency_halo = 10

## General Array Initialization ##

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

for model_number in xrange(0, number_models):

	w_gal.append([])
	mass_gal.append([])
	photons_HI_gal.append([])
	photons_HI_tot_gal.append([])
	SFR_gal.append([])
	sSFR_gal.append([])
        metallicity_gal.append([])
	metallicity_tremonti_gal.append([])
	halo_part_count.append([])

	mass_central.append([])
	photons_HI_central.append([])
	photons_HI_tot_central.append([])

	Mfilt_gnedin.append([])
	Mfilt_sobacchi.append([])

	lyman_alpha.append([])

	L_UV.append([])
	M_UV.append([])

	M_UV_Obs.append([])
	mean_A.append([])

	ejected_fraction.append([])

	fesc_local.append([])

	halo_count.append([])

	galaxy_halo_mass_mean.append([])
	galaxy_halo_mass_std.append([])

dilute_galaxies = 100000
for model_number in xrange(0, number_models):

    	GG, Gal_Desc = ReadScripts.ReadGals_SAGE_DelayedSN(galaxies_filepath_array[model_number], FirstFile[model_number], LastFile[model_number], number_snapshots[model_number], comm)
    	G_Merged, Merged_Desc = ReadScripts.ReadGals_SAGE_DelayedSN(merged_galaxies_filepath_array[model_number], FirstFile[model_number], LastFile[model_number], number_snapshots[model_number], comm)
	#G = GG
	G = ReadScripts.Join_Arrays(GG, G_Merged, Gal_Desc)

	if(simulation_norm[model_number] == 0):
		AllVars.Set_Params_Mysim()
	elif(simulation_norm[model_number] == 3):
		AllVars.Set_Params_Tiamat_extended()
	else:
		print "Simulation norm was set to %d.  This option hasn't been implemented yet." %(simulation_norm[model_number])

	current_fesc_lyman_alpha = fesc_lyman_alpha[model_number]
	current_source_efficiency = source_efficiency[model_number]
	current_halo_cut = halo_cut[model_number]

	def calculate_fesc(fesc_prescription, halo_mass, ejected_fraction, fesc_normalization):
		print "calculating fesc with prescription %d" %(fesc_prescription) 
		if (fesc_prescription == 0):
			fesc = np.full((len(halo_mass)), fesc_normalization)
		elif (fesc_prescription == 1):
			fesc = fesc_normalization[0]*pow(10,halo_mass*fesc_normalization[1])
		elif (fesc_prescription == 2):
			fesc = fesc_normalization[0]*ejected_fraction + fesc_normalization[1]
		'''
		fesc = fesc[~np.isnan(fesc)] # Get rid of any lingering Nans.
		fesc = fesc[fesc < 1] # Get rid of any values above 1. 
		fesc = fesc[fesc > 0] # Get rid of any values below 0. 
		'''

		fesc[np.isnan(fesc)] = 0 # Get rid of any lingering Nans.
		fesc[fesc > 1] = 1.0  # Get rid of any values above 1. 
		fesc[fesc < 0] = 0.0  # Get rid of any values below 0. 
		return fesc

	def calculate_photons(SFR, Z):

		Ngamma_HI = []

		for i in xrange(0, len(SFR)):
			Ngamma_HI_tmp = 0.0
			if (SFR[i] == 0):
				Ngamma_HI_tmp = 0.0	
			elif (Z[i] < 0.0025):
				Ngamma_HI_tmp = SFR[i] + 53.354
			elif (Z[i] >= 0.0025 and Z[i] < 0.006):
				Ngamma_HI_tmp = SFR[i] + 53.290
			elif (Z[i] >= 0.006 and Z[i] < 0.014):
				Ngamma_HI_tmp = SFR[i] + 53.240
			elif (Z[i] >= 0.014 and Z[i] < 0.30):
				Ngamma_HI_tmp = SFR[i] + 53.166
			else:
				Ngamma_HI_tmp = SFR[i] + 53.041 
			if (SFR[i] != 0):
				assert(Ngamma_HI_tmp > 0.0)	
			Ngamma_HI.append(Ngamma_HI_tmp)
		return Ngamma_HI


	for snapshot_idx in xrange(0, len(SnapList[model_number])):
   		 
		current_snap = SnapList[model_number][snapshot_idx]

		w_gal[model_number].append(np.where((G.GridHistory[:, current_snap] != -1) & (G.GridStellarMass[:, current_snap] > 0.0) & (G.GridStellarMass[:, current_snap] < 1e5) & (G.GridCentralGalaxyMass[:, current_snap] < 1e5) & (G.GridSFR[:, current_snap] > 0.0) & (G.LenHistory[:, current_snap] > current_halo_cut))[0])
		current_idx = w_gal[model_number][snapshot_idx]
	
		print "There were %d galaxies for snapshot %d (Redshift %.4f) model %d." %(len(current_idx), current_snap, AllVars.SnapZ[current_snap], model_number)

      		halo_count[model_number].append(G.LenHistory[current_idx, current_snap])
		mass_gal[model_number].append(np.log10(G.GridStellarMass[current_idx, current_snap] * 1.0e10 / AllVars.Hubble_h))

		SFR_gal[model_number].append(np.log10(G.GridSFR[current_idx, current_snap])) # Msun yr^-1.  Log Units.

		halo_part_count[model_number].append(G.LenHistory[current_idx, current_snap])
		metallicity_gal[model_number].append(G.GridZ[current_idx, current_snap])	
		metallicity_tremonti_gal[model_number].append(np.log10(G.GridZ[current_idx, current_snap] / 0.02) + 9.0)

		mass_central[model_number].append(np.log10(G.GridCentralGalaxyMass[current_idx, current_snap] * 1.0e10 / AllVars.Hubble_h))
	
		#Photon_Factor = AllVars.Solar_Mass * current_source_efficiency * AllVars.BaryonFrac / AllVars.Ionized_Mass_H / AllVars.Proton_Mass / (AllVars.Lookback_Time[current_snap] - AllVars.Lookback_Time[current_snap+1]) / AllVars.Sec_Per_Megayear / 1.0e3
		#photons_HI_central[model_number].append(np.log10(10**(mass_central[model_number][snapshot_idx]) * Photon_Factor))
		#photons_HI_tot_central[model_number].append(np.log10(sum(10**photons_HI_central[model_number][snapshot_idx])))

		Mfilt_gnedin[model_number].append(G.MfiltGnedin[current_idx, current_snap])
		Mfilt_sobacchi[model_number].append(G.MfiltSobacchi[current_idx, current_snap])
		
#		lyman_alpha[model_number].append(np.log10(0.68*(1.0 - current_fesc) * current_fesc_lyman_alpha * (AllVars.LymanAlpha_Energy * AllVars.eV_to_erg)) + photons_HI[model_number]



		L_UV[model_number].append(SFR_gal[model_number][snapshot_idx] + 39.927)# Using relationship from STARBURST99, units of erg s^-1 A^-1. Log Units.
		M_UV[model_number].append(AllVars.Luminosity_to_ABMag(L_UV[model_number][snapshot_idx], 1600))

		if (calculate_observed_LF == 1):

	      		M_UV_bins = np.arange(-24, -16, 0.1)
	      		A_Mean = np.zeros((len(MUV_bins)))	

			for j in xrange(0, len(M_UV_bins)):
				beta = calculate_beta(M_UV_bins[j], AllVars.SnapZ[current_snap]) 
				dist = np.random.normal(beta, 0.34, 10000)
				A = 4.43 + 1.99*dist
				A[A < 0] = 0
		
				A_Mean[j] = np.mean(A)


			indices = np.digitize(M_UV[model_number][snap_idz], M_UV_bins) # Bins the simulation magnitude into the MUV bins. Note that digitize defines an index i if bin[i-1] <= x < bin[i] whereas I prefer bin[i] <= x < bin[i+1]
			dust = A_Mean[indices]
			Flux = AllVars.Luminosity_to_Flux(L_UV[model_number][snapshot_idx], 10.0) # Calculate the flux from a distance of 10 parsec, units of erg s^-1 A^-1 cm^-2.  Log Units.  
			Flux_Observed = Flux - 0.4*dust
	
			f_nu = ALlVars.spectralflux_wavelength_to_frequency(10**Flux_Observed, 1600) # Spectral flux desnity in Janksy.
			M_UV_obs[model_number].append(-2.5 * np.log10(f_nu) + 8.90) # AB Magnitude from http://www.astro.ljmu.ac.uk/~ikb/convert-units/node2.html

		ejected_fraction[model_number].append(G.EjectedFraction[current_idx, current_snap])
	
		fesc_local[model_number].append(calculate_fesc(fesc_prescription[model_number], mass_central[model_number][snapshot_idx], ejected_fraction[model_number][snapshot_idx], fesc_normalization[model_number])) 
		
		tmp_mean, tmp_std = Calculate_HaloPartStellarMass(halo_part_count[model_number][snapshot_idx], mass_gal[model_number][snapshot_idx], galaxy_halo_mass_lower[model_number], galaxy_halo_mass_upper[model_number])

		galaxy_halo_mass_mean[model_number].append(tmp_mean)
		galaxy_halo_mass_std[model_number].append(tmp_std)
		 
		photons_HI_gal[model_number].append(calculate_photons(SFR_gal[model_number][snapshot_idx], metallicity_gal[model_number][snapshot_idx])) 


#HaloMassFunction(SnapList, mass_central, simulation_norm, model_tags, "HaloMassFunction_Kink")
#StellarMassFunction(SnapList, mass_gal, simulation_norm, galaxy_halo_mass_mean, model_tags, "18month_SMF") ## PARALLEL COMPATIBLE
#plot_ejectedfraction(SnapList, mass_central, ejected_fraction, model_tags, "18month_Ejected") ## PARALELL COMPATIBLE # Ejected fraction as a function of Halo Map 
plot_fesc(SnapList, fesc_local, mass_gal, mass_central, model_tags, "18month_fesc_diffprescription") ## PARALELL COMPATIBLE 
plot_photoncount(SnapList, photons_HI_gal, fesc_local, mass_central, 125, model_tags, "18month_nion_diffprescription") ## PARALELL COMPATIBLE

#plot_mvir_Ngamma(SnapList, mass_central, photons_HI_gal, fesc_local, model_tags, "Mvir_Ngamma_1file") # Started parallel compatibility. 
#plot_mvir_fesc(SnapList, mass_central, fesc_local, model_tags, "Mvir_fesc")   # Stared parallel compatibility

