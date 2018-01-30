#!/usr/bin/env python
from __future__ import print_function

import sys
print(sys.path)


import numpy as np
from numpy import *
np.set_printoptions(threshold = np.nan, linewidth = 1000000)

import matplotlib
matplotlib.use('Agg')

import os

import pylab as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy as np
from random import sample, seed
from os.path import getsize as getFileSize
import math
import random
import csv
from io import StringIO
from collections import Counter
from matplotlib.colors import LogNorm
import time
from scipy.ndimage.filters import generic_filter as gf
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec
import scipy.integrate as integrate
from numpy.fft import fftn, ifftn
import matplotlib.patheffects as PathEffects
from astropy import units as u
from astropy import cosmology
import matplotlib.ticker as mtick
import itertools
from matplotlib import patches

from mpi4py import MPI

import PlotScripts
import ReadScripts
import AllVars




comm= MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
AllVars.Set_Constants()

AllVars.Set_Constants()
PlotScripts.Set_Params_Plot()
matplotlib.rcdefaults()
plt.rc('axes', color_cycle=['k', 'b', 'r', 'g', 'm', 'c','y', '0.5'], labelsize='x-large')
plt.rc('lines', linewidth='2.0')
# plt.rc('font', variant='monospace')
plt.rc('legend', numpoints=1, fontsize='x-large')
plt.rc('text', usetex=True)

label_size = 20 # 12 for paper.
extra_size = 2
legend_size = 16 # 20 for paper.
paper_lineweight = 3
plt.rc('xtick', labelsize=label_size)
plt.rc('ytick', labelsize=label_size)
plt.rc('text', usetex=True)
tick_interval = 0.25
time_tick_interval = 25
z_plot = np.arange(6, 14)  #Range of redshift we wish to plot.

cosmo = AllVars.Set_Params_Mysim()

time_xlim = [315, 930]
#time_xlim = [round((AllVars.t_BigBang - cosmo.lookback_time(z_plot[-1]).value) * 1.0e3), round((AllVars.t_BigBang - cosmo.lookback_time(z_plot[0]).value) * 1.0e3)]

time_subplot_label = 350 # Location (in 'Time Since Big Bang [Myr^-1]) of subplot identifier (a), (b), (c) etc.

cut_slice = 44

#cosmo = cosmology.FlatLambdaCDM(H0 = AllVars.Hubble_h*100, Om0 = AllVars.Omega_m) 
#t_BigBang = cosmo.lookback_time(100000).value # Lookback time to the Big Bang in Gyr.
output_format = ".pdf"

def linear_growth(z, OM, OL):
        D = (OM*((1+z)**3) / (OM*(1+z) - (OM + OL - 1)*((1+z)**2) + OM))**(4/7)
        return D

def cosmo_density(z, Hubble_h):
        rho = (3*100*100*Hubble_h*Hubble_h*AllVars.km_to_m*AllVars.km_to_m)/(8*np.pi*AllVars.G*AllVars.mpc_to_m*AllVars.mpc_to_m)*(1+z)**3
        return rho

## Number density of hydrogen in m^-3##
def n_HI(z, Hubble_h, OB, Y):
	n_HI = cosmo_density(z, Hubble_h)*OB*(1-Y)/(AllVars.Proton_Mass/1000.0) 
	return n_HI


def stromgren_sphere(Q, T, z, Y):
	# Takes the photon flux (Q, photons per second) with an electron temperature (T, Kelvin) and helium fraction (Y)
	# And returns the radius of a Stromgren sphere at redshift z in Mpc/h.

	Beta = 2.e-16*(T_e)**(-3./4.) # Recombination rate (m^3 s^-1).
	num_HI = n_HI(z, Hubble_h, OB, Y) # Number density of hydrogen (m^3).

	R = ((3.*Q)/(4.*np.pi*num_HI*num_HI*Beta)) ** (1./3.)
	R /= Mpc_m # Stromgren radius in units of Mpc.
	R *= Hubble_h

	print("Stromgren Radius for a flux of %.3e photons per second at a redshift %.3f and an electron temperature of %.3e is %.3e Mpc/h" %(Q, z, T_e, R))
	return R


def T_naught(z, OB, h, OM):
	# T0 brightness temperature from Hutter et al (2016)
	# Output in units of mK.

	T0 = 28.5 * ((1.0+z)/10.0)**(0.5) * OB/0.042 * h/0.73 * (OM/0.24)**(-0.5) 
	return T0

def Hubble_Param(z, h, OM):

    ## Output in units of km/s/Mpc

	H0 = h*100
	H = H0**2 * (OM * (1+z)**3 + (1 - OM))

	return np.sqrt(H)


def plot_single(z, ionized_cells, Ncells, simulation_norm, OutputDir, output_tag):
# Plots the ionization bubbles for a single set of cells.

## Input ##
# z is the redshift we're plotting at. 
# ionized_cells is the array that contains the ionized cells, OutputDir is the path of the output directory.

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*_*z*.*output_format*'


	print("Plotting ionized bubbles.")
	
	ionized_cells = np.log10(1 - ionized_cells)
	#ionized_cells = (1 - ionized_cells)

	ax = plt.subplot(111)

	im = ax.imshow(ionized_cells[:,:,cut_slice:cut_slice+1].mean(axis = -1), interpolation='nearest', origin='low', extent =[0,AllVars.BoxSize,0,AllVars.BoxSize], vmin = -8, vmax = 0, cmap = 'afmhot_r')

	cbar = plt.colorbar(im, ax = ax)
	cbar.set_label(r'$\mathrm{log}_{10}\left(x_\mathrm{HI}\right)$')
				    
	ax.set_xlabel(r'$\mathrm{x}  (h^{-1}Mpc)$')  
	ax.set_ylabel(r'$\mathrm{y}  (h^{-1}Mpc)$')  

	ax.set_xlim([0.0, AllVars.BoxSize]) 
	ax.set_ylim([0.0, AllVars.BoxSize])

	title = r"z = %.3f, $\langle x_{HI} = %.3f \rangle$" %(z, calculate_volume_frac(pow(10,ionized_cells), Ncells))
	ax.set_title(title)

	#outputFile = OutputDir + output_tag + '_z' + str(z) + output_format 
	outputFile = OutputDir + output_tag + output_format 
	plt.savefig(outputFile)  # Save the figure
	print('Saved file to {0}'.format(outputFile))
			
	plt.close()
	
############
'''
def plot_comparison(z, ionized_cells_model1, ionized_cells_model2, nion_model1, nion_model2, OutputDir, model_tags, output_tag):
# Plots ionization bubbles for two realizaiton of the box.  These bubbles will be overlaid on top of each other to allow comparison between the two.

## Input ##
# z is the redshift we're plotting at.
# ionized_cells_model1 (and ionized_cells_model2) are the arrays that contains the ionized cells.
# OutputDir is the path of the output directory.

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*_*z*.*output_format*'

	print("")
	print("Plotting a comparison between the ionized bubbles.")	

	master_cells = np.empty((GridSize,GridSize,GridSize))

	# This block of for loops will decide the stacking of the arrays.
	# In the current implemenation there layers are ionization sources > first ionized_cells > second ionized_cells > no ionization.
	# To change layering take the if statement block and move it up or down.
	for p in range(0, GridSize):
		for l in range(0, GridSize):
			for m in range(0,GridSize):
				if(nion_model1[p,l,m] != 0): 
					master_cells[p,l,m] = -1
				elif(nion_model2[p,l,m] != 0):
					master_cells[p,l,m] = -2
				elif(ionized_cells_model1[p,l,m] == 1):
					master_cells[p,l,m] = 1
				elif(ionized_cells_model2[p,l,m] == 1):
					master_cells[p,l,m] = 0
				else: 
					master_cells[p,l,m] = 2

	master_cells[master_cells < -2] = 2

	ax = plt.subplot(111)

	cmap = cm.get_cmap('afmhot_r', 5)
		
	im = ax.imshow(master_cells[:,:,cut_slice:cut_slice+1].mean(axis = -1), interpolation='none', origin='low', extent =[0,BoxSize,0,BoxSize], vmin = -2, vmax = 2, cmap = cmap)

	cbar = plt.colorbar(im, ax = ax, ticks = [1.6, 0.85, 0.0, -0.85, -1.6]) # Set the tick intervals for the colourbar.
	cbar.ax.set_yticklabels(['None', model_tags[0], model_tags[1], "Sources - " + model_tags[0], "Sources - " + model_tags[1]]) # Colourbar labels.
	cbar.ax.tick_params(axis = 'y', which = 'both', length = 0) # Alter how the ticks appear.
	ax.set_xlabel(r'$\mathrm{x}  (h^{-1}Mpc)$')  
	ax.set_ylabel(r'$\mathrm{y}  (h^{-1}Mpc)$')  

	ax.set_xlim([0.0, BoxSize]) 
	ax.set_ylim([0.0, BoxSize])

	outputFile = OutputDir + output_tag + output_format 	
	plt.savefig(outputFile)  # Save the figure
	print('Saved file to {0}'.format(outputFile))
			
	plt.close()
'''
##############

def plot_sources(z, filepath, nion_tmp, OutputDir, output_tag): 
# Plots the location of the ionizing sources.  At the moment this function only plots if there is a source in the cell and not its photon count.

## Input ##
# z is the redshift we're plotting at. nion is the array that contains the ionizing photons

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*_*z*.*output_format*'

	print("")
	print("Plotting ionizing sources for file {0}".format(filepath))
	#nion_tmp[nion_tmp > 0] = 1 # Simply set to measure if a cell contains a photon source, not its strength.

	ax = plt.subplot(111)

	im = ax.imshow(nion_tmp[:,:,cut_slice:cut_slice+1].mean(axis = -1), interpolation = 'none', origin = 'low', extent = [0,BoxSize,0,BoxSize], vmin = 0, vmax = np.amax(nion_tmp[:,:,cut_slice:cut_slice+1].mean(axis=-1)), cmap = 'afmhot_r')

	cbar = plt.colorbar(im, ax = ax)
	cbar.set_label(r'Photons Emitted in Cell')
				    
	ax.set_xlabel(r'$\mathrm{x}  (h^{-1}Mpc)$')  
	ax.set_ylabel(r'$\mathrm{y}  (h^{-1}Mpc)$')  

	ax.set_xlim([0.0, BoxSize]) 
	ax.set_ylim([0.0, BoxSize])
	
	outputFile = OutputDir + output_tag + '_z' + str(z) + output_format 		
	plt.savefig(outputFile)  # Save the figure
	print('Saved file to {0}'.format(outputFile))	
			
	plt.close()
		
###########

def difference_map(z, ionized_cells_model1, ionized_cells_model2, OutputDir, model_tags, output_tag):
# Plots a difference map between two ionization realizations. This is function compliments 'plot_comparison'.

## Input ##
# z is the redshift we're plotting at.
# ionized_cells_model1 (and ionized_cells_model2) are the arrays that contains the ionized cells.
# OutputDir is the path of the output directory.

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*_*z*.*output_format*'

	print("")
	print("Plotting a difference map between the models.")

	ax = plt.subplot(111)

	cmap = cm.get_cmap('RdGy', 3)
	im = ax.imshow(ionized_cells_model1[:,:,cut_slice:cut_slice+1].mean(axis = -1) - ionized_cells_model2[:,:,cut_slice:cut_slice+1].mean(axis = -1), interpolation = 'none', origin = 'low', extent = [0,BoxSize,0,BoxSize], vmin = -1, vmax = 1, cmap = cmap)

	cbar = plt.colorbar(im, ax = ax, ticks = [0.65, 0.0, -0.65]) # Set the tick intervals for the colourbar.
	cbar.ax.set_yticklabels([model_tags[0], "Agreement", model_tags[1]]) 
	cbar.ax.tick_params(axis = 'y', which = 'both', length = 0) # Alter how the ticks appear.	
				    
	ax.set_xlabel(r'$\mathrm{x}  (h^{-1}Mpc)$')  
	ax.set_ylabel(r'$\mathrm{y}  (h^{-1}Mpc)$')  

	ax.set_xlim([0.0, BoxSize]) 
	ax.set_ylim([0.0, BoxSize])

	outputFile = OutputDir + output_tag + output_format 		
	plt.savefig(outputFile)  # Save the figure
	print('Saved file to {0}'.format(outputFile))	
				
	plt.close()	

###########

def calculate_volume_frac(ionized_cells, Ncell):
# This function will calculate (and return!) the fraction of cells that are ionized.  
# This function was separated from 'plot_global_frac' as the global fractions for each redshifts must be calculated before we can plot them.

## Input ##
# ionized_cells is the array that contains the ionization data.

## Output ##
# The fraction of HI in the grid.

#	HII = 1 - len(ionized_cells[ionized_cells > 0.95])/float(Ncell**3)
    
    assert(len(np.argwhere(np.isnan(ionized_cells))) == 0)
           
    HII = np.mean(ionized_cells) 
    print("There is {0:.4e} fraction of HII.".format(HII))

    return HII

###########


def calculate_mass_frac(ionized_cells, density):
# This function will calculate (and return!) the mass-averaged fraction of cells that are ionized.  
# This function was separated from 'plot_global_frac' as the global fractions for each redshifts must be calculated before we can plot them.

## Input ##
# ionized_cells is the array that contains the ionization data.
# density is the array that contains the density data.

## Output ##
# The mass-averaged fraction of HI in the grid.

    HI = 1 - sum(ionized_cells * density/float(sum(density)))  

    print("")
    print("Mass averaged HII fraction is {0:.4f}".format(HI))

    return HI

###########


def calculate_total_nion(model_name, nion):
# This function will calculate (and return!) the total number of ionizing photons (per second) in the grid.
# This function was separated from 'plot_total_nion' as total nion for each redshifts must be calculated before we can plot them.
# ionized_cells is the array that contains the ionization data.

# Output: The total number of ionizing photons (per second) in the grid.

#        nion_total = sum(nion)/Hubble_h/(BoxSize**3/Hubble_h**3)

	nion_total = sum(nion)

	print("")
	print("The total number of ionizing photons (per second) for the {0} Model is {1:.4e} ({2:.4e} Mpc^-3).".format(model_name, nion_total, nion_total / pow(AllVars.BoxSize/AllVars.Hubble_h, 3)))

	return nion_total


##########




	



#########


def plot_global_frac(ZZ, mass_frac, volume_frac, plot_time, labels, OutputDir, output_tag):
# This function will plot the global HI fraction as a function of redshift. The function accepts an array of arrays for the global HI fraction to allow plotting of multiple realizations on the same graph.

## Input ##
# ZZ is an array containing the redshifts we plot over.
# global_frac is an array of arrays where each array contains the global HI fraction at each redshift.
# labels is an array that contains the labels for each realization we are plotting.

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*_*z*.*output_format*'


    print("")
    print("Plotting the global fraction of HI.")

    legend_labels = []

    fig = plt.figure()

    ax1 = plt.subplot2grid((3,1), (0,0), rowspan =3)
	#ax3 = plt.subplot2grid((3,1), (2,0))

    t = np.empty(len(volume_frac[0]))
	
    for i in range(0, len(volume_frac[0])):
        t[i] = (AllVars.t_BigBang - cosmo.lookback_time(ZZ[i]).value) * 1.0e3 	

    for i in range(0, len(volume_frac)):

        '''
        dt = (cosmo.lookback_time(MC_ZZ[i][0]).value - cosmo.lookback_time(MC_ZZ[i][-1]).value) * 1.0e3
        if (i != 1):
            tmp = r"%s, $\Delta t = %.2f \: \mathrm{Myr}$" %(labels[i], dt)
        else:
            tmp = r"%s, $\Delta t = %.2f \: \mathrm{Myr}$" %(labels[i], dt)
            #tmp = r"%s, $\Delta t = \: ?$" %(labels[i]) 
        '''

        if plot_time == 1:
            ax1.plot(t, 1 - volume_frac[i], color = PlotScripts.colors[i], ls = PlotScripts.linestyles[i], lw = 3, label = labels[i])
        else:
            ax1.plot(ZZ, 1 - volume_frac[i], color = PlotScripts.colors[i], ls = PlotScripts.linestyles[i], lw = 3, label = labels[i])
        #ax3.plot(t, np.subtract(volume_frac[i], volume_frac[0]), color = PlotScripts.colors[i])	
	    
    if plot_time == 1: 
        ax1.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = label_size)
        ax1.xaxis.set_minor_locator(mtick.MultipleLocator(time_tick_interval))
        ax1.set_xlim(time_xlim)

        ax2 = ax1.twiny()

        t_plot = (AllVars.t_BigBang - cosmo.lookback_time(z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
        z_labels = ["$%d$" % x for x in z_plot] # Properly Latex-ize the labels.

        ax2.set_xlabel(r"$z$", size = label_size)
        ax2.set_xlim(time_xlim)
        ax2.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
        ax2.set_xticklabels(z_labels) # But label them as redshifts.

    else:
	
        ax1.set_xlabel(r"$z$", size = label_size)	

	
    leg = ax1.legend(loc='lower left', numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(legend_size + 3)

    ax1.set_ylabel(r"$\langle x_{HI}\rangle$", size = label_size)         
    ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.05))	
    ax1.set_ylim([-0.1, 1.1])
	
	#ax1.set_yscale('log', nonposy='clip')
	#ax1.text(0.075, 0.965, '(c)', horizontalalignment='center', verticalalignment='center', transform = ax1.transAxes)
	
    plt.tight_layout()

    outputFile = OutputDir + output_tag + output_format 			
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))	
			
    plt.close()	

##########

def plot_photo_mean(ZZ, Photo_Mean, Photo_Std, labels, OutputDir, output_tag):

    ax1 = plt.subplot2grid((3,1), (0,0), rowspan =2)
    ax2 = plt.subplot2grid((3,1), (2,0))

    for p in range(0, len(Photo_Mean)):
        ax1.plot(ZZ[1:len(Photo_Mean[p])], np.log10(Photo_Mean[p][1:len(Photo_Mean[p])]), color = PlotScripts.colors[p], label = labels[p])
        ax2.plot(ZZ[1:len(Photo_Mean[p])], np.log10(Photo_Mean[p][1:len(Photo_Mean[p])])/np.log10(Photo_Mean[1][1:len(Photo_Mean[1])]), color = PlotScripts.colors[p])

        leg = ax1.legend(loc='lower left', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts(): 
            t.set_fontsize(legend_size)

    ax2.set_xlabel(r"$\mathrm{z}$", size = label_size + extra_size)
    ax1.set_ylabel(r"$\mathrm{log}_{10} \: \langle \Gamma_\mathrm{HI}\rangle \: \: [\mathrm{s}^{-1}]$", size = label_size + extra_size)
    ax2.set_ylabel(r"$\mathrm{log}_{10} \: \langle \Gamma_\mathrm{HI}\rangle / \mathrm{log}_{10} \: \langle \Gamma_\mathrm{HI}\rangle_{N = 128^3}$", size = label_size - extra_size)

    ax1.get_xaxis().set_visible(False)
    #	plt.axis([6, 12, -16, -9.5])
    ax1.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    ax2.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    ax1.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    ax2.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))

    outputFile = OutputDir + output_tag + output_format 			
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))		
            
    plt.close()	


#########




def plot_total_nion(ZZ, total_nion, simulation_norm, plot_time, labels, OutputDir, output_tag):
# This function will plot the total number of ionizing photons as a function of redshift. The function accepts an array of arrays for the total number of ionizing photons to allow plotting of multiple realizations on the same graph.

## Input ##
# ZZ is an array containing the redshifts we plot over.
# total_nion is an array of arrays where each array contains the total number of ionizing photons at each redshift.
# labels is an array that contains the labels for each realization we are plotting.

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*_*z*.*output_format*'

    print("")
    print("Plotting the total number of ionizing photons") 

    ax1 = plt.subplot2grid((3,1), (0,0), rowspan =3)
    #	ax3 = plt.subplot2grid((3,1), (2,0))

    nion = np.empty((len(total_nion), len(total_nion[0])))

    for model_number in range(0, len(total_nion)):
        if simulation_norm[model_number] == 0:
            cosmo = AllVars.Set_Params_Mysim()
        elif simulation_norm[model_number] == 2:
            cosmo = AllVars.Set_Params_Tiamat()
        elif simulation_norm[model_number] == 6:
            cosmo = AllVars.Set_Params_Kali()
        else:
            raise ValueError("Set simulation norm option in plot_total_nion")
        for snapshot_idx in range(0, len(total_nion[model_number])):
            nion[model_number][snapshot_idx] = np.log10(total_nion[model_number][snapshot_idx] / ((AllVars.BoxSize / AllVars.Hubble_h)**3))	

        print("The sum of ionizing photons for model {0} is {1}".format(model_number, sum(nion[model_number]))) 
        if plot_time == 1:
            ax1.plot((AllVars.t_BigBang- cosmo.lookback_time(ZZ).value) * 1.0e3, nion[model_number], color = PlotScripts.colors[model_number], label = labels[model_number], ls = PlotScripts.linestyles[model_number], lw = 3)
        else:	
            ax1.plot(ZZ, nion[model_number], color = PlotScripts.colors[model_number], label = labels[model_number], ls = PlotScripts.linestyles[model_number], lw = 3)

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


    if plot_time == 1:
        ax1.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr]}$", size = label_size)
        ax1.xaxis.set_minor_locator(mtick.MultipleLocator(time_tick_interval))
        ax1.set_xlim(time_xlim)

        ax1.set_ylabel(r'$\sum f_\mathrm{esc}\dot{N}_\gamma \: [\mathrm{s}^{-1}\mathrm{Mpc}^{-3}]$', fontsize = PlotScripts.global_fontsize)                 
        ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        ax1.set_ylim([48.5, 51.5])

        ax2 = ax1.twiny()

        t_plot = (AllVars.t_BigBang - cosmo.lookback_time(z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
        z_labels = ["$%d$" % x for x in z_plot] # Properly Latex-ize the labels.

        ax2.set_xlabel(r"$z$", fontsize = label_size)
        ax2.set_xlim(time_xlim)
        ax2.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
        ax2.set_xticklabels(z_labels) # But label them as redshifts.
        
    else:
        ax1.set_xlabel(r"$z$", size = label_size)
        ax1.set_ylabel(r"$\mathrm{log}_{10} \: \dot{N}_{\mathrm{HI}} \: [\mathrm{s}^{-1}\mathrm{Mpc}^{-3}]$", fontsize = label_size)
        ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))
        ax1.set_ylim([48.5, 51.5])

        '''
        ax3.set_xlim(time_xlim)
        ax3.set_ylim([-0.3, 0.30])
        ax3.set_yticks(np.arange(-0.20, 0.30, 0.20))
        ax3.yaxis.set_minor_locator(mtick.MultipleLocator(0.05))


        for tick in ax3.yaxis.get_major_ticks():
        tick.label.set_fontsize(label_size - 5)
        '''
    leg = ax1.legend(loc='lower right', numpoints=1,
                     labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(legend_size + 3)

    plt.tight_layout()

    outputFile = OutputDir + output_tag + output_format 			
    plt.savefig(outputFile)  # Save the figure	
    print('Saved file to {0}'.format(outputFile))		
        
    plt.close()	

##########


def plot_density(z, density, OutputDir, output_tag):
# Plots a slice of the density field.

## Input ##
# z is the redshift we're plotting at.
# density is the array that contains the overdensity of each cell.

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*_*z*.*output_format*'

    print("")
    print("Plotting density slice.")
     
    rho_min = np.amin(density)
    rho_max = np.amax(density)

    idx_rho_min = np.argmin(np.ravel(density), axis = 0)
    idx_rho_max = np.argmax(np.ravel(density), axis = 0)

    print("The minimum density for redshift {0:.3f} was {1:.3f}.".format(z, rho_min)) 
    print("The maximum density for redshift {0:.3f} was {1:.3f}.".format(z, rho_max)) 

    print("The minimum density index for redshift {0:.3f} was {1}.".format(z, idx_rho_min))
    print("The maximum density index for redshift {0:.3f} was {1}.".format(z, idx_rho_max)) 

    ax = plt.subplot(111)

    cut_slice = 63
    thickness_cut = 1
    im = ax.imshow(density[:,:,cut_slice:cut_slice+thickness_cut].mean(axis = -1), interpolation='bilinear', origin='low', extent =[0,AllVars.BoxSize,0,AllVars.BoxSize], cmap = 'Purples', vmin = 0.12, vmax = 5) 
    
    cbar = plt.colorbar(im, ax = ax)
    cbar.set_label(r'$\rho/\langle \rho \rangle$')
        
    ax.set_xlabel(r'$\mathrm{y}  (h^{-1}Mpc)$')  
    ax.set_ylabel(r'$\mathrm{x}  (h^{-1}Mpc)$')  

    ax.set_xlim([0.0, AllVars.BoxSize]) 
    ax.set_ylim([0.0, AllVars.BoxSize])

    title = r"$z = %.3f$" %(z)
    ax.set_title(title)
    #ax.text(52.5, 60, r'z = %.1f' %ZZ[i], color = 'black', fontsize = 15)

    outputFile = OutputDir + output_tag + output_format 		
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))			

    plt.close()	

###########

def plot_photofield(z, photofield, OutputDir, output_tag):
# Plots a slice of the HI Photoionization field at a specified redshift.

## Input ##
# z is the redshift we're plotting at.
# photofield is the photoionization rate of each cell. 

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*_*z*.*output_format*'

    photo_min = np.amin(np.log10(photofield))
    photo_max = np.amax(np.log10(photofield))

    print("")
    print("Plotting photoionization slice.")
    print("The minimum photoionization rate for redshift {0:.3f} was {1:.3f}.".format(z, photo_min))
    print("The maximum photoionization rate for redshift {0:.3f} was {1:.3f}.".format(z, photo_max)) 

    print("The minimum photoionization rate index for redshift {0:.3f} was {1:.3f}.".format(z, np.argmin(np.ravel(photofield), axis = 0)))
    print("The maximum photoionization rate index for redshift {0:.3f} was {1:.3f}.".format(z, np.argmax(np.ravel(photofield), axis = 0)))

    photo_slice = photofield[:,:,cut_slice:cut_slice+1].mean(axis = -1)

    ax = plt.subplot(111)

    im = ax.imshow(np.log10(photo_slice), interpolation='bilinear', origin='low', extent =[0,AllVars.BoxSize,0,AllVars.BoxSize], cmap = 'Purples', vmin = -16.22, vmax = -6.8) 

    cbar = plt.colorbar(im, ax = ax)
    cbar.set_label(r'$\mathrm{log}_{10}\Gamma [\mathrm{s}^{-1}]$')
            
    ax.set_xlabel(r'$\mathrm{x}  (h^{-1}Mpc)$')  
    ax.set_ylabel(r'$\mathrm{y}  (h^{-1}Mpc)$')  

    ax.set_xlim([0.0, AllVars.BoxSize]) 
    ax.set_ylim([0.0, AllVars.BoxSize])

    title = r"$z = %.3f$" %(z)
    ax.set_title(title)
    #ax.text(52.5, 60, r'z = %.1f' %ZZ[i], color = 'black', fontsize = 15)

    outputFile = OutputDir + output_tag + output_format 		
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))				

    plt.close()	

############

def plot_nionfield(z, nion, OutputDir, output_tag):
# Plots a slice of the ionizing photon field at a specified redshift.

## Input ##
# z is the redshift we're plotting at.
# photofield is the photoionization rate of each cell. 

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*_*z*.*output_format*'

    if (np.mean(nion) == 0):
        return

    #nion_slice = nion[:,:,cut_slice:cut_slice+1].mean(axis = -1)
    nion_slice = nion[:,:,0:-1].mean(axis = -1)

    nion_min = np.amin(np.log10(nion[nion > 0]))
    nion_max = np.amax(np.log10(nion[nion > 0]))

    print("")
    print("Plotting photoionization slice.")
    print("The minimum ionizing photons rate for redshift {0:.3f} was {1:.3f}.".format(z, nion_min))
    print("The maximum ionizing photons rate for redshift {0:.3f} was {1:.3f}.".format(z, nion_max)) 

    print("The minimum ionizing photons rate index for redshift {0:.3f} was {1:.3f}.".format(z, np.argmin(np.ravel(nion), axis = 0)))
    print("The maximum ionizing photons rate index for redshift {0:.3f} was {1:.3f}.".format(z, np.argmax(np.ravel(nion), axis = 0)))

    ax = plt.subplot(111)

    im = ax.imshow(np.log10(nion_slice), interpolation='bilinear', origin='low', extent =[0,AllVars.BoxSize,0,AllVars.BoxSize], cmap = 'Purples', vmin = nion_min, vmax = nion_max) 


    cbar = plt.colorbar(im, ax = ax)
    cbar.set_label(r'$\mathrm{log}_{10}N_{\gamma} [\mathrm{s}^{-1}]$')
        
    ax.set_xlabel(r'$\mathrm{x}  (h^{-1}Mpc)$')  
    ax.set_ylabel(r'$\mathrm{y}  (h^{-1}Mpc)$')  

    ax.set_xlim([0.0, AllVars.BoxSize]) 
    ax.set_ylim([0.0, AllVars.BoxSize])

    title = r"$z = %.3f$" %(z)
    ax.set_title(title)
    #ax.text(52.5, 60, r'z = %.1f' %ZZ[i], color = 'black', fontsize = 15)

    outputFile = OutputDir + output_tag + output_format 		
    plt.savefig(outputFile)  # Save the figure	
    print('Saved file to {0}'.format(outputFile))				

    plt.close()	

############

def plot_density_numbers(z, density, Ncells, OutputDir, output_tag):
# Plot the numerical value of the overdensity and highlight high values of the overdensities.  Useful if we want to investigate the values around some point.

## Input ##
# z is the redshift we're plotting at.
# density is the array that contains the overdensity of each cell.

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*_*z*.*output_format*'

	# These values will give the bounds of our zoom in. Numerical values are in cMpc.
    lower_x = 64 
    upper_x = 78 
    lower_y = 45 
    upper_y = 58 


    density[density == 0] = 1
    density = np.log10(density)

    step = 1 # The interval we want to plot the values at.
    density_cut = density[:,:,cut_slice:cut_slice+1].mean(axis = -1)
    density_cut = density_cut[lower_x:upper_x:step,lower_y:upper_y:step]
        
    locations_x = np.arange(lower_x,upper_x,step) # Cell locations for each number.
    locations_y = np.arange(lower_y,upper_y,step)

    ax = plt.subplot(111)

    for (p, q), r in np.ndenumerate(density_cut): # Loops over each overdensity.
        if r > 2: # If we have a high overdensity,
            color = 'r' # Want to highlight in a different colour.
        else:
            color = 'k'
        ax.text(locations_y[q]*Ratio, locations_x[p]*Ratio, '{:0.1f}'.format(r), ha='center', va='center', fontsize = 8, color = color) # Write the numerical value.

    ax.set_xlabel(r'$\mathrm{x}  (h^{-1}Mpc)$')  
    ax.set_ylabel(r'$\mathrm{y}  (h^{-1}Mpc)$')  

    ax.set_xlim([(lower_y-1)*BoxSize/Ncells, (upper_y)*BoxSize/Ncells]) 
    ax.set_ylim([(lower_x-1)*BoxSize/Ncells, (upper_x)*BoxSize/Ncells])

    title = r"$z = %.3f$" %(z)

    ax.set_title(title)

    outputFile = OutputDir + output_tag + output_format 
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))				

    plt.close()
##########

def calculate_bubble_MC(z, ionized_cells, Ncell, OutputDir, output_tag):
# Calculate the size of ionized/neutral bubbles by choosing an ionized/neutral cell and measuring the distance to a cell of the opposite phase.

## Input ##
# z is the redshift we are doing the MC walk at.
# ionized_cells is the array that contains the ionization state of a cell.

## Output ##
# The output file will be of the form '*OutputDir*/MC/*output_tag*_*z*.dat'

    def MC_walk(cells, indices, phase, N, Ncell, outfile):
        # Function for the MC walk to calculate bubble size.
    
        ## Input ##
	# cells is the array that contains the ionization field (0 or 1).
	# indices is a length 3 array containing the x,y,z indices of the cells that are neutral (phase = 0) or ionized (phase = 1)
	# phase is 0 if we are starting with neutral or 1 if we are starting with ionized cell.
	# N is the number of random points we select.
    
        ## Output ##
        # A .dat file containing the results of the MC walk (number of steps until it encounters a phase transition).

        radii = []

        for j in range(0,int(N)):
                        
            if (j%20000 == 0):
                print("On walk number {0}".format(j))
    
			## Select a random direction to walk through ##
            direction = random.randint(1,6)

            if direction == 1:
                x = 1
                y = 0
                z = 0
            elif direction == 2:
                x = -1 
                y = 0
                z = 0
            elif direction == 3:
                x = 0
                y = 1
                z = 0
            elif direction == 4:
                x = 0
                y = -1
                z = 0
            elif direction == 5:
                x = 0
                y = 0
                z = 1
            else:
                x = 0
                y = 0
                z = -1

			# Pick the x,y,z coordinates of a random ionized/neutral cell
            random_index = random.randint(0,len(indices[0])-1)
            walk_x = indices[0][random_index]
            walk_y = indices[1][random_index]
            walk_z = indices[2][random_index]

            R = 0
            phase_transition = phase

            while (phase_transition == phase): # While we haven't changed phase yet.
                R += 1 # Keep increasing the radius until we do.
                phase_transition = cells[(walk_x + R*x) % Ncell, (walk_y + R*y) % Ncell, (walk_z + R*z) % Ncell]
                if (phase_transition > 0.8):
                    phase_transition = 1
                else:
                    phase_transition = 0	
                if (R >= Ncell): # If the radius has gone beyond the number of available cells, 
                    phase_transition = (phase + 1) % 2 # We force the change.

            radii.append(R)
			
        np.savetxt(outfile, radii, delimiter = ',')
        print("MC file saved as {0}".format(outfile))

    print("")
    print("Calculating bubble size using MC walk.")	

    N = 1e5

    start_time = time.time()

    ionized_indices= np.where(ionized_cells > 0.8) # Calculate the array indices corresponding to neutral/ionized cells.
    unionized_indices = np.where(ionized_cells < 0.2)

    MCDir = OutputDir + 'MC/' 	
    if not os.path.exists(MCDir):
        os.makedirs(MCDir)
    outfile = MCDir + output_tag + '_z_%.3f.dat' %(z)
    MC_walk(ionized_cells, ionized_indices, 1, N, Ncell, outfile)

    print("MC took {0} seconds.".format(time.time() - start_time))


##########

def plot_bubble_MC(ZZ, fractions_HI, simulation_norm, model_tags, file_tags, Ncell, OutputDir, output_tag):
# Plot the results from the MC walk.  Note that input_tag and labels are arrays that contains the input file tag and label for different realizations.

## Input ##
# ZZ is the redshift range we are plotting the results over.
# upperZ is the upper bound we are plotting to.
# input_tag is an array with the file tags for the input MC file for different realizations #### NOTE: This should be the same as *output_tag* in calculate_bubble_MC
# labels is an array that contains the graph labels for different realizations.

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*.*output_format*'

    print("Plotting results from MC walk.")

    MCDir = OutputDir + 'MC/'

    ax1 = plt.subplot(111)
    for p, q in itertools.product(range(len(fractions_HI)), range(len(file_tags))): # Equivalent to a nested for loop over fractions_HI and file_tags.			
        if simulation_norm[q] == 0:
            cosmo = AllVars.Set_Params_Mysim()
        elif simulation_norm[q] == 2:
            cosmo = AllVars.Set_Params_Tiamat()	
        elif simulation_norm[q] == 5:
            cosmo = AllVars.Set_Param_Britton()
        elif simulation_norm[q] == 6:
            cosmo = AllVars.Set_Params_Kali()
        else:
            print("Have yet to set this simulation parameter.  Should have been caught during initialization thought and no here (in plot_bubble_MC).")
            exit()
        
       

        infile = MCDir + file_tags[q] + '_z_%.3f.dat' %(ZZ[q][p])
      
     
        if (os.path.exists(infile) == False):
            print("Could not find file {0}.  Skipping and moving on".format(infile))
            continue
        fd = open(infile, 'rb')

        print("Plotting Bubble size of file {0}".format(infile))
        R = np.loadtxt(fd)
        print("Maximum radius before scaling is {0} cells.".format(max(R)))
        print("The ratio is: BoxSize = {0:.3f} Mpc/h, Ncell = {1} => 1cell = {2:.3f} Mpc/h".format(AllVars.BoxSize, Ncell[q], AllVars.BoxSize/float(Ncell[q])))		
        R *= AllVars.BoxSize/float(Ncell[q])
        print("Maximum radius after scaling is {0:.4f} Mpc/h.".format(max(R)))

        binwidth = 6*AllVars.BoxSize/float(Ncell[q])

        R_low = binwidth/2 
        if max(R) > AllVars.BoxSize/2:
            R_high = AllVars.BoxSize/2 + binwidth/2 
        else:
            R_high = max(R) + binwidth/2 

        NB = int(np.ceil((R_high - R_low) / binwidth))

        #print(R)
        #print(R_low)
        #print(R_high)
        (counts, binedges) = np.histogram(R, range=(R_low, R_high), bins=NB, density = True)

        # Set the x-axis values to be the centre of the bins
        xaxeshisto = binedges[:-1] + 0.5 * binwidth
        counts *= xaxeshisto
            
        if q == 0:		
            label = r"$\langle x_\mathrm{HI}\rangle = %.2f$" %(fractions_HI[p])
        else:
            label = ""
        ls = PlotScripts.linestyles[q%(len(PlotScripts.linestyles))]
        #ax1.scatter(xaxeshisto, counts, label = label, color = PlotScripts.colors[p%len(PlotScripts.colors)], marker = markers[q%len(markers)])
        ax1.plot(xaxeshisto, counts, label = label, color = PlotScripts.colors[p%len(PlotScripts.colors)], ls = ls) 

    for q in range(0, len(model_tags)):
        ax1.plot(nan, nan, ls = PlotScripts.linestyles[q], label = model_tags[q], color = 'k', lw = 2)

    ax1.set_xscale('log', nonposy='clip')
    ax1.set_xlim([1e-0, 100])
    ax1.set_ylim([0.0, 0.7])

    ax1.set_xlabel(r'$R \: [h^{-1}\mathrm{Mpc}]$', fontsize = label_size)
    ax1.set_ylabel(r'$R \: \: dP/dR$', fontsize = label_size)

    ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.025))
    leg = ax1.legend(loc='upper right', numpoints=1,
         labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(legend_size)

    plt.tight_layout()

    outputFile = OutputDir + output_tag + output_format 
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))				
    plt.close()


##########

def find_max_bubble(MC_ZZ, fractions_HI, MC_mode, model_HI_fractions, model_tags, file_tags, Ncell, OutputDir, output_tag):
#TODO: Update this function header string.


    def acceptance_probability(old_amplitude, new_amplitude, temperature):
        if(new_amplitude > old_amplitude):
            return 1.0
        return exp((new_amplitude - old_amplitude) / temperature)


    print("Finding the maximum and amplitude of the MC bubble distribution.")

    MCDir = OutputDir + 'MC/'

    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)
    ax3 = plt.subplot(223)

    for q in range(0, len(file_tags)):

        peak = []
        amplitude = []


        for p in range(0, len(MC_ZZ)):
            if (MC_mode == 1):
                infile = MCDir + file_tags[q] + '_z_%.3f.dat' %(MC_ZZ[q][p])
            elif (MC_mode == 2): 
                infile = MCDir + file_tags[q] + '_z_%.3f.dat' %(MC_ZZ[p]) 
            if (os.path.exists(infile) == False):
                print("Could not find file {0}.  Skipping and moving on".format(infile))
                continue
            fd = open(infile, 'rb')

            print("Reading in file {0}".format(infile))
            R = np.loadtxt(fd)
            print("Maximum radius before scaling is {0} cells.".format(max(R)))

            R *= BoxSize/Ncell[q] 
            print("Maximum radius after scaling is {1} Mpc/h.".format(max(R)))	

            binwidth = 2*BoxSize/Ncell[q]

            R_low = binwidth/2 
            if max(R) > BoxSize/2:
                R_high = BoxSize/2 + binwidth/2 
            else:
                R_high = max(R) + binwidth/2 		        
       
            NB = np.ceil((R_high - R_low) / binwidth)

            (counts, binedges) = np.histogram(R, range=(R_low, R_high), bins=NB, density = True)

            # Set the x-axis values to be the centre of the bins
            xaxeshisto = binedges[:-1] + 0.5 * binwidth
            counts *= xaxeshisto

            temperature = 1.0
            temperature_min = 0.00001
            cooling_rate = 0.9

            old_peak = NB/2
            old_amplitude = 0.2 

            while(temperature > temperature_min):
                i = 0
                while(i < 100):
                    random_number = np.random.randint(-10, 10)
                    new_peak = old_peak + random_number 
                    new_peak = min(max(new_peak, 0), NB-1) # Truncate the value of the peak to be inside our range.

                    new_amplitude = counts[new_peak]
                    probability = acceptance_probability(old_amplitude, new_amplitude, temperature)

                    #print "old_peak = %d \t new_peak = %d \t old_amplitude = %.4f \t new_amplitde = %.4f \t probability = %.4e \t random_number = %d \t temperature = %.4f" %(old_peak, new_peak, old_amplitude, new_amplitude, probability, random_number, temperature)		

                    if(probability > random.uniform(0,1)):
                        old_peak = new_peak
                        old_amplitude = new_amplitude	
                    i += 1
                temperature *= cooling_rate
                #print "Temperature = %.4f and peak = %.4f Mpc/h" %(temperature, xaxeshisto[new_peak])
					
            print("For fraction {0:.4f} the peak is at {1:.4f} Mpc/h with a probability of {2:.4f}".format(MC_ZZ[p], xaxeshisto[old_peak], old_amplitude))
            peak.append(xaxeshisto[old_peak])
            amplitude.append(old_amplitude)
        if(MC_mode == 1):
            ax1.plot(fractions_HI, peak, label = model_tags[q], color = PlotScripts.colors[q%len(PlotScripts.colors)], ls = PlotScripts.linestyles[q%len(PlotScripts.markers)])
        elif(MC_mode == 2):
            #ax1.plot(fractions_HI[q], peak, label = model_tags[q], color = PlotScripts.colors[q%len(PlotScripts.colors)], ls = PlotScripts.linestyles[q%len(markers)])
            #ax2.plot(fractions_HI[q], amplitude, color = PlotScripts.colors[q%len(PlotScripts.colors)], ls = linestyles[q%len(markers)])
            ax1.plot(MC_ZZ, peak, label = model_tags[q], color = PlotScripts.colors[q%len(PlotScripts.colors)], ls = PlotScripts.linestyles[q%len(PlotScripts.markers)])
            ax2.plot(MC_ZZ, amplitude, color = PlotScripts.colors[q%len(PlotScripts.colors)], ls = PlotScripts.linestyles[q%len(markers)])

		
        ax3.plot(MC_ZZ, model_HI_fractions[q], color = PlotScripts.colors[q%len(PlotScripts.colors)], ls = PlotScripts.linestyles[q%len(PlotScripts.markers)]) 

    #        ax1.set_xscale('log', nonposy='clip')
    #	ax1.set_xlim([0, BoxSize])
    #	ax1.set_ylim([0.0, 0.7])

    ax1.set_yscale('log', nonposy='clip')
    #ax1.set_xscale('log', nonposy='clip')

    #ax1.set_xlabel(r'$x_\mathrm{HI}$', fontsize = label_size)
    ax1.set_xlabel(r"$z$", size = label_size)
    ax1.set_ylabel(r'$\mathrm{Bubble \: Peak} [h^{-1}\mathrm{Mpc}]$', fontsize = label_size)

    #ax2.set_xlabel(r'$x_\mathrm{HI}$', fontsize = label_size)
    ax2.set_xlabel(r"$z$", size = label_size)
    ax2.set_ylabel(r'$\mathrm{Bubble \: Amplitude}$', fontsize = label_size)

    ax2.set_yscale('log', nonposy='clip')
    #ax2.set_xscale('log', nonposy='clip')

    ax3.set_xlabel(r"$z$", size = label_size)
    ax3.set_ylabel(r"$\langle x_{HI}\rangle$", size = label_size)      

    #ax3.set_yscale('log', nonposy='clip')

    plt.tight_layout()

    outputFile = OutputDir + output_tag + output_format 
    plt.savefig(outputFile)  # Save the figure	
    print('Saved file to {0}'.format(outputFile))				
    plt.close()


##########

def Power_Spectrum(ionized_cells, Ncell):
# This function takes the ionized_cells and calculates the associated power spectrum.
# It first determines the dimensions of the bins in k-space.  Then it calculates the 3D Power Spectrum before binning it into 1D.

## Input ##
# ionized_cells is the array that contains the ionization data.
# Ncells is the number of cells along one dimension (Number of cells = Ncell^3).

## Output ##
# k_mid: The middle point for each of the k-space bins.
# PowSpec_1D: The 1-Dimensional Power Spectrum.
# error: The error in PowSpec_1D
    
    ## Calculating the bins in k-space. ##

    mid_point = int(Ncell/2) # Middle of the Grid.
    
    n1 = arange(Ncell)
    n1[1+mid_point:] -= Ncell # Moves to 0, 1, 2, ..., Ncell/2, -Ncell/2, ..., -1  
    n2 = n1**2 
    n_mag = np.sqrt(add.outer(add.outer(n2, n2), n2)).ravel() # Finding the magnitude of k-space vectors.
    
    n_bins = (-1,) + tuple(arange(mid_point-1)+1.5) + (Ncell*2,)
     
    k_bins = digitize(n_mag, n_bins) - 1 # Digitize basically takes an input array (N_Mag) and bins each value according to the binedges defined by N_bins.  The indices of the bin that each element of N_Mag belongs to is returned.
    assert(k_bins.min() == 0)
    assert(k_bins.max() == len(n_bins) -2)

    k_multiplier = 2.0 * np.pi / AllVars.BoxSize # Moving from real-space to k-space.
    
    k_min = (array(n_bins) * k_multiplier)[:-1] # Lower bound for each of the k-space bins.
    k_min[0] = 0
    
    k_max = (array(n_bins) * k_multiplier)[1:] # Upper bound for each of the k-space bins.
    k_max[-1] = mid_point * k_multiplier * sqrt(3.0) 
  
    k_volume = bincount(k_bins) * (k_multiplier)**3 # Volume in k-space (Number of cells of a specified type, divided by the volume of each cell).
    
    ##
    
    transform = ifftn(ionized_cells) # Fourier transform into 3D Power Spectrum.
    PowSpec_3D = square(transform.ravel().real) + square(transform.ravel().imag) # Then take magnitude. 

    v1 = bincount(k_bins, weights=PowSpec_3D) # Takes the 3D Power Spectrum and bins the values into the k-space bins.
    PowSpec_1D = v1 * (1.0/k_volume) # Normalize by the volume of the bins.  Necessary as spheres of larger radii occupy more volume by definition.
    
    # Calculate the error in the 1D Power Spectrum.
    v2 = bincount(k_bins, weights = square(PowSpec_3D)) # Second Order
    v0 = bincount(k_bins) # Zeroth Order
    
    error = sqrt((v2*v0 - v1**2)/(v0-1.0)) / k_volume
    
    k_mid = 0.5 * (k_max + k_min) # Middle point for each of the bins in k-space.

    return k_mid, PowSpec_1D, error

##########

def plot_power(fractions_HI, k_21, Power_21, Error_21, k_XHII, Power_XHII, Error_XHII, fraction_idx_array, model_tags, OutputDir, output_tag, full_debug = 0):
    ## TODO: Update this function header string.

    ax1 = plt.subplot(111)

    if full_debug == 1:
            print("len(k) = {0}".format(len(k_21)))
            print("len(k[0]) = {0}".format(len(k_21[0])))

    for p, q in itertools.product(range(len(k_21)), range(len(k_21[0]))): # Equivalent to a nested for loop over k and k[0]. 
        if p == 0:	
            label = r"$\langle \chi_\mathrm{HI}\rangle = %.2f$" %(fractions_HI[q])
        else:
            label = ""	
    
        if full_debug == 1:
            print("p = {0}".format(p))
            print("len(k[p]) = {0}".format(len(k_21[p])))
            print("q = {0}".format(q))
            print("")
        w = np.where((k_21[p][q] > 9e-2))[0]
        ax1.plot(k_21[p][q][w], np.log10(Power_21[p][q][w]), color = PlotScripts.colors[q], ls = PlotScripts.linestyles[p], label = label, lw = 2, rasterized=True) 
        if (p == 0):
            error = Error_21[p][q][w] * Power_21[p][q][w] 
            ax1.fill_between(k_21[p][q][w], np.log10(Power_21[p][q][w] - error), np.log10(Power_21[p][q][w] + error), color = PlotScripts.colors[q], alpha = 0.3)	

    for p in range(0, len(model_tags)):
        ax1.plot(-1, -5, ls = PlotScripts.linestyles[p], label = model_tags[p], color = 'k', lw = 2)

    ax1.text(0.1, 1.7, r"$\mathrm{Large \: Scales}$", color = 'k', size = PlotScripts.global_fontsize - 2)    
    ax1.text(1.4, 1.7, r"$\mathrm{Small \: Scales}$", color = 'k', size = PlotScripts.global_fontsize - 2)   
    #ax1.text(0.2, 1.0, r"$\mathrm{Error \: Dominated}$", color = 'k', size = PlotScripts.global_fontsize - 2)   
    ax1.arrow(0.2, 1.65, -0.05, 0.0, head_width = 0.01, head_length = 0.005, fc = 'k', ec = 'k') 
    ax1.arrow(1.7, 1.65, 0.7, 0.0, head_width = 0.05, head_length = 0.05, fc = 'k', ec = 'k') 
    #ax1.axvspan(7e-2, 0.5, color = 'k', alpha = 0.5) 
    ax1.set_xlabel(r'$k \: \left[\mathrm{Mpc}^{-1}h\right]$', size = label_size + extra_size)
    #plt.ylabel( r'$P\left(k\right) \: \: \left[\mathrm{Mpc}^3 h^{-3}\right]$', size = label_size + extra_size)
    #ax1.set_ylabel( r'$\bar{\delta T_b^2} \Delta_{21}^2 \left[\mathrm{mK}^2\right]$', size = label_size + extra_size)
    ax1.set_ylabel( r'$\log_{10} \Delta_{21}^2 \left[\mathrm{mK}^2\right]$', size = label_size + extra_size)
    
    ax1.set_xscale('log', nonposy='clip')
    #ax1.set_yscale('log', nonposy='clip')
    leg = ax1.legend(loc='lower right', numpoints=1,
         labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(legend_size)

    ax1.set_xlim([7e-2, 5.5])
    #ax1.set_ylim([np.amin(np.log10(Power[:])), np.amax(np.log10(Power[:]))])	
    ax1.set_ylim([0.9, 1.8])
    #ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    #ax.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))

    plt.tight_layout()

    outputFile = OutputDir + output_tag + output_format 			
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))					

    plt.close()	

    ax1 = plt.subplot(111)
    for p, q in itertools.product(range(len(k_XHII)), range(len(k_XHII[0]))): # Equivalent to a nested for loop over k and k[0]. 
        if p == 0:	
            label = r"$\langle \chi_\mathrm{HI}\rangle = %.2f$" %(fractions_HI[q])
        else:
            label = ""	
    
        if full_debug == 1:
            print("p = {0}".format(p))
            print("len(k[p]) = {0}".format(len(k_XHII[p])))
            print("q = {0}".format(q))
            print("")
        w = np.where((k_XHII[p][q] > 9e-2))[0]
        ax1.plot(k_XHII[p][q][w], np.log10(Power_XHII[p][q][w]), color = PlotScripts.colors[q], ls = PlotScripts.linestyles[p], label = label, lw = 2, rasterized=True) 
   
    
        if (p == 0):
            error = Error_XHII[p][q][w] * Power_XHII[p][q][w]             
            ax1.fill_between(k_XHII[p][q][w], np.log10(Power_XHII[p][q][w] - error), np.log10(Power_XHII[p][q][w] + error), color = PlotScripts.colors[q], alpha = 0.3)	

    for p in range(0, len(model_tags)):
        ax1.plot(-1, -5, ls = PlotScripts.linestyles[p], label = model_tags[p], color = 'k', lw = 2)

    ax1.text(0.1, -1.4, r"$\mathrm{Large \: Scales}$", color = 'k', size = PlotScripts.global_fontsize - 2)    
    ax1.text(1.4, -1.4, r"$\mathrm{Small \: Scales}$", color = 'k', size = PlotScripts.global_fontsize - 2)       
    ax1.arrow(0.2, -1.425, -0.05, 0.0, head_width = 0.01, head_length = 0.005, fc = 'k', ec = 'k') 
    ax1.arrow(1.8, -1.425, 0.7, 0.0, head_width = 0.01, head_length = 0.1, fc = 'k', ec = 'k') 
    ax1.set_xlabel(r'$k \: \left[\mathrm{Mpc}^{-1}h\right]$', size = label_size + extra_size)
    ax1.set_ylabel( r'$\log_{10} \Delta_\mathrm{XHII}^2$', size = label_size + extra_size)
    
    ax1.set_xscale('log', nonposy='clip')
    leg = ax1.legend(loc='lower center', numpoints=1,
         labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(legend_size)

    ax1.set_xlim([7e-2, 5.5])
    ax1.set_ylim([-2.2, -1.3])
    #ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.05))

    plt.tight_layout()

    outputFile = "{0}{1}_XHII{2}".format(OutputDir, output_tag, output_format)
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))					

    plt.close()	


#########


def calculate_reionization_redshift(ionized_cells, density, redshift_array, z, Ncell):

    print("")
    print("Checking to see what cells were reionized at redshift {0:.2f}.".format(z))

    tot_density = 0.0
    tot_density_std = 0.0
    count = 0.0

    for i, j, k in itertools.product(range(Ncell), range(Ncell), range(Ncell)): # Equivalent to nested for loops over x,y, z.
        if (ionized_cells[i,j,k] > 0.9 and redshift_array[i,j,k] == 0):
            tot_density += density[i,j,k]
            count += 1.0

    if (count > 0.0):
        density_z_mean = tot_density/count
        for i, j, k in itertools.product(range(Ncell), range(Ncell), range(Ncell)): # Equivalent to nested for loops over x,y, z.            	
            if (ionized_cells[i,j,k] > 0.9 and redshift_array[i,j,k] == 0):
                redshift_array[i,j,k] = z
                tot_density_std += (density[i,j,k] - density_z_mean)**2 	
		        
        tot_density_std = np.sqrt(tot_density_std / count)
    else:
        density_z_mean = 0.0
        tot_density_std = 0.0

    return redshift_array

##########


def photon_baryon(lowerZ, upperZ, ZZ, total_nion, Hubble_h, OB, Y, labels, OutputDir, output_tag): 


    Baryons = density_baryons * (Mpc_m * 1e2)**3 # Number density of baryons (in cm^-3) to Mpc^-3 
    Baryons *= (BoxSize/Hubble_h)**3 # Total number of baryons in our box

    total_photons = np.zeros((upperZ-lowerZ)*len(total_nion))

    count = 0 

    for p in range(0, len(tota_nion)):
        total_photons[count] = (AllVars.Lookback_Time[0]-AllVars.Lookback_Time[lowerZ]) * total_nion[p][0] * AllVars.Sec_Per_Megayear*1e3
        count += 1
        for i in range(lowerZ+1, upperZ): 
            total_photons[count] = total_photons[count-1] + (AllVars.Lookback_Time[0] - AllVars.Lookback_Time[i])*total_nion[p][count] * AllVars.Sec_Per_Megayear*1e3
            count += 1

    plt.figure()
    for p in range(0, len(total_nion)):
        if (p == 0):
            lower_idx = 0
            upper_idx = len(total_nion[p])
        elif (p == 1):
            lower_idx = len(total_nion[0])
            upper_idx = len(total_nion[0]) + len(total_nion[1])

        plt.scatter(ZZ[0:len(total_nion[p])], np.log10(total_photons[p*len(total_nion[p]):len(total_nion[p])]/Baryons), color = PlotScripts.colors[p], marker = markers[p], label = labels[p])

    plt.xlabel(r"$\mathrm{z}$")
    plt.ylabel(r"$\mathrm{log}_{10} \: \dot{N}_{\mathrm{HI}} \: \: [\mathrm{s}^{-1}]$")

    leg = plt.legend(loc='upper right', numpoints=1,labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(legend_size)

    outputFile = OutputDir + output_tag + output_format 			
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))					    
            
    plt.close()	
    

##########


def analytic_HII(total_nion, redshift, upperZ, snaplist, OutputDir, output_tag, full_debug = 0):

    Q_array = np.empty(upperZ)
    Q_array[0] = 0
    Q = 0

    print("Doing an analytic calculation for ionized hydrogen using the prescription by Lapi 2016.")
    for i in range(0, upperZ-1):

        print("==============================")
        print("REDSHIFT {0:.2f}".format(ZZ[i]))	
        print("==============================")

        ## Using the prescription given by Lapi 2016 ## 

        hydrogen = n_HI(redshift[i], Hubble_h, OB, Y) * 1e-6 * AlVars.cubic_cm_cubic_mpc 
        Nion = total_nion[i] * 1e1 / ((BoxSize/Hubble_h)**3) # Number of ionizing photons in s^-1 Mpc^-3
        one_plus_z = (1 + redshift[i])/float(7)
        t_rec = 3.2 * AllVars.sec_per_Myr * 1e3 * (one_plus_z**(-3)) * (3**(-1)) # Recombination time in seconds. 
        Q_dot = Nion/hydrogen - Q/t_rec
        dt = (AllVars.Lookback_Time[snaplist[i]] - AllVars.Lookback_Time[snaplist[i+1]]) * AllVars.sec_per_Myr * 1e3

        Q = Q + Q_dot*dt
        Q_array[i+1] = Q

        if full_debug == 1:
            print("Hydrogen per Mpc^-3 = {0:.4f}".format(hydrogen))
            print("Nion = {0:.4e} s^-1 Mpc^-3".format(Nion))
            print("one_plus_z {0:.4f}".format(one_plus_z))
            print("t_rec = {0:.4f}".format(t_rec))
            print("Q_dot = {0:.4f}".format(Q_dot))
            print("dt = {0:.4f}".format(dt))
            print("Q = {0:.4f}".format(Q))
    
    plt.figure()

    plt.scatter(redshift[0:upperZ], 1-Q_array, color = 'r', marker = 'x')

    plt.xlabel(r"$\mathrm{z}$")
    plt.ylabel(r"$Q_{HI}$")

    outputFile = OutputDir + output_tag + output_format 			
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))					        
           
    plt.close()	

##########

def plot_redshifts(redshift, ZZ, lowerZ, upperZ, OutputDir, output_tag):

	ax = plt.subplot(111)
	
	im = ax.imshow(redshift[:,:,cut_slice:cut_slice+1].mean(axis = -1), interpolation='bilinear', origin='low', extent =[0,BoxSize,0,BoxSize], cmap = 'Dark2', vmin = ZZ[lowerZ], vmax = ZZ[upperZ-1]) 

	cbar = plt.colorbar(im, ax = ax)
	cbar.set_label(r'$z_{\mathrm{reion}}$')
				    
	ax.set_xlabel(r'$\mathrm{x}  (h^{-1}Mpc)$')  
	ax.set_ylabel(r'$\mathrm{y}  (h^{-1}Mpc)$')  

	ax.set_xlim([0.0, BoxSize]) 
	ax.set_ylim([0.0, BoxSize])

	outputFile = OutputDir + output_tag + output_format 		
	plt.savefig(outputFile)  # Save the figure
	print('Saved file to {0}'.format(outputFile))					        
			
	plt.close()	

##########

def save_redshifts(redshift, OutputDir, output_tag):

	outputFile = OutputDir + output_tag
	redshift.tofile(outputFile)

	print('Saved file to {0}'.format(outputFile))					        	

##########

def plot_density_redshift(ZZ, density_mean, density_std, labels, OutputDir, output_tag): 

    ax = plt.subplot(111)

    for i in range(0, len(density_mean)):
        ax.errorbar(density_mean[i], ZZ, fmt = 'o', xerr = density_std[i], yerr = 0, label = labels[i], color = colours[i])

    ax.set_xlim([0.1, 1.6])
    ax.set_ylim([min(ZZ), max(ZZ)])

    ax.set_xlabel(r'$\Delta$')
    ax.set_ylabel(r'$z_\mathrm{reion}$')

    leg = plt.legend(loc='upper left', numpoints=1,labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(legend_size)


    outputFile = OutputDir + output_tag + output_format 		
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))
            
    plt.close()
	
##########

def plot_density_photons(ZZ, nion, density, count, labels, OutputDir, output_tag): 

    ax = plt.subplot(111)

    for i in range(0, len(nion)):
        w = np.where((nion[i].ravel() > 0.0))
        ax.scatter((density[i].ravel())[w], np.log10(((nion[i].ravel())[w])/((count[i].ravel())[w])), label = labels[i], color = colours[i], alpha = 0.5)

    ax.set_xlim([0.1, 10])
    ax.set_ylim([40, 60])

    ax.set_xlabel(r'$\Delta$')
    ax.set_ylabel(r'$\dot{N}_\Gamma$')

    leg = plt.legend(loc='upper left', numpoints=1,labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(legend_size)

    outputFile = OutputDir + output_tag + output_format 		
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))	
            
    plt.close()

##########

def plot_twentyone(ZZ, ionized_cells, density, OutputDir, output_tag):

    ax = plt.subplot(111)

    T = T_naught(ZZ, OB, Hubble_h, OM)	

    brightness = T * density * ionized_cells
    im = ax.imshow(brightness[:,:,cut_slice:cut_slice+1].mean(axis = -1), interpolation='bilinear', origin='low', extent =[0,BoxSize,0,BoxSize], cmap = 'afmhot_r', norm = PlotScripts.colors.LogNorm(vmin = 1, vmax = 150))

    cbar = plt.colorbar(im, ax = ax)
    cbar.set_label(r'$\delta T_b$')
            
    ax.set_xlabel(r'$\mathrm{x}  (h^{-1}Mpc)$')  
    ax.set_ylabel(r'$\mathrm{y}  (h^{-1}Mpc)$')  

    ax.set_xlim([0.0, BoxSize]) 
    ax.set_ylim([0.0, BoxSize])

    title = r"$z = %.3f$" %(ZZ)
    ax.set_title(title)

    outputFile = OutputDir + output_tag + output_format 		
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))

    plt.close()	


##########

def twentyone_slice(ZZ, z_index, ionized_cells, density, Ncells, brightness_slice):

	T = T_naught(ZZ, OB, Hubble_h, OM)
	
	brightness = T * density * ionized_cells

	for i in range(0, Ncells):
		brightness_slice[z_index, i] = np.mean(brightness[i,i,:])	


	return brightness_slice

##########

def plot_twentyone_slice(ZZ, brightness_slice, Ncells, OutputDir, output_tag):

    ax = plt.subplot(111)
    im = ax.imshow(brightness_slice.T, interpolation='none', extent =[ZZ[-1], ZZ[0],0,BoxSize], cmap = 'winter', origin='low', aspect='auto')

    cbar = plt.colorbar(im, ax = ax)
    cbar.set_label(r'$\delta T_b$')
        
    ax.set_xlabel(r'$z$')
    ax.set_ylabel(r'$\mathrm{x}  (h^{-1}Mpc)$')  

    ax.set_xlim([ZZ[-1], ZZ[0]]) 
    ax.set_ylim([0.0, BoxSize])

    outputFile = OutputDir + output_tag + output_format 		
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))	

    plt.tight_layout()

    plt.close()	
	

##########

def plot_deltat_deltax(ZZ, volume_frac, labels, OutputDir, output_tag):

    print("") 
    print("Plotting the speed of reionization")

    ax = plt.subplot(111)

    t = np.zeros((len(ZZ)-1))
    deltat = np.zeros((len(ZZ)-1))

    deltax_model1 = np.zeros((len(ZZ)-1))
    deltax_model2 = np.zeros((len(ZZ)-1))
    deltax_model3 = np.zeros((len(ZZ)-1))

    plot_frac_model1 = np.zeros((len(ZZ)-1))
    plot_frac_model2 = np.zeros((len(ZZ)-1))
    plot_frac_model3 = np.zeros((len(ZZ)-1))

    for i in range(0, len(ZZ) - 1):
        plot_frac_model1[i] = volume_frac[0][i]
        plot_frac_model2[i] = volume_frac[1][i]	
        plot_frac_model3[i] = volume_frac[1][i]	

        deltat[i] = abs(cosmo.lookback_time(ZZ[i+1]).value  - cosmo.lookback_time(ZZ[i]).value) * 1.0e3

        deltax_model1[i] = volume_frac[0][i+1] - volume_frac[0][i]
        deltax_model2[i] = volume_frac[1][i+1] - volume_frac[1][i]
        deltax_model3[i] = volume_frac[2][i+1] - volume_frac[2][i]

        t[i] = (t_BigBang - cosmo.lookback_time(ZZ[i]).value) * 1.0e3 

    ax.plot(t, deltax_model1/deltat, color = PlotScripts.colors[0], label = labels[0], ls = PlotScripts.linestyles[0])
    ax.plot(t, deltax_model2/deltat, color = PlotScripts.colors[1], label = labels[1], ls = PlotScripts.linestyles[1])
    ax.plot(t, deltax_model3/deltat, color = PlotScripts.colors[2], label = labels[2], ls = PlotScripts.linestyles[2])

    ax.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr}^{-1}]$")
    ax.set_ylabel(r"$\Delta x_\mathrm{HI}/\Delta t \: [\mathrm{Fraction \: Myr}^{-1}]$")

    ax.set_xlim(time_xlim)
    ax.set_ylim([-0.009, 0.001]) 

    leg = plt.legend(loc='lower left', numpoints=1,labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for o in leg.get_texts():  # Reduce the size of the text
        o.set_fontsize(legend_size)

    ax.xaxis.set_minor_locator(mtick.MultipleLocator(25))
    ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.0005))

    ax.text(0.075, 0.965, '(d)', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)	
    ## Creating a second x-axis on the top to track redshift. ##

    ax2 = ax.twiny()

    t_plot = (t_BigBang - cosmo.lookback_time(z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
    z_labels = ["$%d$" % x for x in z_plot] # Properly Latex-ize the labels.

    ax2.set_xlabel(r"$z$")
    ax2.set_xlim(time_xlim)
    ax2.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
    ax2.set_xticklabels(z_labels) # But label them as redshifts.

    ##

    outputFile = OutputDir + output_tag + output_format 			
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))	

    plt.close()	

##########

def plot_deltat_deltaN(ZZ, Nion, labels, OutputDir, output_tag):

    print("") 
    print("Plotting the speed of ionizing photons (with time).")

    ax = plt.subplot(111)

    deltat = np.zeros((len(ZZ)-1))

    deltaN_model1 = np.zeros((len(ZZ)-1))
    deltaN_model2 = np.zeros((len(ZZ)-1))
    deltaN_model3 = np.zeros((len(ZZ)-1))

    t = np.zeros((len(ZZ)-1))

    for i in range(0, len(ZZ) - 1):	
        deltaN_model1[i] = np.log10(Nion[0][i+1]/Nion[0][i]/(AllVars.BoxSize / AllVars.Hubble_h)**3)
        deltaN_model2[i] = np.log10(Nion[1][i+1]/Nion[1][i]/(AllVars.BoxSize / AllVars.Hubble_h)**3)
        deltaN_model3[i] = np.log10(Nion[2][i+1]/Nion[2][i]/(AllVars.BoxSize / AllVars.Hubble_h)**3)

        deltat[i] = (cosmo.lookback_time(ZZ[i]).value  - cosmo.lookback_time(ZZ[i+1]).value) * 1.0e3
        t[i] = (t_BigBang - cosmo.lookback_time(ZZ[i]).value) * 1.0e3 

    ax.plot(t, deltaN_model1/deltat, color = PlotScripts.colors[0], label = labels[0], ls = PlotScripts.linestyles[0])
    ax.plot(t, deltaN_model2/deltat, color = PlotScripts.colors[1], label = labels[1], ls = PlotScripts.linestyles[1])
    ax.plot(t, deltaN_model3/deltat, color = PlotScripts.colors[2], label = labels[2], ls = PlotScripts.linestyles[2])

    ax.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr}^{-1}]$") 
    ax.set_ylabel(r"$\Delta \left(\mathrm{log}_{10} \: \dot{N}_\mathrm{HI}\right)/\Delta t \: \mathrm{[Photons \: Myr}^{-1}\mathrm{Mpc}^{-3}]$")

    ax.xaxis.set_minor_locator(mtick.MultipleLocator(time_tick_interval))
    ax.yaxis.set_minor_locator(mtick.MultipleLocator(0.0005))
    ax.set_xlim(time_xlim)

    leg = plt.legend(loc='upper right', numpoints=1,labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(legend_size)

    ax.text(0.075, 0.965, '(b)', horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)		

    ## Creating a second x-axis on the top to track redshift. ##
    ax2 = ax.twiny()

    t_plot = (t_BigBang - cosmo.lookback_time(z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
    z_labels = ["$%d$" % x for x in z_plot] # Properly Latex-ize the labels.

    ax2.set_xlabel(r"$z$")
    ax2.set_xlim(time_xlim)
    ax2.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
    ax2.set_xticklabels(z_labels) # But label them as redshifts.

    ##

    outputFile = OutputDir + output_tag + output_format 			
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))

    LinearFit = np.polyfit(ZZ, np.log10(Nion[0]), 1)
    print("The linear fit y = a*z + b had coefficients a = {0:.4e} and b = {1:.4e} for {2} model".format(LinearFit[0], LinearFit[1], labels[0]))

    LinearFit = np.polyfit(ZZ, np.log10(Nion[1]), 1)
    print("The linear fit y = a*z + b had coefficients a = {0:.4e} and b = {1:.4e} for {2} model".format(LinearFit[0], LinearFit[1], labels[1]))
            
    LinearFit = np.polyfit(ZZ, np.log10(Nion[2]), 1)
    print("The linear fit y = a*z + b had coefficients a = {0:.4e} and b = {1:.4e} for {2} model".format(LinearFit[0], LinearFit[1], labels[2]))
        
    plt.close()
	
##########


def plot_combined_global_nion(ZZ, total_nion, volume_frac, labels, OutputDir, output_tag):

    print("") 
    print("A four panel plot containg the global fraction and Nion data.")

    ax1 = plt.subplot(221)

    for p in range(0, len(total_nion)):
        nion = np.empty(len(total_nion[p]))
        for i in range(0, len(total_nion[p])):
            nion[i] = np.log10(total_nion[p][i] / ((AllVars.BoxSize / AllVars.Hubble_h)**3))	
	
        ax1.plot((t_BigBang- cosmo.lookback_time(ZZ[0:len(total_nion[p])]).value) * 1.0e3, nion, color = PlotScripts.colors[p], label = labels[p], ls = PlotScripts.linestyles[p])

        ax1.yaxis.set_minor_locator(mtick.MultipleLocator(0.1))

    ## Creating a second x-axis on the top to track redshift. ##
    ax1_top = ax1.twiny()

    t_plot = (t_BigBang - cosmo.lookback_time(z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
    z_labels = ["$%d$" % x for x in z_plot] # Properly Latex-ize the labels.

    ax1_top.set_xlabel(r"$z$")
    ax1_top.set_xlim(time_xlim)
    ax1_top.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
    ax1_top.set_xticklabels(z_labels) # But label them as redshifts.

    ##

    ax2 = plt.subplot(222)

    t = np.zeros((len(ZZ)-1))
    deltat = np.zeros((len(ZZ)-1))

    deltax_model1 = np.zeros((len(ZZ)-1))
    deltax_model2 = np.zeros((len(ZZ)-1))
    deltax_model3 = np.zeros((len(ZZ)-1))

    for i in range(0, len(ZZ) - 1):
        deltat[i] = abs(cosmo.lookback_time(ZZ[i+1]).value  - cosmo.lookback_time(ZZ[i]).value) * 1.0e3

        deltax_model1[i] = volume_frac[0][i+1] - volume_frac[0][i]
        deltax_model2[i] = volume_frac[1][i+1] - volume_frac[1][i]
        deltax_model3[i] = volume_frac[2][i+1] - volume_frac[2][i]

        t[i] = (t_BigBang - cosmo.lookback_time(ZZ[i]).value) * 1.0e3 

    ax2.plot(t, deltax_model1/deltat, color = PlotScripts.colors[0], label = labels[0], ls = PlotScripts.linestyles[0])
    ax2.plot(t, deltax_model2/deltat, color = PlotScripts.colors[1], label = labels[1], ls = PlotScripts.linestyles[1])
    ax2.plot(t, deltax_model3/deltat, color = PlotScripts.colors[2], label = labels[2], ls = PlotScripts.linestyles[2])


    ax2.set_xlim(time_xlim)
    ax2.set_ylim([-0.009, 0.001]) 

    ax2.yaxis.set_minor_locator(mtick.MultipleLocator(0.001))

    ## Creating a second x-axis on the top to track redshift. ##

    ax2_top = ax2.twiny()

    t_plot = (t_BigBang - cosmo.lookback_time(z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
    z_labels = ["$%d$" % x for x in z_plot] # Properly Latex-ize the labels.

    ax2_top.set_xlabel(r"$z$")
    ax2_top.set_xlim(time_xlim)
    ax2_top.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
    ax2_top.set_xticklabels(z_labels) # But label them as redshifts.

    ##

    ##### Bottom Left: Global Fraction Function of Time #####

    ax3 = plt.subplot(223)	

    t = np.empty(len(volume_frac[0]))

    for i in range(0, len(volume_frac[0])):
        t[i] = (t_BigBang - cosmo.lookback_time(ZZ[i]).value) * 1.0e3 	

    for i in range(0, len(volume_frac)):
        dt = (cosmo.lookback_time(MC_ZZ[i][0]).value - cosmo.lookback_time(MC_ZZ[i][-1]).value) * 1.0e3
        tmp = r"%s, $\Delta t = %.2f \: \mathrm{Myr}^{-1}$" %(labels[i], dt)
        ax3.plot(t, volume_frac[i], color = PlotScripts.colors[i], label = tmp, ls = PlotScripts.linestyles[i])
        
    ax3.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr}^{-1}\mathrm{]}$")
    ax3.set_ylabel(r"$\langle x_{HI}\rangle$")

    ax3.xaxis.set_minor_locator(mtick.MultipleLocator(time_tick_interval))
    ax3.yaxis.set_minor_locator(mtick.MultipleLocator(0.05))

    ax3.set_ylim([-0.1, 1.1])
    ax3.set_xlim(time_xlim)	

    ## Creating a second x-axis on the top to track redshift. ##
    ax3_top = ax3.twiny()

    t_plot = (t_BigBang - cosmo.lookback_time(z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
    z_labels = ["$%d$" % x for x in z_plot] # Properly Latex-ize the labels.

    ax3_top.set_xlabel(r"$z$")
    ax3_top.set_xlim(time_xlim)
    ax3_top.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
    ax3_top.set_xticklabels(z_labels) # But label them as redshifts.

    ##

    ax4 = plt.subplot(224)

    t = np.zeros((len(ZZ)-1))
    deltat = np.zeros((len(ZZ)-1))

    deltax_model1 = np.zeros((len(ZZ)-1))
    deltax_model2 = np.zeros((len(ZZ)-1))
    deltax_model3 = np.zeros((len(ZZ)-1))

    plot_frac_model1 = np.zeros((len(ZZ)-1))
    plot_frac_model2 = np.zeros((len(ZZ)-1))
    plot_frac_model3 = np.zeros((len(ZZ)-1))

    for i in range(0, len(ZZ) - 1):
        plot_frac_model1[i] = volume_frac[0][i]
        plot_frac_model2[i] = volume_frac[1][i]	
        plot_frac_model3[i] = volume_frac[1][i]	

        deltat[i] = abs(cosmo.lookback_time(ZZ[i+1]).value  - cosmo.lookback_time(ZZ[i]).value) * 1.0e3

        deltax_model1[i] = volume_frac[0][i+1] - volume_frac[0][i]
        deltax_model2[i] = volume_frac[1][i+1] - volume_frac[1][i]
        deltax_model3[i] = volume_frac[2][i+1] - volume_frac[2][i]

        t[i] = (t_BigBang - cosmo.lookback_time(ZZ[i]).value) * 1.0e3 

    ax4.plot(t, deltax_model1/deltat, color = PlotScripts.colors[0], label = labels[0], ls = PlotScripts.linestyles[0])
    ax4.plot(t, deltax_model2/deltat, color = PlotScripts.colors[1], label = labels[1], ls = PlotScripts.linestyles[1])
    ax4.plot(t, deltax_model3/deltat, color = PlotScripts.colors[2], label = labels[2], ls = PlotScripts.linestyles[2])

    ax4.set_xlabel(r"$\mathrm{Time \: Since \: Big \: Bang \: [Myr}^{-1}]$")
    ax4.set_ylabel(r"$\Delta x_\mathrm{HI}/\Delta t \: [\mathrm{Fraction \: Myr}^{-1}]$")

    ax4.set_xlim(time_xlim)
    ax4.set_ylim([-0.009, 0.001]) 

    ax4.xaxis.set_minor_locator(mtick.MultipleLocator(25))
    ax4.yaxis.set_minor_locator(mtick.MultipleLocator(0.001))

    ## Creating a second x-axis on the top to track redshift. ##
    ax4_top = ax4.twiny()

    t_plot = (t_BigBang - cosmo.lookback_time(z_plot).value) * 1.0e3 # Corresponding Time values on the bottom.
    z_labels = ["$%d$" % x for x in z_plot] # Properly Latex-ize the labels.

    ax4_top.set_xlabel(r"$z$")
    ax4_top.set_xlim(time_xlim)
    ax4_top.set_xticks(t_plot) # Set the ticks according to the time values on the bottom,
    ax4_top.set_xticklabels(z_labels) # But label them as redshifts.

    ##

    plt.tight_layout()
    outputFile = OutputDir + output_tag + output_format 			
    plt.savefig(outputFile)  # Save the figure	
    print('Saved file to {0}'.format(outputFile))


##########


def plot_nine_panel_slices(ZZ, filepaths, GridSizes, simulation_norm, MC_Snaps, fractions_HI, model_tags, OutputDir, output_tag):	
    ## TODO: Update the function header string.

    fig = plt.figure(figsize = (16, 12))
    gs1 = gridspec.GridSpec(3,len(filepaths))
    gs1.update(wspace=0.00, hspace=0.0) # set the spacing between axes. 
    for i in range(1, 3*len(filepaths) + 1):
        ax = plt.subplot(gs1[i-1]) 

        model_index = int((i-1) % len(filepaths)) # Every 3rd step it flips back to the first model.
        count_index = int(floor((i-1)/len(filepaths))) # Every 3rd step it increments to the next fraction we're plotting.	
        snapshot_index = int(MC_Snaps[model_index, count_index]) # Corresponding snapshot that we are loading for this model and HI fraction.
              
        if(simulation_norm[model_index] == 0):
            cosmo = AllVars.Set_Params_Mysim()
        elif(simulation_norm[model_index] == 2):
            cosmo = AllVars.Set_Params_Tiamat()
        elif(simulation_norm[model_index] == 5):
            cosmo = AllVars.Set_Params_Britton()
        elif(simulation_norm[model_index] == 6):
            cosmo = AllVars.Set_Params_Kali()
        else:
            print("Normalization should be < 6 and should have been caught during initialization instead of here (plot_nine_panel_slices).")
            exit()
	
        if (snapshot_index < 10):
            number_tag_anne = '_0%d' %(snapshot_index)
            number_tag_mine = '_00%d' %(snapshot_index)
        else:
            number_tag_anne = '_%d' %(snapshot_index)
            number_tag_mine = '_0%d' %(snapshot_index)

        if snapshot_index < 0 or snapshot_index > len(ZZ):	
            ionized_cells = np.zeros((GridSizes[model_index], GridSizes[model_index], GridSizes[model_index]))
            redshift = 0.00
        elif ZZ[snapshot_index] < 5 or ZZ[snapshot_index] > 20:
            ionized_cells = np.zeros((GridSizes[model_index], GridSizes[model_index], GridSizes[model_index]))
            redshift = 0.00
        else:
            fname = filepaths[model_index] + "_%02d" %(snapshot_index) 
            print("")
            print("Loading in data for {0} Model from file {1}".format(model_tags[model_index], fname))
            fd = open(fname, 'rb')
            ionized_cells = np.fromfile(fd, count = GridSizes[model_index]*GridSizes[model_index]*GridSizes[model_index], dtype = np.float64)	
            ionized_cells.shape = (GridSizes[model_index], GridSizes[model_index], GridSizes[model_index])
            fd.close()

            if len(ionized_cells[ionized_cells == 1]):
                raise ValueError("Arghhh")
			
            ionized_cells = np.log10(1 - ionized_cells)

            #ionized_cells = np.log10(ionized_cells)
            redshift = ZZ[snapshot_index]

        index_cut = int(cut_slice * (AllVars.BoxSize/100.0)*(GridSizes[model_index]/GridSizes[0])) # Wish to cut the box at the same spatial point for all models.  So normalize the index that this corresponds to to model1.
        thickness_cut = int(np.ceil(1 * (AllVars.BoxSize/100.0)*(GridSizes[model_index]/GridSizes[0]))) # We will take a thickness of 1 cell for model1 and then normalize the number of cells we need to get the same thickness for the other models.

        print("For model {0} we begin our cut at index {1} (corresponding to physical position of {2:.4f}) and cut with a thickness of {3} cells.".format(model_index, index_cut, index_cut * AllVars.BoxSize/float(GridSizes[model_index]), thickness_cut))

        im = ax.imshow(ionized_cells[:,:,index_cut:index_cut+thickness_cut].mean(axis = -1), interpolation='none', origin='low', extent =[0,AllVars.BoxSize,0,AllVars.BoxSize], vmin = -8, vmax = 0, cmap = 'afmhot_r')
            
        ax.set_xlim([0.0, AllVars.BoxSize]) 
        ax.set_ylim([0.0, AllVars.BoxSize])

        if (i <= len(filepaths)): 
            title = r"%s" %(model_tags[i-1])
            ax.set_title(title, fontsize = PlotScripts.global_labelsize) 

        if (model_index == 0): 
            tmp = (r"$\langle \chi_\mathrm{HI}\rangle = %.2f$") %(1 - fractions_HI[count_index])
            ax.text(-0.65,0.5, tmp, transform = ax.transAxes, fontsize = PlotScripts.global_labelsize) 

        tmp = r"$z = %.2f$" %(redshift)
        txt = ax.text(0.6,0.9, tmp, transform = ax.transAxes, fontsize = PlotScripts.global_labelsize, color = 'k')
        txt.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])
        plt.draw()
        frame1 = plt.gca()
        frame1.axes.get_xaxis().set_visible(False)
        frame1.axes.get_yaxis().set_visible(False)
        ax.set_aspect('equal')
	
    cax = fig.add_axes([0.90, 0.1, 0.03, 0.8])
    cbar = fig.colorbar(im, cax=cax, ticks = np.arange(-8.0, 1.0, 1.0))
    cbar.ax.set_ylabel(r'$\mathrm{log}_{10}\left(\chi_\mathrm{HI}\right)$', rotation = 90, size = PlotScripts.global_labelsize)
    cbar.ax.tick_params(labelsize = PlotScripts.global_legendsize + 10)
    fig.subplots_adjust(right = None, hspace = 0.0, wspace = 0.0)

    outputFile = OutputDir + output_tag + output_format 
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))	
            
    plt.close()

##########

def plot_nine_panel_photHI(ZZ, filepaths, GridSizes, MC_Snaps, fractions_HI, model_tags, OutputDir, output_tag):	
    #TODO: Update function header string.
    
    fig = plt.figure(figsize = (16,12))
    gs1 = gridspec.GridSpec(3,3)
    gs1.update(wspace=0.00, hspace=0.0) # set the spacing between axes. 
    for i in range(1, 9 + 1):
        ax = plt.subplot(gs1[i-1]) 

        model_index = int((i-1) % 3) # Every 3rd step it flips back to the first model.
        count_index = int(floor((i-1)/3)) # Every 3rd step it increments to the next fraction we're plotting.	
        snapshot_index = MC_Snaps[model_index, count_index] # Corresponding snapshot that we are loading for this model and HI fraction.

        if (snapshot_index < 10):
            number_tag_anne = '_0%d' %(snapshot_index)
            number_tag_mine = '_00%d' %(snapshot_index)
        else:
            number_tag_anne = '_%d' %(snapshot_index)
            number_tag_mine = '_0%d' %(snapshot_index)

	
        if snapshot_index < 0 or snapshot_index > 150:
            ionized_cells = np.zeros((GridSizes[model_index], GridSizes[model_index], GridSizes[model_index]))
        elif ZZ[snapshot_index] < 5 or ZZ[snapshot_index] > 20 :	
            ionized_cells = np.zeros((GridSizes[model_index], GridSizes[model_index], GridSizes[model_index]))
        else:
            fname = filepaths[model_index] + number_tag_anne 
            print("")
            print("Loading in data for {0} Model from file {1}".format(model_tags[model_index], fname))
            fd = open(fname, 'rb')
            ionized_cells = np.fromfile(fd, count = GridSizes[model_index]*GridSizes[model_index]*GridSizes[model_index], dtype = np.float64)	
            ionized_cells.shape = (GridSizes[model_index], GridSizes[model_index], GridSizes[model_index])
            fd.close()

        index_cut = int(cut_slice * (GridSizes[model_index]/GridSizes[0])) # Wish to cut the box at the same spatial point for all models.  So normalize the index that this corresponds to to model1.
        thickness_cut = int(1 * (GridSizes[model_index]/GridSizes[0])) # We will take a thickness of 1 cell for model1 and then normalize the number of cells we need to get the same thickness for the other models.
        print("For model {0} we begin our cut at index {1} (corresponding to physical position of {2:.4f}) and cut with a thickness of {3} cells.".format(model_index, index_cut, index_cut * BoxSize/float(GridSizes[model_index]), thickness_cut))

        photo_cut = ionized_cells[:,:,index_cut:index_cut+thickness_cut].mean(axis = -1)
        print("For model {0} At Redshift {1:.2f} the log-mean photoionization rate is {2:.4f}".format(model_tags[model_index], ZZ[snapshot_index], np.log10(np.mean(photo_cut))))
        im = ax.imshow(np.log10(photo_cut), interpolation='none', origin='low', extent =[0,BoxSize,0,BoxSize], vmin = -15, vmax = -8, cmap = 'Purples')
            
        ax.set_xlim([0.0, BoxSize]) 
        ax.set_ylim([0.0, BoxSize])

        if (i == 1 or i == 2 or i == 3):
            title = r"%s" %(model_tags[i-1])
            ax.set_title(title, fontsize = legend_size + 6)
        if (i == 1 or i == 4 or i == 7):
            tmp = (r"$x_\mathrm{HI} = %.3f$") %(fractions_HI[count_index])
            ax.text(10.0,0.5, tmp, transform = ax.transAxes, fontsize = legend_size + 6)

        frame1 = plt.gca()
        frame1.axes.get_xaxis().set_visible(False)
        frame1.axes.get_yaxis().set_visible(False)
        ax.set_aspect('equal')

    cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    cbar = fig.colorbar(im, cax=cax, ticks = np.arange(-15, -7))
    cbar.ax.set_ylabel(r"$\log_{10} \: \Gamma \: [\mathrm{s}^{-1}]$",  rotation = 90, fontsize = legend_size)		
    cbar.ax.tick_params(labelsize = legend_size - 2)

    outputFile = OutputDir + output_tag + output_format 
    plt.savefig(outputFile)  # Save the figure	
    print('Saved file to {0}'.format(outputFile))	
            
    plt.close()

##########


### TODO: Use mass-averaged ionization fraction instead of volume-averaged.
def plot_optical_depth(ZZ, volume_frac, model_tags, OutputDir, output_tag):

    print("Plotting the optical depth")	
    def integrand(z):
        H = Hubble_Param(z, AllVars.Hubble_h, AllVars.Omega_m) / (AllVars.pc_to_m * 1.0e6 / 1.0e3) # Hubble Parameter in Mpc / s / Mpc. 
        return (((1 + z)**2) / H)
 
    tau_04 = integrate.quad(integrand, 0, 4)[0] # Tau for z = 0 to 4.
    tau_04 *= (1 + 2*AllVars.Y/(4 * (1-AllVars.Y)))

    tau_46 = integrate.quad(integrand, 4, ZZ[-1])[0] # Tau for z = 4 to 6.
    tau_46 *= (1 + AllVars.Y/(4* (1-AllVars.Y)))

    tau_06 = tau_04 + tau_46

    ax = plt.subplot(111)
    ax.fill_between(np.arange(0, 20, 0.01), 0.058 - 0.012, 0.058 + 0.012, color = 'k', alpha = 0.3) 
    ax.axhline(y = 0.058, xmin = 0, xmax = 20, color = 'k', alpha = 0.3)

    ax.text(12, 0.055, r"$\mathrm{Planck \: 2016}$")

    for p in range(0, len(volume_frac)):
        tau = np.empty((len(ZZ)))

        for i in range(len(ZZ)-1, -1, -1):
            if p < len(volume_frac):
                XHII = volume_frac[p][i]
            elif ZZ[i] < 8:
                XHII = 0
            else:
                XHII = 1
	
            H = Hubble_Param(ZZ[i], AllVars.Hubble_h, AllVars.Omega_m) / (AllVars.pc_to_m * 1.0e6 / 1.0e3) # Hubble Parameter in Mpc / s / Mpc.	
            numerator = ((1 + ZZ[i]) **2) * XHII 

            if (i == len(ZZ) - 1):
                tau[i] = tau_06 + (( numerator / H) * (ZZ[i] - ZZ[-1]) * (1 + AllVars.Y/(4 * (1 - AllVars.Y)))) 
            else:
                tau[i] = tau[i+1] + (( numerator / H) * (ZZ[i] - ZZ[i+1]) * (1 + AllVars.Y/(4 * (1-AllVars.Y)))) 
			
            # Some debug print statements on Github. #
           
            ''' 
            if (i == len(ZZ) - 1):
                print("z = {0:.2f} \t 1-Volume_frac = {1:.4e} \t numerator = {2:.4e}\t tau[{3}] = {4:.4e} \t tau[{5}] = {6:.4e} \t nHI = {7:.4e} \t H = {8:.4e} \t numerator / H = {10:4e}".format(ZZ[i], 1 - volume_frac[0][i], numerator, i, tau[i], i, tau[i], n_HI(ZZ[i], AllVars.Hubble_h, AllVars.Omega_b, AllVars.Y), H, numerator / H))
            else: 
                print("z = {0:.2f} \t 1-Volume_frac = {1:.4e} \t numerator = {2:.4e}\t tau[{3}] = {4:.4e} \t tau[{5}] = {6:.4e} \t nHI = {7:.4e} \t H = {8:.4e} \t numerator / H = {10:4e}".format(ZZ[i], 1 - volume_frac[0][i], numerator, i+1, tau[i+1], i, tau[i], n_HI(ZZ[i], AllVars.Hubble_h, AllVars.Omega_b, AllVars.Y), H, numerator / H)) 
            ''' 


        if p < len(volume_frac):
            label = model_tags[p]
        else:
            label = label_planck
        tau *= n_HI(0, AllVars.Hubble_h, AllVars.Omega_b, AllVars.Y) * AllVars.c_in_ms * AllVars.Sigmat
        print("tau {0}".format(tau))
        ax.plot(ZZ, tau, label = label, color = PlotScripts.colors[p], ls = PlotScripts.linestyles[p], lw = paper_lineweight)	

    ax.set_xlabel(r"$z$")
    ax.set_ylabel(r"$\tau$")

    ax.set_xlim([6, 13])
    ax.set_ylim([0.04, 0.08])

    leg = ax.legend(loc='upper left', numpoints=1,labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(legend_size)

    outputFile = OutputDir + output_tag + output_format 
    plt.savefig(outputFile)	
    print('Saved file to {0}'.format(outputFile))	
            
    plt.close()


##

def hoshen_kopelman(ionized_cells):

	## Just a quick function that replicates the !! behaviour in C.
	## If the input value (a) is != 0 it returns 1 otherwise it returns 0.


	print("Running the Hosen-Kopelman Algorithm")

	def double_not(a):  
		if(a != 0):
			return 1
		else:
			return 0

	l = len(ionized_cells)
	m = len(ionized_cells)
	n = len(ionized_cells)
	

	max_labels = l*m*n / 2 # The maximum number of discrete ionized regions.
	labels = np.zeros((max_labels), dtype = np.int)

	test = np.zeros((l, m, n), dtype = np.int32)
	for i in range(0, l):
		for j in range(0, m):
			for k in range(0, n):
				if (ionized_cells[i][j][k] > 0.8):
					test[i][j][k] = 1
				else:
					test[i][j][k] = 0
	
	def make_set(label_number):
		label_number += 1	
		assert(label_number < max_labels), "The current_label value is greater than the maximum label limit."
		labels[label_number] = label_number 
		return label_number 

	def find(x):
		y = x
		while (labels[y] != y):
			y = labels[y]

		while(labels[x] != x):
			z = labels[x]
			labels[x] = y
			x = z

		return y

	def union(x, y):
		labels[find(x)] = find(y)
		return find(y)

	def union_3D(x, y, z):
		labels[find(x)] = find(y)
		labels[find(z)] = find(y)
		return find(y) 

	current_label = 0

	for i in range(0, l):
		for j in range(0, m):
			for k in range(0, l):

				if(test[i][j][k] == 1):

					if(i == 0):
						up = 0
					else:
						up = test[i-1][j][k]

					if(j == 0):
						left = 0
					else:
						left = test[i][j-1][k]

					if(k == 0):
						back = 0
					else:
						back = test[i][j][k-1]

					tmp = double_not(left) + double_not(up) + double_not(back) # Since the labels can be greater than 1, double_not returns 1 if the value is >= otherwise it returns 0.
					# So it will return 1 if there is a cluster at the position.

					if(tmp == 0): # We have a new cluster
						test[i][j][k] = make_set(current_label)
						current_label += 1
	
					elif(tmp == 1): # Part of an existing cluster
						test[i][j][k] = max(up, left, back)

					elif(tmp == 2): # Joins two existing clusters together
						if(up == 0):
							test[i][j][k] = union(left, back)
						elif(left == 0):
							test[i][j][k] = union(up, back)
						elif(back == 0):
							test[i][j][k] = union(up, left)
					elif(tmp == 3): # Joins three existing clusters together
						test[i][j][k] = union_3D(up, left, back)

	new_labels = np.zeros((max_labels), dtype = np.int)
	size_bubbles = np.zeros((max_labels), dtype = np.int32)
	new_labels_len = 0

	print("Finished labelling all the cells. Starting the relabelling.")

	for i in range(0, l):
		for j in range(0, m):
			for k in range(0, n):
				if(test[i][j][k] > 0):
					x = find(test[i][j][k])
					if(new_labels[x] == 0):
						new_labels_len += 1
						new_labels[x] = new_labels_len

					test[i][j][k] = new_labels[x]
					size_bubbles[test[i][j][k]] += 1

	total_clusters = new_labels_len		

	print("There was {0} clusters found".format(total_clusters))

	largest_bubble = max(size_bubbles)
	ionization_fraction = 1 - calculate_volume_frac(ionized_cells, 128)
	num_cells_ionized = 128**3 * ionization_fraction 

	print("The largest bubble contains {0} cells. This is {1:.4f} of the entire box.  Ionization Fraction = {2:.4f}. The number of cells ionized is {3}. Num Cells Inside Bubble/Num Cells Ionized = {4:.4f}".format(largest_bubble, float(largest_bubble)/float(l*m*n), ionization_fraction, num_cells_ionized, float(largest_bubble)/float(num_cells_ionized))) 

	return float(largest_bubble)/float(num_cells_ionized)

##


def plot_hoshen(order_parameter, HI_frac_model, model_tags, OutputDir, output_tag):

    ax1 = plt.subplot(111)

    print("Plotting the Hoshen-Koppelman results.")

    for p in range(0, len(HI_frac_model)):
        ax1.plot(np.subtract(1, HI_frac_model[p]), order_parameter[p], lw = 3, color = PlotScripts.colors[p], label = model_tags)

    ax1.set_xlabel(r"$x_\mathrm{HII}$", size = label_size)
    ax1.set_ylabel(r"$P_\mathrm{inf}/x_\mathrm{HI}$", size = label_size)      

    leg = ax1.legend(loc='upper right', numpoints=1,
             labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(legend_size)

    plt.tight_layout()

    outputFile = OutputDir + output_tag + output_format 
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))		
    plt.close()

##

################################################
########## Functions for calculations ##########
################################################

def modes_to_pspec(modes, boxsize):
    """
    From a given set of fourier modes, compute the (binned) power spectrum
    with errors.

    modes   - (n,n,n) array (from ifftn of overdensity values)
    boxsize - size of box (e.g. in Mpc/h)

    returns

    kmid    - k values of power spec
    pspec   - estimate of power spec at k
    perr    - error on estimate
    """

    ngrid = modes.shape[0]
    assert(modes.shape==(ngrid, ngrid,ngrid))
    kmin, kmax, kbins, kvol = powerspec_bins(ngrid, boxsize)

    wts = square(modes.ravel().real) + square(modes.ravel().imag)

    v1 = bincount(kbins, weights=wts)
    powerspec = v1 * (1.0 / kvol)

    # work out error on power spectrum
    v2 = bincount(kbins, weights=square(wts))
    v0 = bincount(kbins)
    p_err = sqrt((v2*v0 - v1*v1)/(v0-1)) / kvol

    kmid_bins = 0.5 * (kmin+kmax)

    return kmid_bins, powerspec, p_err

def powerspec_bins(ngrid, boxsize):
    """
    find power spectrum bins to for a cubic grid of size ngrid^3 of fourier modes.
    Assumes the FFT convention of 0, ..., n/2, -n/2+1, ..., -1
    ngrid   - num cells on side of cube
    boxsize - size of the box in real space

    returns kmin, kmax, kbins, kvol
    kmin  - the lower bound of the bin
    kmax  - the upper bound of the bin
    kbins - index (0, ..., m) of the bin of each cell
    kvol  - the volume in k space (number of cells of that type divided by k vol of each cell)
    """

    mid = int(ngrid/2)
    # find the magnitude of the indices (i.e. ix**2+iy**2+iz**2 in the FFT convention)
    n1 = arange(ngrid)
    n1[1+mid:] -= ngrid
    n2 = square(n1)
    nmag = sqrt(add.outer(add.outer(n2, n2), n2)).ravel()

    nbins = (-1,) + tuple(arange(mid-1)+1.5) + (ngrid*2,)
    #print 'nbins', nbins
    kbins = digitize(nmag, nbins) - 1
    assert(kbins.min()==0)
    assert(kbins.max()==len(nbins)-2)

    # multiplier to go to k-space
    kmult = 2.0 * pi / boxsize

    kmin = (array(nbins) * kmult)[:-1]
    kmin[0] = 0

    kmax = (array(nbins) * kmult)[1:]
    kmax[-1] = mid * kmult * sqrt(3.0)

    kvol = bincount(kbins) * (kmult * kmult * kmult)
    return kmin, kmax, kbins, kvol



def calculate_power_spectrum(ionized_cells, density, GridSize, mode):

    if (mode == 0): # Calculating 21cm power spectrum
        print("Calculating the 21cm power spectrum.")
        meanion = np.mean(ionized_cells)
        Tb = (1.0 - ionized_cells)*density
        modes = ifftn(Tb)

    elif (mode == 1): # Calculating XHII power spectrum

        print("Calculating the ionized power spectrum.") 
        modes = ifftn(ionized_cells)
    else:
        print("Currently only support calculation of the 21cm power spectrum (mode = 0) and XHII power spectrum (mode = 1).")
        exit()   
 
    kmid_bins, powerspec, p_err = modes_to_pspec(modes, boxsize=AllVars.BoxSize/AllVars.Hubble_h)
    return (kmid_bins, powerspec, p_err)
    
    '''	
    Brightness = 27.0*(1.0-ionized_cells)*density*((1+ZZ[snapshot_idx])/10.0 * 0.15/(AllVars.Omega_m*AllVars.Hubble_h*AllVars.Hubble_h))**(1/2) * (AllVars.Omega_m*AllVars.Hubble_h*AllVars.Hubble_h / 0.023)

    Mean_Brightness = np.mean(Brightness)
    delta = Brightness/Mean_Brightness - 1.0	

    k_Bins, PowSpec, Error = Power_Spectrum(delta, GridSize)
   
    return (k_Bins , PowSpec, Error, Mean_Brightness) 
    '''         


##########################################################################

def check_fractions(ionized_fraction, HI_fraction_high, HI_fraction_low):
    for i in range(0, len(HI_fraction_high)):
        if (ionized_fraction >= HI_fraction_low[i] and ionized_fraction <= HI_fraction_high[i]):
            return i

    return -1

##########################################################################

##########

if __name__ == '__main__':

    ###########################   
    
    PlotScripts.Set_Params_Plot()

    ## Kali ## 

    #model_tags = [r"$f_\mathrm{esc} = 0.35$", r"$f_\mathrm{esc} \: \propto \: M_\mathrm{H}^{-1}$", r"$f_\mathrm{esc} \: \propto \: M_\mathrm{H}$", r"$f_\mathrm{esc} \: \propto \: \mathrm{Quasar \: Activity}$", r"$f_\mathrm{esc} \: \propto \: f_\mathrm{ej}$"]
    #model_tags = [r"$f_\mathrm{esc} = 0.35$", r"$f_\mathrm{esc} \: \propto \: M_\mathrm{H}^{-1}$", r"$f_\mathrm{esc} \: \propto \: M_\mathrm{H}$"]
    #model_tags = [r"$f_\mathrm{esc} = \mathrm{Constant}$", r"$f_\mathrm{esc} \: \propto \: f_\mathrm{ej}$", r"$df_\mathrm{esc} / dM_\mathrm{H} > 0$"]
    model_tags = [r"$f_\mathrm{esc} = \mathrm{Constant}$", r"$f_\mathrm{esc} \: \propto \: f_\mathrm{ej}$", r"$f_\mathrm{esc} \: \propto \: \mathrm{Quasar \: Activity}$"]
    #model_tags = [r"Ejectd", r"Anne MHPos", r"Anne MHNeg"]
    #model_tags = [r"$f_\mathrm{esc} = 0.35$", r"$f_\mathrm{esc} \: \propto \: f_\mathrm{ej}$", r"$f_\mathrm{esc} \: \propto \: \mathrm{Quasar \: Activity}$"]
    #model_tags = [r"$f_\mathrm{esc} = 0.35$", r"$f_\mathrm{esc} \: \propto \: M_\mathrm{H}$",r"$\mathrm{Quasar}$", r"$f_\mathrm{esc} = 0.4f_\mathrm{ej} + 0.2$"]
    #model_tags = [r"$f_\mathrm{esc} = 0.35$", r"$f_\mathrm{esc} \: \propto \: M_\mathrm{H}$",r"$\mathrm{Quasar}$"]

    output_tags = [r"Constant", r"Ejected", "MHPos"]
#    output_tags = [r"Delayed_Constant", r"Delayed_Quasar", r"Delayed_fej"]
#    output_tags = [r"Delayed_Constant", r"Delayed_MHneg", r"Delayed_MHPos", r"Delayed_Quasar", r"Delayed_fej"]
#    output_tags = [r"Constant", r"MHPos", r"Quasar"]
#    output_tags = [r"Densfield_full_Snap20"]
#    output_tags = [r"densgrid_snapshot020_full"]
 
    number_models = 3

    simulation_model1 = 6 # Which simulation are we using?
    # 0 : Mysim (Manodeep's original simulation run).
    # 2 : Tiamat (down to z = 5).
    # 4 : Simfast21 (deprecated).
    # 5 : Britton's. 
    # 6 : Kali

    model = 'Kali_ejected_quasar'

    GridSize_model1 = 256
        
    precision_model1 = 2 # 0 for int reading, 1 for float, 2 for double.

    ## IRA
    filepath_model1 = "/lustre/projects/p004_swin/jseiler/kali/grids/cifog/XHII_photHI1_sphere_fesc0.35"
    filepath_model2 = "/lustre/projects/p004_swin/jseiler/kali/grids/cifog/XHII_photHI1_fescMHneg"  
    filepath_model3 = "/lustre/projects/p004_swin/jseiler/kali/grids/cifog/XHII_photHI1_fescMHpos_sphere"
    filepath_model4 = "/lustre/projects/p004_swin/jseiler/kali/grids/cifog/XHII_photHI1_fescquasar_sphere"
    filepath_model5 = "/lustre/projects/p004_swin/jseiler/kali/grids/cifog/XHII_photHI1_sphere_fej"

    ## Delayed
    filepath_model6 = "/lustre/projects/p004_swin/jseiler/kali/grids/cifog/XHII_delayed_photHI1_sphere_fesc0.35"
    filepath_model7 = "/lustre/projects/p004_swin/jseiler/kali/grids/cifog/XHII_delayed_photHI1_fescMHneg_sphere"  
    filepath_model8 = "/lustre/projects/p004_swin/jseiler/kali/grids/cifog/XHII_delayed_photHI1_fescMHpos_sphere"
    filepath_model9 = "/lustre/projects/p004_swin/jseiler/kali/grids/cifog/XHII_delayed_photHI1_fescquasar_sphere"
    filepath_model10 = "/lustre/projects/p004_swin/jseiler/kali/grids/cifog/XHII_delayed_photHI1_sphere_fej"

    ## Anne
    filepath_model12 = "/lustre/projects/p004_swin/jseiler/kali/grids/cifog/XHII_delayed_photHI1_fescMHpos_anne_sphere"
    filepath_model13 = "/lustre/projects/p004_swin/jseiler/kali/grids/cifog/XHII_delayed_photHI1_ejected_new_sphere"
    

    ## IRA        
    filepath_nion_model1 = "/lustre/projects/p004_swin/jseiler/kali/grids/kali_fiducial_grid256_densfieldgridding_fesc0.35_HaloPartCut32_nionHI"    
    filepath_nion_model2 = "/lustre/projects/p004_swin/jseiler/kali/grids/Galaxies_kali_fiducial_grid256_densfieldgridding_MH_1.000e+08_0.70_1.000e+12_0.20_nionHI"
    filepath_nion_model3 = "/lustre/projects/p004_swin/jseiler/kali/grids/kali_fiducial_grid256_densfieldgridding_MH_1.000e+08_0.20_1.000e+12_0.70_nionHI"
    filepath_nion_model4 = "/lustre/projects/p004_swin/jseiler/kali/grids/kali_fiducial_grid256_densfieldgridding_quasar_0.20_1.00_1_nion_HI"
    filepath_nion_model5 = "/lustre/projects/p004_swin/jseiler/kali/grids/kali_fiducial_grid256_densfieldgridding_Ejected_alpha0.400beta0.200_nionHI"

    ## Delayed
    filepath_nion_model6 = "/lustre/projects/p004_swin/jseiler/kali/grids/kali_delayedSN_grid256_densfieldgridding_fesc0.35_HaloPartCut32_nionHI"    
    filepath_nion_model7 = "/lustre/projects/p004_swin/jseiler/kali/grids/kali_delayedSN_grid256_densfieldgridding_MH_1.000e+08_0.70_1.000e+12_0.20_nionHI"
    filepath_nion_model8 = "/lustre/projects/p004_swin/jseiler/kali/grids/kali_delayedSN_grid256_densfieldgridding_MH_1.000e+08_0.20_1.000e+12_0.70_nionHI"
    filepath_nion_model9 = "/lustre/projects/p004_swin/jseiler/kali/grids/kali_delayedSN_grid256_densfieldgridding_quasar_0.20_1.00_1_nion_HI"
    filepath_nion_model10 = "/lustre/projects/p004_swin/jseiler/kali/grids/kali_delayedSN_grid256_densfieldgridding_Ejected_alpha0.400beta0.200_nionHI"

    ## Anne
    filepath_nion_model11 = "/lustre/projects/p004_swin/jseiler/kali/grids/kali_delayedSN_grid256_densfieldgridding_Anne_MH_6.310e+09_0.99_3.160e+11_0.05_nionHI"
    filepath_nion_model12 = "/lustre/projects/p004_swin/jseiler/kali/grids/kali_delayedSN_grid256_densfieldgridding_Anne_MH_3.960e+09_0.00_3.160e+11_0.95_nionHI"
    filepath_nion_model13 = "/lustre/projects/p004_swin/jseiler/kali/grids/kali_delayedSN_grid256_densfieldgridding_Ejected_alpha1.000beta0.000_nionHI"

    filepath_density_model1 = "/lustre/projects/p134_swin/jseiler/kali/density_fields/averaged/"
        
    filepath_photofield_model1 = "/lustre/projects/p004_swin/jseiler/kali/grids/cifog/photHI_photHI1_fesc0.35"
    filepath_photofield_model2 = "/lustre/projects/p004_swin/jseiler/kali/grids/cifog/photHI_photHI1_fescMHneg"  
    

    simulation_norm = [simulation_model1, simulation_model1, simulation_model1, simulation_model1, simulation_model1]
    precision_array = [precision_model1, precision_model1, precision_model1, precision_model1, precision_model1]
    GridSize_array = [GridSize_model1, GridSize_model1, GridSize_model1, GridSize_model1, GridSize_model1]
    ionized_cells_filepath_array = [filepath_model6, filepath_model13, filepath_model9]
    #ionized_cells_filepath_array = [filepath_model1, filepath_model2, filepath_model3, filepath_model4, filepath_model5]
    #ionized_cells_filepath_array = [filepath_model6, filepath_model7, filepath_model8, filepath_model9, filepath_model10]
    #nion_filepath_array = [filepath_nion_model6, filepath_nion_model7, filepath_nion_model8, filepath_nion_model9, filepath_nion_model10] 
    nion_filepath_array = [filepath_nion_model6, filepath_nion_model13, filepath_nion_model9]
    #nion_filepath_array = [filepath_nion_model1, filepath_nion_model3, filepath_nion_model4, filepath_nion_model5] 
    #nion_filepath_array = [filepath_nion_model1, filepath_nion_model2, filepath_nion_model3, filepath_nion_model4, filepath_nion_model5] 
    density_filepath_array = [filepath_density_model1, filepath_density_model1, filepath_density_model1]
    photofield_filepath_array = [filepath_photofield_model1, filepath_photofield_model2]
    
    ###########################      
    
    cosmo = AllVars.Set_Params_Kali()

    OutputDir = "./ionization_plots/" + model + '/'
    if not os.path.exists(OutputDir):
        os.makedirs(OutputDir)

    #snaplist = np.arange(0, len(AllVars.SnapZ)) 
    #snaplist = np.arange(0, 70)
    snaplist = np.arange(27, 94) 
    #snaplist = np.arange(20, 50)

    ZZ = np.empty(len(snaplist))
    
    for i in range(0, len(snaplist)):
        ZZ[i] = AllVars.SnapZ[snaplist[i]]  
    z_index = 0 

    lowerZ = 0 
    upperZ = len(ZZ) 

    volume_frac_array = np.zeros((number_models, upperZ - lowerZ), dtype = np.float64) 
    mass_frac_array = np.zeros((number_models, upperZ - lowerZ), dtype = np.float64)
    nion_total_array = np.zeros((number_models, upperZ - lowerZ), dtype = np.float64) 
    photo_mean_array = np.zeros((number_models, upperZ - lowerZ), dtype = np.float64)
    photo_std_array = np.zeros((number_models, upperZ - lowerZ), dtype = np.float64)
    hoshen_array = np.zeros((number_models, upperZ - lowerZ), dtype = np.float64)

    wavenumber_array = [] 
    powerspectra_array = [] 
    powerspectra_error_array = [] 

    wavenumber_XHII_array = []
    powerspectra_XHII_array = []
    powerspectra_XHII_error_array = []        

    fraction_idx_array = []

    for model_number in range(0, number_models):
        wavenumber_array.append([])
        powerspectra_array.append([])
        powerspectra_error_array.append([])

        wavenumber_XHII_array.append([])
        powerspectra_XHII_array.append([])
        powerspectra_XHII_error_array.append([])

        fraction_idx_array.append([]) 

    #fractions_HI = [0.90, 0.01]
    #delta_HI = [0.01, 0.01]

    #fractions_HI = np.arange(0.75, 0.25, -0.05) 
    #delta_HI = np.full(len(fractions_HI), 0.025)
  
    fractions_HI = [0.75, 0.25]
    delta_HI = [0.03, 0.03]

    HI_fraction_high = np.add(fractions_HI, delta_HI) 
    HI_fraction_low= np.subtract(fractions_HI, delta_HI) 

    MC_ZZ = np.full((number_models, len(fractions_HI)), -1.0)
    MC_Snaps = np.zeros((number_models, len(fractions_HI)))

    do_MC = 0 
    calculate_MC = 0 # 0 to NOT calculate the MC bubbles, 1 to calculate them at the HI fractions specified by fractions_HI, 2 to calculate them at the snapshots given by calculate_MC_snaps. 
    plot_MC = 0 # 0 is nothing, 1 to plot bubble properties at specified HI fractions, 2 to plot them at specified snapshots.
    calculate_MC_snaps = np.arange(lowerZ, upperZ, 1)
    MC_ZZ_snaps = [ZZ[i] for i in calculate_MC_snaps] 
 
    calculate_power = 1 # 0 to NOT calculate the power spectra, 1 to calculate it. 

    do_hoshen = 0
    
    for snapshot_idx in range(lowerZ+rank, upperZ, size):

        do_power_array = np.zeros((number_models))	
        ionized_cells_array = []
        density_array = [] 
        nion_array = [] 
        photofield_array = [] 

        print("==============================")
        print("REDSHIFT {0:.2f}".format(ZZ[snapshot_idx]))	
        print("==============================")

        if (snapshot_idx < 10):
            number_tag_anne = '_0%d' %(snapshot_idx)
            number_tag_mine = '_00%d' %(snapshot_idx)
        else:
            number_tag_anne = '_%d' %(snapshot_idx)
            number_tag_mine = '_0%d' %(snapshot_idx)

        ##########################################################################

        for model_number in range(0, number_models): 

            if (simulation_norm[model_number] == 0):
                cosmo = AllVars.Set_Params_Mysim()
            elif (simulation_norm[model_number] == 2):
                cosmo = AllVars.Set_Params_Tiamat()
            elif (simulation_norm[model_number] == 4):
                print("Attempted to set the parameters to Simfast21.  This is so old.")
                exit()
                cosmo = AllVars.Set_Params_Simfast21()
            elif (simulation_norm[model_number] == 5):
                cosmo = AllVars.Set_Params_Britton()
            elif (simulation_norm[model_number] == 6):
                cosmo = AllVars.Set_Params_Kali()
            else:
                print("The current value for simulation norm ({0}) is not implmented yet.".format(simulation_norm[model_number]))
                exit()
            GridSize_model = GridSize_array[model_number]
            precision_model = precision_array[model_number]

                #######################################
                ##### Reading in all the files. ####### 
                #######################################


            ionized_cells_path = ionized_cells_filepath_array[model_number] + "_%02d" %(snapshot_idx) 
            ionized_cells_array.append(ReadScripts.read_binary_grid(ionized_cells_path, GridSize_model, precision_model)) 

            nion_path = "{0}_{1:03d}".format(nion_filepath_array[model_number], snaplist[snapshot_idx])
            nion_array.append(ReadScripts.read_binary_grid(nion_path, GridSize_model, 1) * 1.0e50) 

            density_path = "{0}snap{1:03d}.dens.dat".format(density_filepath_array[model_number], snaplist[snapshot_idx])  
            #density_path = "{0}.dens.dat".format(density_filepath_array[model_number])  
            density_array.append(ReadScripts.read_binary_grid(density_path, GridSize_model, precision_array[model_number]))

            #photofield_path = "{0}_{1:02d}".format(photofield_filepath_array[model_number], snaplist[snapshot_idx]) 
            #photofield_array.append(ReadScripts.read_binary_grid(photofield_path, GridSize_model, 2)) 

            tag = "{0}_{1:03d}".format(output_tags[model_number], snaplist[snapshot_idx])                  
            #plot_single(ZZ[snapshot_idx], ionized_cells_array[model_number], GridSize_array[model_number], simulation_norm[model_number], OutputDir, tag) 

            volume_frac_array[model_number][snapshot_idx] = calculate_volume_frac(ionized_cells_array[model_number], GridSize_array[model_number])
            nion_total_array[model_number][snapshot_idx] = calculate_total_nion(model_tags[model_number], nion_array[model_number])
            volume_frac_model = volume_frac_array[model_number][snapshot_idx]	

            if(do_hoshen == 1):
                hoshen_array[model_number][snapshot_idx] = hoshen_kopelman(ionized_cells_array[model_number])	
        
            #plot_nionfield(ZZ[snapshot_idx], nion_array[model_number], OutputDir, "Nion_" + output_tags[model_number] + '_' + str(snaplist[snapshot_idx]))
            #density_array[model_number] = density_array[model_number] + 1
            #plot_density(ZZ[snapshot_idx], density_array[model_number], OutputDir, "Density_" + output_tags[model_number] + '_' + str(snaplist[snapshot_idx]))
            #plot_density_numbers(ZZ[i], density_array[model_number], OutputDir, "DensityNumbers" + str(i))

            fraction_idx = check_fractions(volume_frac_model, HI_fraction_high, HI_fraction_low) # Checks the current ionization fraction with the fractions that we wanted to do extra stuff at.
            

            if(fraction_idx != -1): 
                if(MC_ZZ[model_number, fraction_idx] == -1):	
                    MC_ZZ[model_number, fraction_idx] = ZZ[snapshot_idx]
                    MC_Snaps[model_number, fraction_idx] = snapshot_idx

                    print("Model {0} reached x_HI = {1:.3f} at z = {2:.2f}".format(model_number, 1 - fractions_HI[fraction_idx], ZZ[snapshot_idx]))

                    do_MC = 1
                    do_power_array[model_number] = 1
          
            if((do_MC == 1 and calculate_MC == 1) or (calculate_MC == 2 and snapshot_idx == calculate_MC_snaps[snapshot_idx])):
         
                print("snapshot_idx {0}".format(snapshot_idx))
                print("Model_number {0}".format(model_number))
                print("ZZ[snapshot_idx {0}".format(ZZ[snapshot_idx]))
                #print("ionized_cells_array[model_number] {0}".format(ionized_cells_array[model_number]))
                print("GridSize_array[model_number] {0}".format(GridSize_array[model_number]))	
        
                calculate_bubble_MC(ZZ[snapshot_idx], ionized_cells_array[model_number], GridSize_array[model_number], OutputDir, output_tags[model_number])
                do_MC = 0

            if(do_power_array[model_number] == 1 and calculate_power == 1):
            
                T0 = T_naught(ZZ[snapshot_idx], AllVars.Omega_b, AllVars.Hubble_h, AllVars.Omega_m)
                #tmp_k, tmp_PowSpec, tmp_Error, Mean_Brightness = calculate_power_spectrum(ionized_cells_array[model_number], density_array[model_number], GridSize_array[model_number])
                tmp_k, tmp_PowSpec, tmp_Error = calculate_power_spectrum(ionized_cells_array[model_number], density_array[model_number], GridSize_array[model_number], 0)
                wavenumber_array[model_number].append(tmp_k)
                #powerspectra_array[model_number].append(Mean_Brightness**2 * tmp_PowSpec * tmp_k**3 * 2.0 * np.pi * np.pi)
                powerspectra_array[model_number].append(T0*T0 * tmp_PowSpec * tmp_k**3 * 2.0 * np.pi)
                powerspectra_error_array[model_number].append(tmp_Error)
                fraction_idx_array[model_number].append(fraction_idx)

                tmp_k, tmp_PowSpec, tmp_Error = calculate_power_spectrum(ionized_cells_array[model_number], density_array[model_number], GridSize_array[model_number], 1)
                wavenumber_XHII_array[model_number].append(tmp_k)
                powerspectra_XHII_array[model_number].append(tmp_PowSpec * tmp_k**3 * 2.0 * np.pi)
                powerspectra_XHII_error_array[model_number].append(tmp_Error)

                do_power_array[model_number] = 0

            #plot_photofield(ZZ[i], photofield_array[model_number], OutputDir, "PhotHIField_" + output_tags[model_number] + str(i))

            #photo_mean_array[model_number][snapshot_idx] = np.mean(photofield_array[model_number][photofield_array[model_number] != 0])
            #photo_std_array[model_number][snapshot_idx] = np.std(photofield_array[model_number][photofield_array[model_number] != 0])
         
        #print "This snapshot has index %d with lookback time %.4f (Gyr)" %(i, cosmo.lookback_time(ZZ[i]).value)

    
    #for model_number in range(0, number_models):
    #    if (MC_ZZ[model_number][-1] < 5 or MC_ZZ[model_number][-1] > 1000):
    #        MC_ZZ[model_number][-1] = ZZ[-1]


    volume_frac_array_final = np.zeros((number_models, upperZ - lowerZ), dtype = np.float64) 
    comm.Reduce([volume_frac_array, MPI.DOUBLE], [volume_frac_array_final, MPI.DOUBLE], op = MPI.SUM, root = 0)

    mass_frac_array_final = np.zeros((number_models, upperZ - lowerZ), dtype = np.float64)
    comm.Reduce([mass_frac_array, MPI.DOUBLE], [mass_frac_array_final, MPI.DOUBLE], op = MPI.SUM, root = 0)
   
    nion_total_array_final = np.zeros((number_models, upperZ - lowerZ), dtype = np.float64)
    comm.Reduce([nion_total_array, MPI.DOUBLE], [nion_total_array_final, MPI.DOUBLE], op = MPI.SUM, root = 0)

    photo_mean_array_final = np.zeros((number_models, upperZ - lowerZ), dtype = np.float64)
    comm.Reduce([photo_mean_array, MPI.DOUBLE], [photo_mean_array_final, MPI.DOUBLE], op = MPI.SUM, root = 0)

    photo_std_array_final = np.zeros((number_models, upperZ - lowerZ), dtype = np.float64)
    comm.Reduce([photo_std_array, MPI.DOUBLE], [photo_std_array_final, MPI.DOUBLE], op = MPI.SUM, root = 0)
    
    hoshen_array_final = np.zeros((number_models, upperZ - lowerZ), dtype = np.float64)
    comm.Reduce([hoshen_array, MPI.DOUBLE], [hoshen_array_final, MPI.DOUBLE], op = MPI.SUM, root = 0)

    MC_Snaps_final = np.zeros((number_models, len(fractions_HI)))
    comm.Reduce([MC_Snaps, MPI.DOUBLE], [MC_Snaps_final, MPI.DOUBLE], op = MPI.MAX, root = 0)

    MC_Snaps_final = np.zeros((number_models, len(fractions_HI)))
    comm.Reduce([MC_Snaps, MPI.DOUBLE], [MC_Snaps_final, MPI.DOUBLE], op = MPI.MAX, root = 0)
   
    MC_ZZ_final = np.zeros((number_models, len(fractions_HI)))
    comm.Reduce([MC_ZZ, MPI.DOUBLE], [MC_ZZ_final, MPI.DOUBLE], op = MPI.MAX, root = 0)

    if (rank == 0): 
        #plot_global_frac(ZZ, mass_frac_array_final, volume_frac_array_final, 1, model_tags, OutputDir, "GlobalFraction_time")
        #plot_total_nion(ZZ, nion_total_array_final, simulation_norm, 1, model_tags, OutputDir, "Nion_z")
        #plot_optical_depth(ZZ, volume_frac_array_final, model_tags, OutputDir, "OpticalDepth")

        #plot_nine_panel_slices(ZZ, ionized_cells_filepath_array, GridSize_array, simulation_norm, MC_Snaps_final, fractions_HI, model_tags, OutputDir, "3PanelSlice")

        if(plot_MC == 1):
            plotting_MC_ZZ = np.zeros(np.shape(MC_ZZ_final))
            plotting_HI = np.zeros(np.shape(fractions_HI))
            for i in range((number_models)):
                for j in range(1, len(MC_ZZ_final)+1):
                    print(plotting_MC_ZZ)
                    print(MC_ZZ_final[i][-j%len(MC_ZZ_final)])
                    plotting_MC_ZZ[i][j-1] = MC_ZZ_final[i][-j%len(MC_ZZ_final)]
#            plotting_MC_ZZ = MC_ZZ_final
            plotting_HI = fractions_HI
        elif(plot_MC == 2):
            plotting_MC_ZZ = MC_ZZ_snaps 
            plotting_HI = volume_frac_array 

    #plot_bubble_MC(plotting_MC_ZZ, plotting_HI, simulation_norm, model_tags, output_tags, GridSize_array, OutputDir, "BubbleSizes") 
    plot_power(fractions_HI, wavenumber_array, powerspectra_array, powerspectra_error_array, wavenumber_XHII_array, powerspectra_XHII_array, powerspectra_XHII_error_array, fraction_idx_array, model_tags, OutputDir, "PowerSpectrum_Anne")

    #print "Duration of reionization for Model %s is %.4f Myr (%.4f Gyr - %.4f Gyr)" %(model_tags[0], (cosmo.lookback_time(MC_ZZ[0][0]).value - cosmo.lookback_time(MC_ZZ[0][-1]).value) * 1.0e3, cosmo.lookback_time(MC_ZZ[0][0]).value, cosmo.lookback_time(MC_ZZ[0][-1]).value)
    #print "Duration of reionization for Model %s is %.4f Myr (%.4f Gyr - %.4f Gyr)" %(model_tags[1], (cosmo.lookback_time(MC_ZZ[1][0]).value - cosmo.lookback_time(MC_ZZ[1][-1]).value) * 1.0e3, cosmo.lookback_time(MC_ZZ[1][0]).value, cosmo.lookback_time(MC_ZZ[1][-1]).value)
    #print "Duration of reionization for Model %s is %.4f Myr (%.4f Gyr - %.4f Gyr)" %(model_tags[2], (cosmo.lookback_time(MC_ZZ[2][0]).value - cosmo.lookback_time(MC_ZZ[2][-1]).value) * 1.0e3, cosmo.lookback_time(MC_ZZ[2][0]).value, cosmo.lookback_time(MC_ZZ[2][-1]).value)
    #plot_redshifts(redshift_array_model1, ZZ, lowerZ, upperZ, OutputDir, "zreiond3_" + output_tags[0])
    #plot_redshifts(redshift_array_model2, ZZ, lowerZ, upperZ, OutputDir, "zreiond3_" + output_tags[1])

    #save_redshifts(redshift_array_model1, OutputDir, "ReionizationRedshift_" + output_tags[0])
    #save_redshifts(redshift_array_model2, OutputDir, "ReionizationRedshift_" + output_tags[1])
    #save_redshifts(redshift_array_model3, OutputDir, "ReionizationRedshift_" + output_tags[2])

    #plot_hoshen(hoshen_array, volume_frac_array, model_tags, OutputDir, "Hoshen")
    #find_max_bubble(plotting_MC_ZZ, plotting_HI, plot_MC, volume_frac_array,  model_tags, output_tags, GridSize_array, OutputDir, "MaxBubble_withNB_z_ylog")

    #photon_baryon(lowerZ, upperZ, ZZ, [nion_total_model1, nion_total_model2], Hubble_h, OB, Y, model_tags, OutputDir, "nion_total2")
    #analytic_HII(nion_total_model1, ZZ, upperZ, snaplist, OutputDir, "Q_Analytic")	
    #plot_photo_mean(ZZ, photo_mean_array, photo_std_array, model_tags, OutputDir, "Mean_Photo")

    #plot_deltat_deltax(ZZ, volume_frac_array, model_tags, OutputDir, "ReionizationSpeed")
    #plot_deltat_deltaN(ZZ, nion_array, model_tags, OutputDir, "NionSpeed")
    #plot_combined_global_nion(ZZ, nion_array, volume_frac_array, model_tags, OutputDir, "Combined")


    #plot_nine_panel_photHI(ZZ, [filepath_photofield_model1, filepath_photofield_model2, filepath_photofield_model3], [GridSize_model1, GridSize_model2, GridSize_model3], MC_Snaps, fractions_HI, model_tags, OutputDir, "9PanelSlices_PhotHI")
    
