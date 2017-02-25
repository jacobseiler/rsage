import matplotlib
matplotlib.use('Agg')

import os
import numpy as np
import pylab as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from numpy import *
from random import sample, seed
from os.path import getsize as getFileSize
import math
import random
import csv
from io import StringIO
#np.set_printoptions(threshold=np.nan)
from collections import Counter
from matplotlib.colors import LogNorm
import time
from scipy.ndimage.filters import generic_filter as gf
from matplotlib.ticker import MultipleLocator


import matplotlib.ticker as mtick
import PlotScripts
import ReadScripts
import AllVars


G = 6.674e-11 # Gravitational constant.
mp_kg = 1.673e-27 # Mass proton in kg.
Mpc_per_km = 3.24e-20 # How many Mpc in a km.
Mpc_m = 3.0857e22 # Mpc in m
Cubic_cm_Cubic_Mpc = 2.93799e73 # cm^-3 to Mpc_3
km_m = 1e3 # km in m
sec_per_Myr = 3.154e13 # How many seconds are in a Myr. 
nb = 2.5e-7 # Number density of baryons in cm^-3.

def linear_growth(z, OM, OL):
        D = (OM*((1+z)**3) / (OM*(1+z) - (OM + OL - 1)*((1+z)**2) + OM))**(4/7)
        return D

def cosmo_density(z, Hubble_h):
        rho = (3*100*100*Hubble_h*Hubble_h*km_m*km_m)/(8*np.pi*G*Mpc_m*Mpc_m)*(1+z)**3
        return rho

## Number density of hydrogen in m^-3##
def n_HI(z, Hubble_h, OB, Y):
	n_HI = cosmo_density(z, Hubble_h)*OB*(1-Y)/mp_kg 
	return n_HI


def stromgren_sphere(Q, T, z, Y):
	# Takes the photon flux (Q, photons per second) with an electron temperature (T, Kelvin) and helium fraction (Y)
	# And returns the radius of a Stromgren sphere at redshift z in Mpc/h.

	Beta = 2.e-16*(T_e)**(-3./4.) # Recombination rate (m^3 s^-1).
	num_HI = n_HI(z, Hubble_h, OB, Y) # Number density of hydrogen (m^3).

	R = ((3.*Q)/(4.*np.pi*num_HI*num_HI*Beta)) ** (1./3.)
	R /= Mpc_m # Stromgren radius in units of Mpc.
	R *= Hubble_h

	print "Stromgren Radius for a flux of %.3e photons per second at a redshift %.3f and an electron temperature of %.3e is %.3e Mpc/h" %(Q, z, T_e, R)
	return R

def T_naught(z, OB, h, OM):
	# T0 brightness temperature from Hutter et al (2016)
	# Output in units of mK.

	T0 = 28.5 * ((1.0+z)/10.0)**(0.5) * OB/0.042 * h/0.73 * (OM/0.24)**(-0.5) 
	return T0

matplotlib.rcdefaults()
plt.rc('axes', color_cycle=['k', 'b', 'r', 'g', 'm', 'c','y', '0.5'], labelsize='x-large')
plt.rc('lines', linewidth='2.0')
# plt.rc('font', variant='monospace')
plt.rc('legend', numpoints=1, fontsize='x-large')
plt.rc('text', usetex=True)

label_size = 12
extra_size = 2
plt.rc('xtick', labelsize=label_size)
plt.rc('ytick', labelsize=label_size)
plt.rc('text', usetex=True)
tick_interval = 0.25
np.set_printoptions(formatter={'float': lambda x: "{0:0.10e}".format(x)})
 
colours = ['r', 'b', 'g', 'c', 'm', 'k', 'r', 'b', 'g', 'c', 'm', 'k']
markers = ['x', 'o', 'x', 'o', 'D', 's']
linestyles = ['-', '--', '-.', ':']

Output_Format = ".png"

GridSize = 256
BoxSize = 100.0 # Mpc/h
Ratio = BoxSize/GridSize 
Hubble_h = 0.678 
BaryonFrac = 0.17 # Baryon Fraction 
OM = 0.308 # Omega Matter
OB = OM*BaryonFrac # Omega Baryons
OL = 1 - OM # Omega Lambda
Y = 0.24 # Helium Fraction
cut_slice = 44

AllVars.Set_Params_Mysim()
AllVars.Set_Constants()

output_format = ".png"
output_label = "halos_most_new"

def plot_single(z, filepath, ionized_cells, OutputDir, output_tag):
# Plots the ionization bubbles for a single set of cells.

## Input ##
# z is the redshift we're plotting at, filepath is the path of the input data.
# ionized_cells is the array that contains the ionized cells, OutputDir is the path of the output directory.

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*_*z*.*output_format*'

	print
	print "Plotting ionized bubbles for file %s." %(filepath)
	
	ionized_cells = 1 - ionized_cells

	ax = plt.subplot(111)

	im = ax.imshow(ionized_cells[:,:,cut_slice:cut_slice+1].mean(axis = -1), interpolation='none', origin='low', extent =[0,BoxSize,0,BoxSize], vmin = 0, vmax = 1, cmap = 'afmhot_r')

	cbar = plt.colorbar(im, ax = ax)
	cbar.set_label(r'$1 - Q$')				
				    
	ax.set_xlabel(r'$\mathrm{x}  (h^{-1}Mpc)$')  
	ax.set_ylabel(r'$\mathrm{y}  (h^{-1}Mpc)$')  

	ax.set_xlim([0.0, BoxSize]) 
	ax.set_ylim([0.0, BoxSize])

	global_frac = sum(ionized_cells[ionized_cells != 1])/float(GridSize**3)

	title = r"z = %.3f, $\langle x_{HI} = %.3f \rangle$" %(z, 1-global_frac)
	ax.set_title(title)

	print "Number of ionized cells is %d giving an ionized Hydrogen fraction %.3f" %(len(ionized_cells[ionized_cells != 1]), global_frac) 

	#ax.text(52.5, 60, r'z = %.1f' %z, color = 'white', fontsize = 15)

	#outputFile = OutputDir + output_tag + '_z' + str(z) + output_format 
	outputFile = OutputDir + output_tag + output_format 
	plt.savefig(outputFile)  # Save the figure
	print 'Saved file to', outputFile
			
	plt.close()
	
############

def plot_comparison(z, ionized_cells_model1, ionized_cells_model2, nion_model1, nion_model2, OutputDir, model_tags, output_tag):
# Plots ionization bubbles for two realizaiton of the box.  These bubbles will be overlaid on top of each other to allow comparison between the two.

## Input ##
# z is the redshift we're plotting at.
# ionized_cells_model1 (and ionized_cells_model2) are the arrays that contains the ionized cells.
# OutputDir is the path of the output directory.

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*_*z*.*output_format*'

	print
	print "Plotting a comparison between the ionized bubbles." 	

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
	print 'Saved file to', outputFile
			
	plt.close()

##############

def plot_sources(z, filepath, nion_tmp, OutputDir, output_tag): 
# Plots the location of the ionizing sources.  At the moment this function only plots if there is a source in the cell and not its photon count.

## Input ##
# z is the redshift we're plotting at. nion is the array that contains the ionizing photons

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*_*z*.*output_format*'

	print
	print "Plotting ionizing sources for file %s" %(filepath)
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
	print 'Saved file to', outputFile
			
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

	print
	print "Plotting a difference map between the models."

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
	print 'Saved file to', outputFile
			
	plt.close()	

###########

def calculate_volume_frac(ionized_cells):
# This function will calculate (and return!) the fraction of cells that are ionized.  
# This function was separated from 'plot_global_frac' as the global fractions for each redshifts must be calculated before we can plot them.

## Input ##
# ionized_cells is the array that contains the ionization data.

## Output ##
# The fraction of HI in the grid.

	HI = 1 - sum(ionized_cells)/float(GridSize**3)
        print "There is %.4f fraction of HII." %(HI)

	return HI

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

        print 
        print "Mass averaged HII fraction is %.4f" %(HI)

	return HI

###########


def calculate_total_nion(model_name, nion):
# This function will calculate (and return!) the total number of ionizing photons (per second) in the grid.
# This function was separated from 'plot_total_nion' as total nion for each redshifts must be calculated before we can plot them.
# ionized_cells is the array that contains the ionization data.

# Output: The total number of ionizing photons (per second) in the grid.

#        nion_total = sum(nion)/Hubble_h/(BoxSize**3/Hubble_h**3)
	nion_total = sum(nion)

	print
	print "The total number of ionizing photons (per second) for the %s Model is %.4e." %(model_name, nion_total)

	return nion_total


##########




	



#########


def plot_global_frac(ZZ, mass_frac, volume_frac, labels, OutputDir, output_tag):
# This function will plot the global HI fraction as a function of redshift. The function accepts an array of arrays for the global HI fraction to allow plotting of multiple realizations on the same graph.

## Input ##
# ZZ is an array containing the redshifts we plot over.
# global_frac is an array of arrays where each array contains the global HI fraction at each redshift.
# labels is an array that contains the labels for each realization we are plotting.

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*_*z*.*output_format*'


	print 
	print "Plotting the global fraction of HI."

	ax = plt.subplot(111)

	colors = ['r', 'b', 'g', 'c', 'm']
	markers = ['x', 'o', '^', 's', 'D']

        top_label1, = ax.plot(np.NaN, np.NaN, '-', color = 'none', label = labels[0])
	volume_label1 = ax.scatter(ZZ[0:len(volume_frac[0])], volume_frac[0], color = colors[0], marker = 'x', label = 'Volume Averaged')
        mass_label1 = ax.scatter(ZZ[0:len(mass_frac[0])], mass_frac[0], color = colors[0], marker = 'o', label = 'Mass Averaged')

        top_label2, = ax.plot(np.NaN, np.NaN, '-', color = 'none', label = labels[1])
	volume_label2 = ax.scatter(ZZ[0:len(volume_frac[1])], volume_frac[1], color = colors[1], marker = 'x', label = 'Volume Averaged')
        mass_label2 = ax.scatter(ZZ[0:len(mass_frac[1])], mass_frac[1], color = colors[1], marker = 'o', label = 'Mass Averaged')
        
	ax.set_xlabel(r"$z$")
	ax.set_ylabel(r"$\langle x_{HI}\rangle$")

        labels = [top_label1, volume_label1, mass_label1, top_label2, volume_label1, mass_label2]   
    
#        leg = plt.legend([top_label1, volume_label1, mass_label1, top_label2, volume_label2, mass_label2], [H.get_label() for H in labels])#, loc=2, numpoints=1, labelspacing=0.1)
        leg = plt.legend(handles = [top_label1, volume_label1, mass_label1, top_label2, volume_label2, mass_label2], loc = 4, scatterpoints=1, labelspacing = 0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')

    	ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    	ax.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
	outputFile = OutputDir + output_tag + output_format 			
	plt.savefig(outputFile)  # Save the figure
	print 'Saved file to', outputFile
			
	plt.close()	

##########

def plot_photo_mean(ZZ, Photo_Mean, Photo_Std, labels, OutputDir, output_tag):
	
	ax = plt.subplot(111)
	for p in xrange(0, len(Photo_Mean)):

		plt.scatter(ZZ[0:len(Photo_Mean[p])], np.log10(Photo_Mean[p]), color = colours[p], marker = markers[p], label = labels[p])

	plt.xlabel(r"$\mathrm{z}$", size = label_size + extra_size)
        plt.ylabel(r"$\mathrm{log}_{10} \: \langle \Gamma_\mathrm{HI}\rangle \: \: [\mathrm{s}^{-1}]$")

        leg = plt.legend(loc='upper right', numpoints=1,
                         labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')


	plt.axis([6, 12, -16, -9.5])
    	ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    	ax.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))

	outputFile = OutputDir + output_tag + output_format 			
	plt.savefig(outputFile)  # Save the figure
	print 'Saved file to', outputFile
			
	plt.close()	


#########




def plot_total_nion(ZZ, total_nion, labels, OutputDir, output_tag):
# This function will plot the total number of ionizing photons as a function of redshift. The function accepts an array of arrays for the total number of ionizing photons to allow plotting of multiple realizations on the same graph.

## Input ##
# ZZ is an array containing the redshifts we plot over.
# total_nion is an array of arrays where each array contains the total number of ionizing photons at each redshift.
# labels is an array that contains the labels for each realization we are plotting.

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*_*z*.*output_format*'

	print 
	print "Plotting the total number of ionizing photons" 

	colors = ['r', 'b', 'g', 'c', 'm']


	ax = plt.subplot(111)

	for p in xrange(0, len(total_nion)):

		plt.scatter(ZZ[0:len(total_nion[p])], np.log10(total_nion[p]), color = colors[p], marker = markers[p], label = labels[p])

	plt.xlabel(r"$\mathrm{z}$")
        plt.ylabel(r"$\mathrm{log}_{10} \: \dot{N}_{\mathrm{HI}} \: \: [\mathrm{s}^{-1}]$")

        leg = plt.legend(loc='upper right', numpoints=1,
                         labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')


    	ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    	ax.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))

	outputFile = OutputDir + output_tag + output_format 			
	plt.savefig(outputFile)  # Save the figure
	print 'Saved file to', outputFile
			
	plt.close()	

##########


def plot_density(z, density, OutputDir, output_tag):
# Plots a slice of the density field.

## Input ##
# z is the redshift we're plotting at.
# density is the array that contains the overdensity of each cell.

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*_*z*.*output_format*'

	print
	print "Plotting density slice."
	print "The minimum density for redshift %.3f was %.3f." %(z, np.amin(density)) 
	print "The maximum density for redshift %.3f was %.3f." %(z, np.amax(density)) 

	print "The minimum density index for redshift %.3f was %.3f." %(z, np.argmin(np.ravel(density), axis = 0))
	print "The maximum density index for redshift %.3f was %.3f." %(z, np.argmax(np.ravel(density), axis = 0))

	print len(density[density == 0])

	ax = plt.subplot(111)
	
	im = ax.imshow(density[:,:,cut_slice:cut_slice+1].mean(axis = -1), interpolation='bilinear', origin='low', extent =[0,BoxSize,0,BoxSize], cmap = 'Purples', norm = colors.LogNorm(vmin = 0.12, vmax = 50)) 

	cbar = plt.colorbar(im, ax = ax)
	cbar.set_label(r'$\rho/\langle \rho \rangle$')
				    
	ax.set_xlabel(r'$\mathrm{x}  (h^{-1}Mpc)$')  
	ax.set_ylabel(r'$\mathrm{y}  (h^{-1}Mpc)$')  

	ax.set_xlim([0.0, BoxSize]) 
	ax.set_ylim([0.0, BoxSize])

	title = r"$z = %.3f$" %(z)
        ax.set_title(title)
	#ax.text(52.5, 60, r'z = %.1f' %ZZ[i], color = 'black', fontsize = 15)

	outputFile = OutputDir + output_tag + output_format 		
	plt.savefig(outputFile)  # Save the figure
	print 'Saved file to', outputFile
			
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

	print
	print "Plotting photoionization slice."
	print "The minimum photoionization rate for redshift %.3f was %.3f." %(z, photo_min)
	print "The maximum photoionization rate for redshift %.3f was %.3f." %(z, photo_max) 

	print "The minimum photoionization rate index for redshift %.3f was %.3f." %(z, np.argmin(np.ravel(photofield), axis = 0))
	print "The maximum photoionization rate index for redshift %.3f was %.3f." %(z, np.argmax(np.ravel(photofield), axis = 0))

	photo_slice = photofield[:,:,cut_slice:cut_slice+1].mean(axis = -1)

	ax = plt.subplot(111)
	
	im = ax.imshow(np.log10(photo_slice), interpolation='bilinear', origin='low', extent =[0,BoxSize,0,BoxSize], cmap = 'Purples', vmin = -16.22, vmax = -6.8) 

	cbar = plt.colorbar(im, ax = ax)
	cbar.set_label(r'$\mathrm{log}_{10}\Gamma [\mathrm{s}^{-1}]$')
				    
	ax.set_xlabel(r'$\mathrm{x}  (h^{-1}Mpc)$')  
	ax.set_ylabel(r'$\mathrm{y}  (h^{-1}Mpc)$')  

	ax.set_xlim([0.0, BoxSize]) 
	ax.set_ylim([0.0, BoxSize])

	title = r"$z = %.3f$" %(z)
        ax.set_title(title)
	#ax.text(52.5, 60, r'z = %.1f' %ZZ[i], color = 'black', fontsize = 15)

	outputFile = OutputDir + output_tag + output_format 		
	plt.savefig(outputFile)  # Save the figure
	print 'Saved file to', outputFile
			
	plt.close()	

############

def plot_nionfield(z, nion, OutputDir, output_tag):
# Plots a slice of the ionizing photon field at a specified redshift.

## Input ##
# z is the redshift we're plotting at.
# photofield is the photoionization rate of each cell. 

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*_*z*.*output_format*'

	nion_slice = nion[:,:,cut_slice:cut_slice+50].mean(axis = -1)
	
	nion_min = np.amin(np.log10(nion[nion > 1e20]))
	nion_max = np.amax(np.log10(nion[nion > 1e20]))

	print
	print "Plotting photoionization slice."
	print "The minimum ionizing photons rate for redshift %.3f was %.3f." %(z, nion_min)
	print "The maximum ionizing photons rate for redshift %.3f was %.3f." %(z, nion_max) 

	print "The minimum ionizing photons rate index for redshift %.3f was %.3f." %(z, np.argmin(np.ravel(nion), axis = 0))
	print "The maximum ionizing photons rate index for redshift %.3f was %.3f." %(z, np.argmax(np.ravel(nion), axis = 0))

	ax = plt.subplot(111)
	
	im = ax.imshow(np.log10(nion_slice), interpolation='bilinear', origin='low', extent =[0,BoxSize,0,BoxSize], cmap = 'Purples', vmin = nion_min, vmax = nion_max) 

	cbar = plt.colorbar(im, ax = ax)
	cbar.set_label(r'$\mathrm{log}_{10}N_{\gamma} [\mathrm{s}^{-1}]$')
				    
	ax.set_xlabel(r'$\mathrm{x}  (h^{-1}Mpc)$')  
	ax.set_ylabel(r'$\mathrm{y}  (h^{-1}Mpc)$')  

	ax.set_xlim([0.0, BoxSize]) 
	ax.set_ylim([0.0, BoxSize])

	title = r"$z = %.3f$" %(z)
        ax.set_title(title)
	#ax.text(52.5, 60, r'z = %.1f' %ZZ[i], color = 'black', fontsize = 15)

	outputFile = OutputDir + output_tag + output_format 		
	plt.savefig(outputFile)  # Save the figure
	print 'Saved file to', outputFile
			
	plt.close()	

############

def plot_density_numbers(z, density, OutputDir, output_tag):
# Plot the numerical value of the overdensity and highlight high values of the overdensities.  Useful if we want to investigate the values around some point.

## Input ##
# z is the redshift we're plotting at.
# density is the array that contains the overdensity of each cell.

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*_*z*.*output_format*'

	# These values will give the bounds of our zoom in. Numerical values are in cMpc.
	#lower_x = int(55.5*GridSize/BoxSize)
	#upper_x = int(62.5*GridSize/BoxSize)
	#lower_y = int(42.5*GridSize/BoxSize)
	#upper_y = int(55.5*GridSize/BoxSize) 

	#lower_x = 64 
	#upper_x = 78 
	#lower_y = 50 
	#upper_y = 64 

	lower_x = 64 
	upper_x = 78 
	lower_y = 45 
	upper_y = 58 


	density[density == 0] = 1
	density = np.log10(density)

	step = 1 # The interval we want to plot the values at.
	density_cut = density[:,:,cut_slice:cut_slice+1].mean(axis = -1)
	density_cut = density_cut[lower_x:upper_x:step,lower_y:upper_y:step]

        print np.amax(density_cut) 

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

	ax.set_xlim([(lower_y-1)*BoxSize/GridSize, (upper_y)*BoxSize/GridSize]) 
	ax.set_ylim([(lower_x-1)*BoxSize/GridSize, (upper_x)*BoxSize/GridSize])

	title = r"$z = %.3f$" %(z)

        ax.set_title(title)

	outputFile = OutputDir + output_tag + output_format 
	plt.savefig(outputFile)  # Save the figure
	print 'Saved file to', outputFile

	plt.close()
##########

def calculate_bubble_MC(z, ionized_cells, OutputDir, output_tag):
# Calculate the size of ionized/neutral bubbles by choosing an ionized/neutral cell and measuring the distance to a cell of the opposite phase.

## Input ##
# z is the redshift we are doing the MC walk at.
# ionized_cells is the array that contains the ionization state of a cell.

## Output ##
# The output file will be of the form '*OutputDir*/MC/*output_tag*_*z*.dat'

	def MC_walk(cells, indices, phase, N, outfile):
	# Function for the MC walk to calculate bubble size.
    
        ## Input ##
	# cells is the array that contains the ionization field (0 or 1).
	# indices is a length 3 array containing the x,y,z indices of the cells that are neutral (phase = 0) or ionized (phase = 1)
	# phase is 0 if we are starting with neutral or 1 if we are starting with ionized cell.
	# N is the number of random points we select.
    
        ## Output ##
        # A .dat file containing the results of the MC walk (number of steps until it encounters a phase transition).

		radii = []

		for j in xrange(0,int(N)):
                        
                        if (j%20000 == 0):
                            print "On walk number %d" %(j)
    
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
				phase_transition = cells[(walk_x + R*x) % GridSize, (walk_y + R*y) % GridSize, (walk_z + R*z) % GridSize]
				if (phase_transition > 0.8):
					phase_transition = 1
				else:
					phase_transition = 0	
				if (R >= BoxSize):
					phase_transition = (phase + 1) % 2

			radii.append(R)
			
		np.savetxt(outfile, radii, delimiter = ',')
		print "MC file saved as", outfile

	print
	print "Calculating bubble size using MC walk."	

	N = 1e5

	start_time = time.time()

	ionized_indices= np.where(ionized_cells > 0.8) # Calculate the array indices corresponding to neutral/ionized cells.
	unionized_indices = np.where(ionized_cells < 0.2)

	MCDir = OutputDir + 'MC/' 	
	if not os.path.exists(MCDir):
		os.makedirs(MCDir)
	outfile = MCDir + output_tag + '_z_%.3f.dat' %(ZZ[i])
	MC_walk(ionized_cells, ionized_indices, 1, N, outfile)

	print("MC took %s seconds." % (time.time() - start_time))


##########

def plot_bubble_MC(ZZ, model_tags, file_tags, OutputDir, output_tag):
# Plot the results from the MC walk.  Note that input_tag and labels are arrays that contains the input file tag and label for different realizations.

## Input ##
# ZZ is the redshift range we are plotting the results over.
# upperZ is the upper bound we are plotting to.
# input_tag is an array with the file tags for the input MC file for different realizations #### NOTE: This should be the same as *output_tag* in calculate_bubble_MC
# labels is an array that contains the graph labels for different realizations.

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*.*output_format*'

	print "Plotting results from MC walk."
	colors = ['r', 'b', 'g', 'm', 'c', 'k', 'y', 'b']
	linestyles = ['-', '--', '-.', ':']

	MCDir = OutputDir + 'MC/'
	
	ax = plt.subplot(111)
	for p in xrange(0, len(ZZ)):
		for q in xrange(0, len(file_tags)): 

			infile = MCDir + file_tags[q] + '_z_%.3f.dat' %(ZZ[p]) 
			fd = open(infile, 'rb')

			print("Plotting Bubble size of file %s") %(infile)
			R = np.loadtxt(fd)
			print("Maximum radius before scaling is %d cells.") %(max(R))
		
			R *= BoxSize/GridSize 
			print("Maximum radius after scaling is %.4f Mpc/h.") %(max(R))	

			binwidth = 1*BoxSize/GridSize

			R_low = 0.5*BoxSize/GridSize
            		if max(R) > BoxSize/2:
		                R_high = BoxSize/2 + 0.5*BoxSize/GridSize
			else:
		                R_high = max(R) + 0.5*BoxSize/GridSize
		        R_high = max(R) + 0.5*BoxSize/GridSize
			NB = (R_high - R_low) / binwidth
                
			(counts, binedges) = np.histogram(R, range=(R_low, R_high), bins=NB, density = True)

			# Set the x-axis values to be the centre of the bins
			xaxeshisto = binedges[:-1] + 0.5 * binwidth
            		counts *= xaxeshisto
            
                        print max(counts)
			if q == 0:
				label = "z = %.3f" %(ZZ[p])
			else:
				label = ""
			plt.plot(xaxeshisto, counts, ls = linestyles[q%len(linestyles)], label = label, color = colors[p%len(colors)])	

	for q in xrange(0, len(model_tags)):
		plt.plot(-1, -5, ls = linestyles[q], label = model_tags[q], color = 'k')

        plt.xscale('log', nonposy='clip')
	plt.axis([0, BoxSize, 0.0, 0.7])

	plt.xlabel(r'$R \: (h^{-1} \mathrm{Mpc})$')
	plt.ylabel(r'$R \: \: dP/dR$')

	leg = plt.legend(loc='upper right', numpoints=1,
			 labelspacing=0.1)
	leg.draw_frame(False)  # Don't want a box frame
	for t in leg.get_texts():  # Reduce the size of the text
		   t.set_fontsize('medium')

    	ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    	ax.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))



	outputFile = OutputDir + output_tag + output_format 
	plt.savefig(outputFile)  # Save the figure
	print 'Saved file to', outputFile
	plt.close()


##########

def Power_Spectrum(ionized_cells):
# This function takes the ionized_cells and calculates the associated power spectrum.
# It first determines the dimensions of the bins in k-space.  Then it calculates the 3D Power Spectrum before binning it into 1D.

## Input ##
# ionized_cells is the array that contains the ionization data.

## Output ##
# k_mid: The middle point for each of the k-space bins.
# PowSpec_1D: The 1-Dimensional Power Spectrum.
# error: The error in PowSpec_1D
    
    ## Calculating the bins in k-space. ##

    mid_point = GridSize/2 # Middle of the Grid.
    
    n1 = arange(0, GridSize)
    n1[1+mid_point:] -= GridSize
    n2 = n1**2
    n_mag = np.sqrt(add.outer(add.outer(n2, n2), n2)).ravel() # Finding the magnitude of k-space vectors.
    
    n_bins = (-1,) + tuple(arange(mid_point-1)+1.5) + (BoxSize*2,)
    
    k_bins = digitize(n_mag, n_bins) - 1 # Digitize basically takes an input array (N_Mag) and bins each value according to the binedges defined by N_bins.  The indices of the bin that each element of N_Mag belongs to is returned.
    
    k_multiplier = 2.0 * np.pi / BoxSize # Moving from real-space to k-space.
    
    k_min = (array(n_bins) * k_multiplier)[:-1] # Lower bound for each of the k-space bins.
    k_min[0] = 0
    
    k_max = (array(n_bins) * k_multiplier)[1:] # Upper bound for each of the k-space bins.
    k_max[-1] = mid_point * k_multiplier * sqrt(3.0) 

    k_volume = bincount(k_bins) * (k_multiplier)**3 # Volume in k-space (Number of cells of a specified type, divided by the volume of each cell).
    
    ##
    
    PowSpec_3D = (np.abs(np.fft.fftn(ionized_cells))**2).ravel() # Fourier transform into 3D Power Spectrum.

    PowSpec_1D = bincount(k_bins, weights=PowSpec_3D) # Takes the 3D Power Spectrum and bins the values into the k-space bins.
    PowSpec_1D /= k_volume # Normalize by the volume of the bins.  Necessary as spheres of larger radii occupy more volume by definition.
    
    # Calculate the error in the 1D Power Spectrum.
    v2 = bincount(k_bins, weights = (PowSpec_3D**2)) # Second Order
    v0 = bincount(k_bins) # Zeroth Order
    
    error = sqrt((v2*v0 - PowSpec_1D**2)/(v0-1)) / k_volume
    
    k_mid = 0.5 * (k_max + k_min) # Middle point for each of the bins in k-space.
#    PowSpec_1D *= k_mid**3 / 2 / (np.pi**2) / (BoxSize**3) 

    return k_mid, PowSpec_1D, error

##########

def plot_power(ZZ, k, Power, Error, labels, OutputDir, output_tag):

	ax = plt.subplot(111)
	for p in xrange(0, len(k)):
		print Power[p]
		plt.scatter(k[p], Power[p], color = colours[p], marker = markers[0], label = labels[p]) 

	for p in xrange(0, len(ZZ)):
		label = "z = %.3f" %(ZZ[p])
		plt.scatter(1e100, 1e100, color = 'k', marker = markers[0], label = label) 

	plt.xlabel(r'$k \: \left[\mathrm{Mpc}^{-1}h\right]$', size = label_size + extra_size)
	#plt.ylabel( r'$P\left(k\right) \: \: \left[\mathrm{Mpc}^3 h^{-3}\right]$', size = label_size + extra_size)
	plt.ylabel( r'$\bar{\delta T_b^2} \Delta_{21}^2 \left[\mathrm{mK}^2\right]$', size = label_size + extra_size)

        plt.xscale('log', nonposy='clip')
        plt.yscale('log', nonposy='clip')
        leg = plt.legend(loc='upper left', numpoints=1,
                         labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize('medium')

	plt.axis([0.5e-1, 1e1, np.amin(Power[:]), np.amax(Power[:])])	
    	#ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    	#ax.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))

	outputFile = OutputDir + output_tag + output_format 			
	plt.savefig(outputFile)  # Save the figure
	print 'Saved file to', outputFile
			
	plt.close()	





#########


def Reionization_Redshift(ionized_cells, density, redshift_array, z):

    print
    print "Checking to see what cells were reionized at redshift %.2f." %(z)

    tot_density = 0.0
    tot_density_std = 0.0
    count = 0.0

    for i in xrange(0, GridSize):
        for j in xrange(0, GridSize):
            for k in xrange(0, GridSize):
                if (ionized_cells[i,j,k] > 0.9 and redshift_array[i,j,k] == 0):
		    tot_density += density[i,j,k]
		    count += 1.0

    if (count > 0.0):
	density_z_mean = tot_density/count
    	for i in xrange(0, GridSize):
            for j in xrange(0, GridSize):
            	for k in xrange(0, GridSize):
                    if (ionized_cells[i,j,k] > 0.9 and redshift_array[i,j,k] == 0):
                    	redshift_array[i,j,k] = z
                        tot_density_std += (density[i,j,k] - density_z_mean)**2 	
		        
	tot_density_std = np.sqrt(tot_density_std / count)
    else:
	density_z_mean = 0.0
	tot_density_std = 0.0

    return redshift_array, density_z_mean, tot_density_std

##########


def photon_baryon(lowerZ, upperZ, ZZ, total_nion, Hubble_h, OB, Y, labels, OutputDir, output_tag): 


    Baryons = nb * (Mpc_m * 1e2)**3 # Number density of baryons (in cm^-3) to Mpc^-3 
    print Baryons
    Baryons *= (BoxSize/Hubble_h)**3 # Total number of baryons in our box
    print Baryons

    total_photons = np.zeros((upperZ-lowerZ)*len(total_nion))

    count = 0 


    for p in xrange(0, len(tota_nion)):
        total_photons[count] = (AllVars.Lookback_Time[0]-AllVars.Lookback_Time[lowerZ]) * total_nion[p][0] * AllVars.Sec_Per_Megayear*1e3
        count += 1
        for i in xrange(lowerZ+1, upperZ): 
#        total_photons[count] = total_photons[count-1] + (AllVars.Lookback_Time[i-1] - AllVars.Lookback_Time[i]) * total_nion[0][count] * AllVars.Sec_Per_Megayear*1e3
#        print AllVars.Lookback_Time[i-1] - AllVars.Lookback_Time[i]
            total_photons[count] = total_photons[count-1] + (AllVars.Lookback_Time[0] - AllVars.Lookback_Time[i])*total_nion[p][count] * AllVars.Sec_Per_Megayear*1e3
            count += 1

    plt.figure()

    colors = ['r', 'b', 'g', 'c', 'm']
    markers = ['x', 'o', '^', 's', 'D']

    for p in xrange(0, len(total_nion)):

        if (p == 0):
            lower_idx = 0
            upper_idx = len(total_nion[p])
        elif (p == 1):
            lower_idx = len(total_nion[0])
            upper_idx = len(total_nion[0]) + len(total_nion[1])

        plt.scatter(ZZ[0:len(total_nion[p])], np.log10(total_photons[p*len(total_nion[p]):len(total_nion[p])]/Baryons), color = colors[p], marker = markers[p], label = labels[p])

    plt.xlabel(r"$\mathrm{z}$")
    plt.ylabel(r"$\mathrm{log}_{10} \: \dot{N}_{\mathrm{HI}} \: \: [\mathrm{s}^{-1}]$")

    leg = plt.legend(loc='upper right', numpoints=1,labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')

    outputFile = OutputDir + output_tag + output_format 			
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
			
    plt.close()	
    

##########


def analytic_HII(total_nion, redshift, upperZ, snaplist, OutputDir, output_tag):


    Q_array = np.empty(upperZ)

    Q_array[0] = 0
    Q = 0

    for i in xrange(0, upperZ-1):

        print "=============================="
        print "REDSHIFT %.3f" %(ZZ[i])	
        print "=============================="


	## Using the prescription given by Lapi 2016 ## 
	
	hydrogen = n_HI(redshift[i], Hubble_h, OB, Y) * 1e-6 * Cubic_cm_Cubic_Mpc 
	#hydrogen = 2e-7 * (OB*(Hubble_h**2)/0.022) * Cubic_cm_Cubic_Mpc # Number density of Hydrogen in Mpc^-3 
	print "Hydrogen per Mpc^-3", hydrogen
	Nion = total_nion[i] * 1e1 / ((BoxSize/Hubble_h)**3) # Number of ionizing photons in s^-1 Mpc^-3
	print "Nion = %.4e s^-1 Mpc^-3" %(Nion)

	one_plus_z = (1 + redshift[i])/float(7)
	print "one_plus_z", one_plus_z

	t_rec = 3.2 * sec_per_Myr * 1e3 * (one_plus_z**(-3)) * (3**(-1)) # Recombination time in seconds. 
	print "t_rec", t_rec

	Q_dot = Nion/hydrogen - Q/t_rec

	print "Q_dot",Q_dot

	dt = (AllVars.Lookback_Time[snaplist[i]] - AllVars.Lookback_Time[snaplist[i+1]]) * sec_per_Myr * 1e3
	print len(AllVars.Lookback_Time)
	print "dt", dt

	Q = Q + Q_dot*dt

	Q_array[i+1] = Q
	print "Q", Q


    print 
    print "Plotting the analytic Q" 

    plt.figure()

    colors = ['r', 'b', 'g', 'c', 'm']
    markers = ['x', 'o', '^', 's', 'D']


    plt.scatter(redshift[0:upperZ], 1-Q_array, color = 'r', marker = 'x')

    plt.xlabel(r"$\mathrm{z}$")
    plt.ylabel(r"$Q_{HI}$")

    outputFile = OutputDir + output_tag + output_format 			
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
           
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
	print 'Saved file to', outputFile
			
	plt.close()	

##########

def save_redshifts(redshift, OutputDir, output_tag):


	outputFile = OutputDir + output_tag
	redshift.tofile(outputFile)

	print 'Saved redshift to file', outputFile 


##########

def plot_density_redshift(ZZ, density_mean, density_std, labels, OutputDir, output_tag): 

	ax = plt.subplot(111)

	for i in xrange(0, len(density_mean)):
		ax.errorbar(density_mean[i], ZZ, fmt = 'o', xerr = density_std[i], yerr = 0, label = labels[i], color = colours[i])

	ax.set_xlim([0.1, 1.6])
	ax.set_ylim([min(ZZ), max(ZZ)])

	ax.set_xlabel(r'$\Delta$')
	ax.set_ylabel(r'$z_\mathrm{reion}$')

    	leg = plt.legend(loc='upper left', numpoints=1,labelspacing=0.1)
    	leg.draw_frame(False)  # Don't want a box frame
   	for t in leg.get_texts():  # Reduce the size of the text
        	t.set_fontsize('medium')


	outputFile = OutputDir + output_tag + output_format 		
	plt.savefig(outputFile)  # Save the figure
	print 'Saved file to', outputFile
			
	plt.close()
	
##########

def plot_density_photons(ZZ, nion, density, count, labels, OutputDir, output_tag): 

	ax = plt.subplot(111)

	for i in xrange(0, len(nion)):
		w = np.where((nion[i].ravel() > 0.0))
		
		ax.scatter((density[i].ravel())[w], np.log10(((nion[i].ravel())[w])/((count[i].ravel())[w])), label = labels[i], color = colours[i], alpha = 0.5)
		#ax.scatter((density[i].ravel())[w], np.log10(((nion[i].ravel())[w])), label = labels[i], color = colours[i])

	ax.set_xlim([0.1, 10])
	ax.set_ylim([40, 60])

	ax.set_xlabel(r'$\Delta$')
	ax.set_ylabel(r'$\dot{N}_\Gamma$')

    	leg = plt.legend(loc='upper left', numpoints=1,labelspacing=0.1)
    	leg.draw_frame(False)  # Don't want a box frame
   	for t in leg.get_texts():  # Reduce the size of the text
        	t.set_fontsize('medium')


	outputFile = OutputDir + output_tag + output_format 		
	plt.savefig(outputFile)  # Save the figure
	print 'Saved file to', outputFile
			
	plt.close()

##########

def plot_twentyone(ZZ, ionized_cells, density, OutputDir, output_tag):

	ax = plt.subplot(111)

	T = T_naught(ZZ, OB, Hubble_h, OM)	
	
	brightness = T * density * ionized_cells
	im = ax.imshow(brightness[:,:,cut_slice:cut_slice+1].mean(axis = -1), interpolation='bilinear', origin='low', extent =[0,BoxSize,0,BoxSize], cmap = 'afmhot_r', norm = colors.LogNorm(vmin = 1, vmax = 150))

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
	print 'Saved file to', outputFile
			
	plt.close()	


##########

def twentyone_slice(ZZ, z_index, ionized_cells, density, brightness_slice):

	T = T_naught(ZZ, OB, Hubble_h, OM)
	
	brightness = T * density * ionized_cells

	for i in xrange(0, GridSize):
		brightness_slice[z_index, i] = brightness[i,i,i]	


	return brightness_slice

##########

def plot_twentyone_slice(ZZ, brightness_slice, OutputDir, output_tag):

	ax = plt.subplot(111)

	print brightness_slice
	print len(brightness_slice)
	print ZZ[-1]
	print ZZ[0]
	im = ax.imshow(brightness_slice.T, interpolation='bilinear', extent =[ZZ[-1], ZZ[0],0,BoxSize], cmap = 'afmhot_r', origin='low', aspect='auto')

	cbar = plt.colorbar(im, ax = ax)
	cbar.set_label(r'$\delta T_b$')
				    
	ax.set_xlabel(r'$z$')
	ax.set_ylabel(r'$\mathrm{x}  (h^{-1}Mpc)$')  

	ax.set_xlim([ZZ[-1], ZZ[0]]) 
	ax.set_ylim([0.0, GridSize])

	outputFile = OutputDir + output_tag + output_format 		
	plt.savefig(outputFile)  # Save the figure
	print 'Saved file to', outputFile
			
    	plt.tight_layout()

	plt.close()	
	

##########

if __name__ == '__main__':

    import os

    output_tags = ["512_noreion_reionmine", "512_noreion_reionmine"]
    model_tags = [r"$512 NoReion$", r"$512 Noreion_Reionmine$"]

    print "Model 1 is %s." %(output_tags[0])
    print "Model 2 is %s." %(output_tags[1])

    model = 'fescMH'

    OutputDir = "./ionization_plots/" + model + '/'
    if not os.path.exists(OutputDir):
        os.makedirs(OutputDir) 
    filepath_model1 = "/lustre/projects/p004_swin/jseiler/anne_output_january/XHII_noreion_512_fesc0.40_DRAGONSPhotHI"
    filepath_model2 = "/lustre/projects/p004_swin/jseiler/anne_output_january/XHII_noreion_512_fesc0.40_DRAGONSPhotHI"
    filepath_nion_model1 = "/lustre/projects/p004_swin/jseiler/SAGE_output/1024/grid/January_input/Galaxies_noreion_z5.000_fesc0.15_nionHI" 
    filepath_nion_model2 = "/lustre/projects/p004_swin/jseiler/SAGE_output/1024/grid/January_input/Galaxies_noreion_z5.000_alpha3.62beta-0.43_nionHI"
    filepath_density_model1 = "/lustre/projects/p004_swin/jseiler/densfield_output/512/256/snap_50.dens.dat"
    filepath_density_model2 = "/lustre/projects/p004_swin/jseiler/densfield_output/512/256/snap_50.dens.dat" 
    filepath_photofield_model1 = "/lustre/projects/p004_swin/jseiler/anne_output_january/PhotHI_noreion_fesc0.15_DRAGONSPhotHI"
    filepath_photofield_model2 = "/lustre/projects/p004_swin/jseiler/anne_output_january/PhotHI_noreion_512_fesc0.40_DRAGONSPhotHI"
    filepath_count_model1 = "/lustre/projects/p004_swin/jseiler/SAGE_output/1024/grid/January_input/Galaxies_noreion_z5.000_fesc0.15_count" 
    filepath_count_model2 = "/lustre/projects/p004_swin/jseiler/SAGE_output/1024/grid/January_input/Galaxies_noreion_z5.000_fesc0.15_count" 
    
    snaplist = np.arange(24, 78)

    ZZ = np.empty(len(snaplist))
    
    for i in xrange(0, len(snaplist)):
	ZZ[i] = AllVars.SnapZ[snaplist[i]]  
    z_index = 0 

    lowerZ = 0 
    upperZ = len(ZZ) 
    volume_frac_model1 = []
    volume_frac_model2 = []

    mass_frac_model1 = []
    mass_frac_model2 = []

    nion_total_model1 = []
    nion_total_model2 = []
    
    PowerSpectra_model1 = []
    k_model1 = []
    PowerSpectra_Error_model1 = []

    PowerSpectra_model2 = []
    k_model2 = []
    PowerSpectra_Error_model2 = []

    calculate_power = 0

    Redshift_Power = []

    Photo_Mean_model1 = []
    Photo_Mean_model2 = []

    Photo_Std_model1 = []
    Photo_Std_model2 = []

    redshift_array_model1 = np.zeros((GridSize, GridSize, GridSize))
    redshift_array_model2 = np.zeros((GridSize, GridSize, GridSize))

    density_z_mean_model1 = np.zeros((upperZ - lowerZ))
    density_z_mean_model2 = np.zeros((upperZ - lowerZ))

    density_z_std_model1 = np.zeros((upperZ - lowerZ))
    density_z_std_model2 = np.zeros((upperZ - lowerZ))

    brightness_slice = np.zeros((upperZ - lowerZ, GridSize), dtype = np.float32)

    MC_lower = 20 # Not Snapshot number, but rather the index of the snapshot in snaplist.  
    MC_upper = 40
    MC_snaps = [45, 50, 55, 60]
    MC_ZZ = np.empty(len(MC_snaps))
    for i in xrange(0, len(MC_snaps)):
	MC_ZZ[i] = AllVars.SnapZ[MC_snaps[i]]

    print "Plotting MC for Redshifts", MC_ZZ

    calculate_MC = 0 # 0 to NOT calculate the MC bubbles, 1 to calculate them. 

    for i in xrange(lowerZ, upperZ):

        print "=============================="
        print "REDSHIFT %.3f" %(ZZ[i])	
        print "=============================="

        if (i < 10):
            number_tag_anne = '_0%d' %(i)
            number_tag_mine = '_00%d' %(i)
        else:
            number_tag_anne = '_%d' %(i)
            number_tag_mine = '_0%d' %(i)

        fname_ionized_model1 = filepath_model1 + number_tag_anne 
        print
        print "Loading in data for %s Model from file %s" %(model_tags[0], fname_ionized_model1)
        fd = open(fname_ionized_model1, 'rb')
        ionized_cells_model1 = np.fromfile(fd, count = GridSize*GridSize*GridSize, dtype = np.float64)	
        ionized_cells_model1.shape = (GridSize, GridSize, GridSize)
        fd.close()


        fname_ionized_model2 = filepath_model2 + number_tag_anne 
        print "Loading in data for %s Model from file %s." %(model_tags[1], fname_ionized_model2)
        fd = open(fname_ionized_model2, 'rb')
        ionized_cells_model2 = np.fromfile(fd, count = GridSize*GridSize*GridSize, dtype = np.float64)	
        ionized_cells_model2.shape = (GridSize, GridSize, GridSize)
        fd.close()

        fname_nion_model1 = filepath_nion_model1 + number_tag_mine 
        print "Loading Nion data for %s Model from file %s." %(model_tags[0], fname_nion_model1)
        fd = open(fname_nion_model1, 'rb')
        nion_model1 = np.fromfile(fd, count = GridSize*GridSize*GridSize, dtype = np.float64)
        nion_model1.shape = (GridSize, GridSize, GridSize)
	fd.close()
	
        fname_nion_model2 = filepath_nion_model2 + number_tag_mine 
        print "Loading Nion data for %s Model from file %s." %(model_tags[1], fname_nion_model2)
        fd = open(fname_nion_model2, 'rb')
        nion_model2 = np.fromfile(fd, count = GridSize*GridSize*GridSize, dtype = np.float64)
        nion_model2.shape = (GridSize, GridSize, GridSize)
	fd.close()


        fname_density = filepath_density_model1 #+ number_tag_mine 
        print "Loading density data for %s Model from file %s." %(model_tags[0], fname_density)
        fd = open(fname_density, 'rb')
        density_model1 = np.fromfile(fd, count = GridSize*GridSize*GridSize, dtype = np.float64)
        density_model1.shape = (GridSize, GridSize, GridSize)
	fd.close()

        fname_density = filepath_density_model2 #+ number_tag_mine
        print "Loading density data for %s Model from file %s." %(model_tags[1], fname_density)
        fd = open(fname_density, 'rb')
        density_model2 = np.fromfile(fd, count = GridSize*GridSize*GridSize, dtype = np.float64)
        density_model2.shape = (GridSize, GridSize, GridSize)
	fd.close() 

        fname_photofield = filepath_photofield_model1 + number_tag_anne
        print "Loading photoionization data for %s Model from file %s." %(model_tags[0], fname_photofield)
        fd = open(fname_photofield, 'rb')
        photofield_model1 = np.fromfile(fd, count = GridSize*GridSize*GridSize, dtype = np.float64)
        photofield_model1.shape = (GridSize, GridSize, GridSize)
	fd.close()

        fname_photofield = filepath_photofield_model2 + number_tag_anne
        print "Loading photoionization data for %s Model from file %s." %(model_tags[1], fname_photofield)
        fd = open(fname_photofield, 'rb')
        photofield_model2 = np.fromfile(fd, count = GridSize*GridSize*GridSize, dtype = np.float64)
        photofield_model2.shape = (GridSize, GridSize, GridSize)
	fd.close()

        fname_count = filepath_count_model1 + number_tag_mine
        print "Loading counts for %s Model from file %s." %(model_tags[0], fname_count)
        fd = open(fname_count, 'rb')
        count_model1 = np.fromfile(fd, count = GridSize*GridSize*GridSize, dtype = np.int32)
        count_model1.shape = (GridSize, GridSize, GridSize)
	fd.close()

        fname_count = filepath_count_model2 + number_tag_mine
        print "Loading counts for %s Model from file %s." %(model_tags[1], fname_count)
        fd = open(fname_count, 'rb')
        count_model2 = np.fromfile(fd, count = GridSize*GridSize*GridSize, dtype = np.int32)
        count_model2.shape = (GridSize, GridSize, GridSize)
	fd.close()


        ##########################################################################
        ## Clean up/reduce the data; this will be fed into the other functions. ##
        ##########################################################################

        ionized_cells_clean_model1 = np.zeros((GridSize, GridSize, GridSize))
        ionized_cells_clean_model2 = np.zeros((GridSize, GridSize, GridSize))

        #w1 = np.where(ionized_cells_model1 > 0.5)
        #ionized_cells_clean_model1[w1] = 1
        #w2 = np.where(ionized_cells_model2 > 0.5)
        #ionized_cells_clean_model2[w2] = 1

        ##########################################################################

	if (i == 30 and calculate_power == 1):	
		## Model 1 ##
		Brightness_model1 = 27.0*(1.0-ionized_cells_model1)*density_model1*((1+ZZ[i])/10.0 * 0.15/(OM*Hubble_h*Hubble_h))**(1/2) * (OB*Hubble_h*Hubble_h / 0.023)

		Mean_Brightness_model1 = np.mean(Brightness_model1)
		delta_model1 = Brightness_model1/Mean_Brightness_model1 - 1	

	        (k_Bins , PowSpec, Error) = Power_Spectrum(delta_model1)
        	PowerSpectra_model1.append(PowSpec * k_Bins**3 / (2*np.pi*np.pi*BoxSize**3)) 
        	k_model1.append(k_Bins)
        	PowerSpectra_Error_model1.append(Error)

		## Model 2 ##
		Brightness_model2 = 27.0*(1.0-ionized_cells_model2)*density_model2*((1+ZZ[i])/10.0 * 0.15/(OM*Hubble_h*Hubble_h))**(1/2) * (OB*Hubble_h*Hubble_h / 0.023)
		Mean_Brightness_model2 = np.mean(Brightness_model2)
		delta_model2 = Brightness_model2/Mean_Brightness_model2 - 1	

        	(k_Bins , PowSpec, Error) = Power_Spectrum(delta_model2)
        	PowerSpectra_model2.append(PowSpec * k_Bins**3 / (2*np.pi*np.pi*BoxSize**3)) 
        	k_model2.append(k_Bins)
        	PowerSpectra_Error_model2.append(Error)

		Redshift_Power.append(ZZ[i])

		print "Mean_Brightness_model1", Mean_Brightness_model1
		print "Mean_Brightness_model2", Mean_Brightness_model2


       
      #  plot_single(ZZ[i], fname_ionized_model1, ionized_cells_model1, OutputDir, output_tags[0] + str(i)) 
       # plot_single(ZZ[i], fname_ionized_model2, ionized_cells_model2, OutputDir, output_tags[1] + str(i)) 
        #plot_comparison(ZZ[i], ionized_cells_clean_model1, ionized_cells_clean_model2, nion_model1, nion_model2, OutputDir, model_tags, "Comparison" + str(i))
        #plot_sources(ZZ[i], fname_nion_model1, nion_model1, OutputDir, model_tags[0] + "_nion")
        #plot_sources(ZZ[i], fname_nion_model2, nion_model2, OutputDir, model_tags[1] + "_nion")

        #difference_map(ZZ[i], ionized_cells_clean_model1, ionized_cells_clean_model2, OutputDir, model_tags, "Difference_Map" + str(i))

        #volume_frac_model1.append(calculate_volume_frac(ionized_cells_model1))
        #volume_frac_model2.append(calculate_volume_frac(ionized_cells_model2))
        #mass_frac_model1.append(calculate_mass_frac(ionized_cells_clean_model1, density_model1))
        #mass_frac_model2.append(calculate_mass_frac(ionized_cells_clean_model2, density_model2))

        nion_total_model1.append(calculate_total_nion(model_tags[0], nion_model1))
        nion_total_model2.append(calculate_total_nion(model_tags[1], nion_model2))
        
	#plot_nionfield(ZZ[i], nion_model1, OutputDir, "Nion" + output_tags[0] + '_' + str(i))
	#plot_nionfield(ZZ[i], nion_model1, OutputDir, "Nion" + output_tags[1] + '_' + str(i))

        plot_density(ZZ[i], density_model1, OutputDir, "Density_" + output_tags[0] + str(i))
        plot_density(ZZ[i], density_model2, OutputDir, "Density_" + output_tags[1] + str(i))
        
        #plot_density_numbers(ZZ[i], density_model1, OutputDir, "DensityNumbers" + str(i))
	
	if (i > MC_lower and i < MC_upper and calculate_MC == 1):
	        calculate_bubble_MC(ZZ[i], ionized_cells_model1, OutputDir, output_tags[0])
        	calculate_bubble_MC(ZZ[i], ionized_cells_model2, OutputDir, output_tags[1])

        #redshift_array_model1, density_z_mean_model1[i], density_z_std_model1[i] = Reionization_Redshift(ionized_cells_model1, density_model1, redshift_array_model1, ZZ[i])
        #redshift_array_model2, density_z_mean_model2[i], density_z_std_model2[i] = Reionization_Redshift(ionized_cells_model2, density_model2, redshift_array_model2, ZZ[i])

	#plot_photofield(ZZ[i], photofield_model1, OutputDir, "PhotHIField_" + output_tags[0] + '_' + str(i))
	#plot_photofield(ZZ[i], photofield_model2, OutputDir, "PhotHIField_" + output_tags[1] + '_' + str(i))
	#plot_density_numbers(ZZ[i], photofield_model1, OutputDir, "PhotHIField_Numbers_" + output_tags[0] + '_' + str(i))
	#plot_density_numbers(ZZ[i], photofield_model2, OutputDir, "PhotHIField_Numbers_" + output_tags[1] + '_' + str(i))

	#Photo_Mean_model1.append(np.mean(photofield_model1))
 	#Photo_Mean_model2.append(np.mean(photofield_model2))
	#Photo_Std_model1.append(np.std(photofield_model1))
	#Photo_Std_model2.append(np.std(photofield_model2))
        #plot_density_photons(ZZ[i], [nion_model1, nion_model2], [density_model1, density_model2], [count_model1, count_model2], output_tags, OutputDir, "density_photons_mean" + str(i))
	#plot_twentyone(ZZ[i], ionized_cells_model1, density_model1, OutputDir, "21cm" + output_tags[0] + str(i))

	#brightness_slice = twentyone_slice(ZZ[i], z_index, ionized_cells_model1, density_model1, brightness_slice)	
	
	z_index += 1

    #plot_redshifts(redshift_array_model1, ZZ, lowerZ, upperZ, OutputDir, "zreiond3_" + output_tags[0])
    #plot_redshifts(redshift_array_model2, ZZ, lowerZ, upperZ, OutputDir, "zreiond3_" + output_tags[1])

    #save_redshifts(redshift_array_model1, OutputDir, "ReionizationRedshift_" + output_tags[0])
    #save_redshifts(redshift_array_model2, OutputDir, "ReionizationRedshift_" + output_tags[1])
 
    #plot_power(Redshift_Power, [k_model1, k_model2], [PowerSpectra_model1, PowerSpectra_model2], [PowerSpectra_Error_model1, PowerSpectra_Error_model2], model_tags, OutputDir, "PowerSpectrum")

    #plot_bubble_MC(MC_ZZ, model_tags, output_tags, OutputDir, "BubbleSizes")
    #plot_global_frac(ZZ, [volume_frac_model1, volume_frac_model2], [mass_frac_model1, mass_frac_model2], model_tags, OutputDir, "GlobalFraction")
    plot_total_nion(ZZ, [nion_total_model1, nion_total_model2], model_tags, OutputDir, "nion_total")
    #photon_baryon(lowerZ, upperZ, ZZ, [nion_total_model1, nion_total_model2], Hubble_h, OB, Y, model_tags, OutputDir, "nion_total2")
    #analytic_HII(nion_total_model1, ZZ, upperZ, snaplist, OutputDir, "Q_Analytic")	
    #plot_photo_mean(ZZ, [Photo_Mean_model1, Photo_Mean_model2], [Photo_Std_model1, Photo_Std_model2], model_tags, OutputDir, "Mean_Photo")
    #plot_density_redshift(ZZ, [density_z_mean_model1, density_z_mean_model2], [density_z_std_model1, density_z_std_model2], model_tags, OutputDir, "density_redshift_" + output_tags[0])

    #plot_twentyone_slice(ZZ, brightness_slice, OutputDir, "21cm_Slice") 
