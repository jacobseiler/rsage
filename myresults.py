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

import matplotlib.ticker as mtick
import PlotScripts
import ReadScripts
import AllVars

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
 
colours = ['r', 'b', 'g', 'm', 'c', 'k', 'y']
markers = ['x', 'o', 'x', 'o', 'D', 's']
linestyles = ['-', '--', '-.', ':']

Output_Format = ".png"


def multiply(n):
	total = 1
	for i in xrange(0, len(n)):
		total *= n[i]
		print "%.4e" %(total)
	return total

def Calculate_Histogram(Data, Bin_Width, Weights):

# This calculates the counts and Bin_Edges for a given set of data.

## Input ##
# Data is an array containing the data to be binned.
# Bin_Width is the width of the bins.
# Weights is either 0 or 1.  0: Denotes that the histogram should be a frequency (count) histogram. 1: Denotes that the histogram should be a probability histogram.

## Output ##
# Counts: The count (either frequency or probability) in each bin.
# Bin_Edges: The location of the edges of the bins.
# Bin_Middle: The middle of the bins.

    mi = np.floor(min(Data)) - 10*Bin_Width
    ma = np.floor(max(Data)) + 10*Bin_Width
    NB = (ma - mi) / Bin_Width

    print "The maximum of the data is %.4e" %(max(Data))
    print "The minimum of the data is %.4e" %(min(Data))

    print "The total of the data is %.4e" %(sum(Data))

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

def Plot_Spectrum():

    SFR = [0.1, 0.3, 0.9, 2.7, 8.1, 24.3]
    Lambda_tot = []
    Total_tot = []
    Stellar_tot = []
    Labels = []
    for i in xrange(0, len(SFR)):

    	fname = "/Users/100921091/Desktop/sage-forked/sage-1/STARBURST/output/output_files/SFR_%.1f_Z_33.spectrum1" %(SFR[i])
	fd = open(fname, 'rb')
	t, Lambda, Total, Stellar, Nebular = np.loadtxt(fd, dtype = np.float32, skiprows=6, unpack=True)  
    
	Lambda_tot.append(Lambda[1222:2442])
	Total_tot.append(Total[1222:2442])
	Stellar_tot.append(Stellar[1222:2442])
	Labels.append("%.1f" %(SFR[i]))


    print Stellar_tot[0]
    print Stellar_tot[1]
    print Stellar_tot[1]/Stellar_tot[0]

    PlotScripts.Plot_Line(Lambda_tot, Stellar_tot, [0], [1,0, ["Galaxies"]], Labels, [80, 2e6, 28, 42], [1, 0], ['Lambda', 'L'], 1, 'Stellar2', '.png') 

def Plot_Ngamma_Ratio():

    SFR = np.arange(0.01, 1, 0.01)
    Ngamma = []
    delta_SFR = np.empty(len(SFR)-1)
    delta_Ngamma = np.empty(len(SFR)-1)

    for i in xrange(0, len(SFR)):
        fname = "/Users/100921091/Desktop/sage-forked/sage-1/STARBURST/output/output_files/SFR_%.2f_Z_35.quanta1" %(SFR[i])
        fd = open(fname, 'rb')
        t, N_HI, L_HI, N_HeI, L_HeI, N_HeII, L_HeII, L = np.loadtxt(fd, dtype = np.float32, skiprows=7, unpack=True)	

        Ngamma.append(N_HI[-1])

    for i in xrange(1, len(SFR)):
        delta_SFR[i-1] = SFR[i]/float(SFR[i-1])
        delta_Ngamma[i-1] = Ngamma[i]/float(Ngamma[i-1])

    SFR = np.log10(SFR)

    PlotScripts.Fit_Linear([SFR], [Ngamma], ["Galaxies"], [min(SFR), max(SFR), min(Ngamma), max(Ngamma)], [r'log(SFR)', r'log(N$\gamma$)'], 1, 'SFR_Ngamma_Z35', '.png') 
#    PlotScripts.Plot_Scatter([SFR], [Ngamma], ["Galaxies"], [min(SFR), max(SFR), min(Ngamma), max(Ngamma)], [1,0], ['Delta SFR', 'Delta Ngamma'], 1, 'deltaSFR_deltaNgamma', '.png') 


def Photon_Totals(Simulation, Redshift, Photons, Mysim_Len): 


    title = ["MySim 512 Halo", "MySim 512 Galaxy", "Millennium 1024 Halo", "Millennium 1024 Galaxy"]

    if Simulation == 0:
        OutputFile = 'TotalPhotons_Millennium'
    elif Simulation == 1:
        OutputFile = 'TotalPhotons_MySim'
    elif Simulation == 2:
        OutputFile = 'TotalPhotons_Both'

    print Redshift
    print Photons

    PlotScripts.Plot_Scatter(Redshift, Photons, title, [4.5, 15.5, 53, 58.3], [0,0], ['z', r'Log Ionizing Photons [s$^{-1}$]'], 1, 2, OutputFile, '.png')  

def StellarMassFunction(Simulation, Redshift, Mass, HaloPartStellarMass, MySim_Len):

    title = []
    Normalization = []

    ## Plot Parameters ##
    binwidth = 0.1
    Observations = 1
    Output_Tag = "512_vs_1024"
    Frequency = 0 # 0 for a frequency (count) histogram, 1 for a probaility histogram.
    ##

    ## Normalization and Titles for MySim ##
    for i in xrange(0, MySim_Len): 
        tmp = '512: z = %.2f' %(Redshift[i])
        title.append(tmp)

        AllVars.Set_Params_Mysim()
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        Normalization.append(norm)


    ## Normalization and Titles for MySim2 ## 
    for i in xrange(0, MySim_Len): 
        tmp = '1024: z = %.2f' %(Redshift[i])
        title.append(tmp)

        AllVars.Set_Params_Mysim()
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        Normalization.append(norm)

### Plotting ###

    plt.figure()  
    ax = plt.subplot(111)  

### Plots the Histograms ###

    for i in xrange(0, len(Mass)):

        (counts, Bin_Edges, Bin_Middle) = Calculate_Histogram(Mass[i], binwidth, Frequency)
        if (i < MySim_Len):
            ls = '-'
        else:
            ls = '--'
        plt.plot(Bin_Middle, counts / Normalization[i], colours[i], linestyle = ls, label = title[i])

##

### Draws a vertical line to denote lower bounds for what is an 'acceptable' Stellar Mass ### 

    plt.axvline(x = np.log10(HaloPartStellarMass), ymin = 0, ymax = 10, linestyle = '-.')

## 

    plt.yscale('log', nonposy='clip')

    plt.axis([6, 11.5, 1e-6, 1e-1])

    ax.set_xlabel(r'$\log_{10}\ M_{\mathrm{*}}\ (M_{\odot})$')
    ax.set_ylabel(r'$\Phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$')
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

### If we want to put observations on the figure ###

    if Observations == 1:
    	Gonzalez_z6 = np.array([[7.77, -2.0956, -1.8596, -2.3539],
                                [8.27, -2.1742, -1.9494, -2.4101],
                                [8.77, -2.5674, -2.3876, -2.7921],
                                [9.27, -2.8483, -2.6573, -3.0843],
                                [9.77, -3.5787, -3.3764, -3.8258],
                                [10.27, -4.3202, -4.0281, -4.5674]], dtype = np.float32)

                                #plt.errorbar(Gonzalez_z6[:,0], 10**Gonzalez_z6[:,1], yerr= (10**Gonzalez_z6[:,3], 10**Gonzalez_z6[:,2]) , alpha=0.75, lw=1.0, marker='o', ls='none', label = 'Gonzalez 2011, z = 6', color = 'cyan')


        Gonzalez_z7 = np.array([[7.75, -2.1828, -1.7463, -2.5858],
                                [8.26, -2.25, -1.8694, -2.2631],
                                [8.77, -2.7425, -2.3731, -3.1231],
                                [9.27, -3.0672, -2.6753, -3.4142],
                                [9.76, -3.8731, -3.4831, -4.2537]], dtype = np.float32)

#plt.errorbar(Gonzalez_z7[:,0], 10**Gonzalez_z7[:,1], yerr= (10**Gonzalez_z7[:,3], 10**Gonzalez_z7[:,2]) , alpha=0.75, lw=1.0, marker='o', ls='none'    , label = 'Gonzalez 2011, z = 7', color = 'magenta')

        Song_z6 = np.array([[7.25, -1.47, -1.47 + 0.35, -1.47 - 0.23],
                            [7.75, -1.81, -1.81 + 0.23, -1.81 - 0.28],
                            [8.25, -2.26, -2.26 + 0.21, -2.26 - 0.16],
                            [8.75, -2.65, -2.65 + 0.15, -2.65 - 0.15],
                            [9.25, -3.14, -3.14 + 0.12, -3.14 - 0.11],
                            [9.75, -3.69, -3.69 + 0.12, -3.69 - 0.13],
                            [10.25, -4.27, -4.27 + 0.38, -4.27 - 0.86]], dtype = np.float32)

        plt.errorbar(Song_z6[:,0], 10**Song_z6[:,1], yerr= (10**Song_z6[:,1] - 10**Song_z6[:,3], 10**Song_z6[:,2] - 10**Song_z6[:,1]), xerr = 0.25, alpha = 1.0, lw=2.0, marker='o', ls='none', label = 'Song 2015, z = 6', color = '#bd0026')

        Song_z7 = np.array([[7.25, -1.63, -1.63 + 0.54, -1.63 - 0.54],
                            [7.75, -2.07, -2.07 + 0.45, -2.07 - 0.41],
                            [8.25, -2.49, -2.49 + 0.38, -2.49 - 0.32],
                            [8.75, -2.96, -2.96 + 0.32, -2.96 - 0.30],
                            [9.25, -3.47, -3.47 + 0.32, -3.47 - 0.35],
                            [9.75, -4.11, -4.11 + 0.41, -4.11 - 0.57],
                            [10.25, -4.61, -4.61 + 0.72, -4.61 - 0.82],
                            [10.75, -5.24, -5.24 + 0.90, -5.25 - 0.57]], dtype = np.float32)

        plt.errorbar(Song_z7[:,0], 10**Song_z7[:,1], yerr= (10**Song_z7[:,1] - 10**Song_z7[:,3], 10**Song_z7[:,2] - 10**Song_z7[:,1]), xerr = 0.25, alpha=0.75, lw=1.0, marker='o', ls='none', label = 'Song 2015, z = 7', color = 'blue')

        Song_z8 = np.array([[7.25, -1.73, -1.73 + 1.01, -1.73 - 0.84],
                            [7.75, -2.28, -2.28 + 0.84, -2.28 - 0.64],
                            [8.25, -2.88, -2.88 + 0.75, -2.88 - 0.57],
                            [8.75, -3.45, -3.45 + 0.57, -3.45 - 0.60],
                            [9.25, -4.21, -4.21 + 0.63, -4.21 - 0.78],
                            [9.75, -5.31, -5.31 + 1.01, -5.31 - 1.64]], dtype = np.float32)


        plt.errorbar(Song_z8[:,0], 10**Song_z8[:,1], yerr= (10**Song_z8[:,1] - 10**Song_z8[:,3], 10**Song_z8[:,2] - 10**Song_z8[:,1]), xerr = 0.25, alpha=0.75, lw=1.0, marker='o', ls='none', label = 'Song 2015, z = 8', color = 'green')

    leg = plt.legend(loc=1, numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')

    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()
 

#####

def CentralGalaxy_Comparison(Simulation, Redshift, Mass, Photons):

    min_mass_H = 9
    max_mass_H = 12
    min_mass_G = 9
    max_mass_G = 12
    min_Photons = 49
    max_Photons = 53

    if Simulation == 0:
        OutputFile = 'CentralStellar_Photons_Fiducial_Millennium'
    elif Simulation == 1:
        OutputFile = 'CentralStellar_Photons_Fiducial_MySim_NoReion'

    PlotScripts.Plot_Scatter_SixPanel(Mass, Photons, 0, 1, ['Halos', 'Galaxies'], [min_mass_H, max_mass_H, min_Photons, max_Photons], [r'Log Halo Mass [$M_{\odot}$]', r'Log Ionizing Photons [s$^{-1}$]'], 2, Redshift, OutputFile, '.png')
    #PlotScripts.Plot_HeatMap_SixPanel(Mass, Photons, 0, 1, ['Halos', 'Galaxies'], [min_mass_H, max_mass_H, min_Photons, max_Photons], [r'Log Halo Mass [$M_{\odot}$]', r'Log Ionizing Photons [s$^{-1}$]'], 2, Redshift, OutputFile, '.png')

def CentralGalaxy_Comparison_Difference(Simulation,Redshift, Mass, Photons):

    min_mass_H = 9
    max_mass_H = 12
    min_mass_G = 9
    max_mass_G = 12

    if Simulation == 0:
        OutputFile = 'CentralStellar_Photons_Fiducial_Millennium_Difference'
    elif Simulation == 1:
        OutputFile = 'CentralStellar_Photons_Fiducial_MySim_Difference'

    PlotScripts.Plot_Scatter_SixPanel(Mass, Photons, 1, 1, ['Halos', 'Galaxies'], [min_mass_H, max_mass_H, -2, 1], [r'Log Halo Mass [$M_{\odot}$]', r'Log (Halo/Galaxy Ionizing Photons [s$^{-1}$]'], 2, Redshift, OutputFile, '.png')

##

def CentralGalaxy_Projection(Simulation, Redshift, Mass, Photons):

    if Simulation == 0:
        OutputFile = 'CentralStellar_Photons_Fiducial_Millennium_HistProjection'
    elif Simulation == 1:
        OutputFile = 'CentralStellar_Photons_Fiducial_MySim_HistProjection'

    title = []

    for i in xrange(0, len(Redshift)):
        tmp = 'z = %.2f: 512' %(Redshift[i])
        title.append(tmp)

    for i in xrange(0, len(Redshift)):
        tmp = 'z = %.2f: 1024' %(Redshift[i])
        title.append(tmp)

    PlotScripts.Plot_Scatter_TwoHists(Mass, Photons, title, [9, 12, 51, 54], [0, 0.30, 0, 0.30], 0.2, 1, 1, [r'Log Halo Mass [$M_{\odot}$]', r'Log Ionizing Photons [s$^{-1}$]'], 0, OutputFile, '.png') 

##

def Metallicity(Simulation, Redshift, Mass, Z):

    title = []
    for i in xrange(0, len(Redshift)):
        tmp = 'z = %.2f' %(Redshift[i])
        title.append(tmp)

    PlotScripts.Plot_Scatter(Mass, Z, title, [6, 11, 7.0, 9.5], [0, 0], [r'$\log_{10}\ M_{\mathrm{*}}\ (M_{\odot})$', r'$12\ +\ \log_{10}[\mathrm{O/H}]$'], 4, 1, "Z", ".png")
    PlotScripts.Plot_Scatter_TwoHists(Mass, Z, title, [6, 11, 7.0, 9.5], [0, 0.1, 0, 0.15], 0.1, 1, 1, [r'$\log_{10}\ M_{\mathrm{*}}\ (M_{\odot})$', r'$12\ +\ \log_{10}[\mathrm{O/H}]$'], 1, 'Z_Hist', '.png')

##

def HaloMassFunction(Simulation, Redshift, Mass, Mysim_Len): 

    title = []
    normalization = []
    binwidth = 0.1

    
    for i in xrange(0, Mysim_Len): 
        tmp = 'Mysim 512: z = %.2f' %(Redshift[i])
        title.append(tmp)

        AllVars.Set_Params_Mysim()
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        normalization.append(norm)
 
    for i in xrange(0, Mysim_Len): 
        tmp = 'Mysim 1024: z = %.2f' %(Redshift[i])
        title.append(tmp)

        AllVars.Set_Params_Mysim()
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        normalization.append(norm)

    '''    
    for i in xrange(Mysim_Len, len(Redshift)):

        tmp = 'Millennium: z = %.2f' %(Redshift[i])
        title.append(tmp)

        AllVars.Set_Params_MiniMill()
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        normalization.append(norm)
    ''' 
    if Simulation == 0:
        OutputFile = 'HaloMassFunction_Millennium'
    elif Simulation == 1:
        OutputFile = 'HaloMassFunction_Comparison'
    elif Simulation ==2:   
        OutputFile = 'HaloMassFunction_Comparison'

    PlotScripts.Plot_Histogram(Mass, title, [7, 13.5, 1e-6, 1e-1], binwidth, normalization, 0, [0, 1], [r'$\log_{10}\ M_{\mathrm{Vir}}\ (M_{\odot})$', r'$\Phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$'], 1, 0, OutputFile, '.png')

##

def Central_Galaxy_Projection(Simulation, Redshift, Mass, Photons):

    if Simulation == 0:
        OutputFile = 'CentralStellar_Photons_Fiducial_Millennium_HistProjection'
    elif Simulation == 1:
        OutputFile = 'CentralStellar_Photons_Fiducial_MySim_HistProjection'



    PlotScripts.Plot_Scatter_TwoHists(Mass, Photons, title, [9, 12, 52, 55], [0, 0.25, 0, 0.25], 0.2, 1, 1, [r'Log Halo Mass [$M_{\odot}$]', r'Log Ionizing Photons [s$^{-1}$]'], 0, OutputFile, '.png') 



##

def FilteringMass(Simulation, SnapListZ, Gnedin, Sobacchi, Gnedinv2, Mvir, MySim_Len):

    title = ["Gnedin", "Sobacchi", "Gnedin Masstouse", "Mvir"]
    Normalization = []
 

    ## Plot Parameters ##
    binwidth = 0.1
    Observations = 0
    Output_Tag = "Filtering_Mass_to_Use"

    ##

    Mean_Array = []
    Std_Array = []

    Gnedin_Mean = np.zeros((len(SnapListZ)))
    Gnedin_Std = np.zeros((len(SnapListZ)))

    Gnedin_Masstouse_Mean = np.zeros((len(SnapListZ)))
    Gnedin_Masstouse_Std = np.zeros((len(SnapListZ)))


    Sobacchi_Mean = np.zeros((len(SnapListZ)))
    Sobacchi_Std = np.zeros((len(SnapListZ)))

    Mvir_Mean = np.zeros((len(SnapListZ)))
    Mvir_Std = np.zeros((len(SnapListZ)))


    print "Sobacchi", Sobacchi
    print "SNAPSHOTS", SnapListZ
    for i in xrange(0, len(SnapListZ)):
	Gnedin_Mean[i] = np.mean(Gnedin[i])
	Gnedin_Std[i] = np.std(Gnedin[i])

	w = np.where(Sobacchi[i] > 1.0)[0]

	Sobacchi_Mean[i] = np.mean(Sobacchi[i][w])
	Sobacchi_Std[i] = np.std(Sobacchi[i][w])

	Mvir_Mean[i] = np.mean(10**Mvir[i][w])
	Mvir_Std[i] = np.std(10**Mvir[i][w])
 
	Gnedin_Masstouse_Mean[i] = np.mean(Gnedinv2[i])
	Gnedin_Masstouse_Std[i] = np.std(Gnedinv2[i])

    fig = plt.figure()
    ax = plt.subplot(111)
          
    ax.plot(SnapListZ, np.log10(Gnedin_Mean), color = 'r', linestyle = '--', label = title[0])
#    ax.fill_between(SnapListZ, np.log10(Gnedin_Mean)-0.434*Gnedin_Std/Gnedin_Mean, np.log10(Gnedin_Mean)+0.434*Gnedin_Std[i]/Gnedin_Mean, alpha = 0.5, color = 'r')	
    print "Gnedin_Mean", Gnedin_Mean
#    print "Gnedin_Std", Gnedin_Std

    print "Sobacchi_Mean", np.log10(Sobacchi_Mean)
    print "Sobacchi_Std", np.log10(Sobacchi_Std)

    ax.plot(SnapListZ, np.log10(Sobacchi_Mean), color = 'b', linestyle = '-', label = title[1])
    ax.fill_between(SnapListZ, np.log10(Sobacchi_Mean)-0.434*Sobacchi_Std/Sobacchi_Mean, np.log10(Sobacchi_Mean)+0.434*Sobacchi_Std/Sobacchi_Mean, alpha = 0.5, color = 'b')	

    print "Mvir_Mean", np.log10(Mvir_Mean)
    print "Mvir_Std", np.log10(Mvir_Std)
    ax.plot(SnapListZ, np.log10(Mvir_Mean), color = 'k', linestyle = '-', label = r"$M_{\mathrm{vir}}$")
    ax.fill_between(SnapListZ, np.log10(Mvir_Mean)-0.434*Mvir_Std/Mvir_Mean, np.log10(Mvir_Mean)+0.434*Mvir_Std/Mvir_Mean, alpha = 0.5, color = 'k')

    ax.plot(SnapListZ, np.log10(Gnedin_Masstouse_Mean), color = 'b', linestyle = '-', label = title[2])

    ax.set_xlabel(r"z", fontsize = label_size + extra_size)
    ax.set_ylabel(r"$\mathrm{log}_{10} \: M_\odot$", fontsize = label_size + extra_size)

    ax.set_xlim([5.5,13])
    ax.set_ylim([6,13])

    ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    ax.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))

    leg = plt.legend(loc= 1, numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')


    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()


##

def StarFormationRate(Simulation, SnapListZ, SFR1, SFR2):

    title = [r"$\mathrm{Noreion Reionmine}$", r"$\mathrm{Noreion}$"]
 
    ## Plot Parameters ##
    binwidth = 0.1
    Observations = 0
    Output_Tag = "SFR"

    ##

    SFR1_Mean = np.zeros((len(SnapListZ))) 
    SFR1_Std = np.zeros((len(SnapListZ)))

    SFR2_Mean = np.zeros((len(SnapListZ)))
    SFR2_Std = np.zeros((len(SnapListZ)))

    print "SFR1", SFR1
    print "SFR2", SFR2

    for i in xrange(0, len(SnapListZ)):

	SFR1_Mean[i] = np.mean(10**SFR1[i])
	SFR1_Std[i] = np.std(10**SFR1[i])

	SFR2_Mean[i] = np.mean(10**SFR2[i])
	SFR2_Std[i] = np.std(10**SFR2[i])

 
    fig = plt.figure()
    ax = plt.subplot(111)
          
    ax.plot(SnapListZ, np.log10(SFR1_Mean), color = 'r', linestyle = '-', label = title[0])
    ax.fill_between(SnapListZ, np.log10(SFR1_Mean)-0.434*SFR1_Std/SFR1_Mean, np.log10(SFR1_Mean)+0.434*SFR1_Std[i]/SFR1_Mean, alpha = 0.5, color = 'r')	

    print "SFR1_Mean", np.log10(SFR1_Mean)
    print "SFR1_Std", np.log10(SFR1_Std)

    ax.plot(SnapListZ, np.log10(SFR2_Mean), color = 'b', linestyle = '-', label = title[1])
    ax.fill_between(SnapListZ, np.log10(SFR2_Mean)-0.434*SFR2_Std/SFR2_Mean, np.log10(SFR2_Mean)+0.434*SFR2_Std[i]/SFR2_Mean, alpha = 0.5, color = 'b')	

    print "SFR2_Mean", np.log10(SFR2_Mean)
    print "SFR2_Std", np.log10(SFR2_Std)

    ax.set_xlabel(r"z", fontsize = label_size + extra_size)
    ax.set_ylabel(r"$\mathrm{SFR} \: [M_\odot\mathrm{yr}^{-1}]$", fontsize = label_size + extra_size)

    ax.set_xlim([5.5,13])
    ax.set_ylim([-3,2.5])

    ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    ax.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))

    leg = plt.legend(loc= 1, numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')


    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()


##

def PhotonsStellarMass(Simulation, SnapListZ, Mass, Photons):

	Output_Tag = "PhotonsStellarMass"

	bound1 = 7 
	bound2 = 8
	bound3 = 9
	bound4 = 10
	bound5 = 11  

	bounds = [7, 8, 9, 10, 100]
	labels = [r"$7 \leq M_* \le 8$", r"$8 \leq M_* \le 9$", r"$9 \leq M_* \le 10$", r"$M_* \geq 10$"]

	photons_mean = []
	photons_std = []

	for i in xrange(0, len(bounds)-1):
		for j in xrange(0, len(SnapListZ)):
			w = np.where((Mass[j] >= bounds[i]) & (Mass[j] < bounds[i+1]))[0]
			print "Redshift %.4f, bounds %d-%d" %(SnapListZ[j], bounds[i], bounds[i+1])
			print "len(w)", len(w)
			print "W", w
			print "Photons[j][w]", Photons[j][w]

			if (len(w) > 0):
				photons_mean.append(Sum_Log(Photons[j][w]))
#				photons_std.append(Std_Log(Photons[j][w], Sum_Log(Photons[j][w])/len(w)))
			else:
				photons_mean.append(nan)
				photons_std.append(nan)

			print "photons_mean", photons_mean[i*len(SnapListZ) + j]

	## Photons are now in an array of [Bin1, Bin1, Bin1, Bin2, Bin2, Bin2, Bin3, Bin3, Bin3]
	##                                SNAP1,SNAP2,SNAP3,SNAP1,SNAP2,SNAP3,SNAP1,SNAP2,SNAP3

	ax = plt.subplot(111)

	for i in xrange(0, len(bounds)-1):
		print "For bounds ", bounds[i], "-", bounds[i+1], "photons_mean = ", photons_mean[i*len(SnapListZ):(i+1)*len(SnapListZ)]
		low_index = i*len(SnapListZ)
		high_index = (i+1)*len(SnapListZ)

		ax.plot(SnapListZ, np.log10(photons_mean[low_index:high_index]), color = colours[i], label = labels[i]) 
#    		ax.fill_between(SnapListZ, np.log10(photons_mean[low_index:high_index])-0.434*photons_std[low_index:high_index]/photons_mean[low_index:high_index], np.log10(photons_mean[low_index:high_index])+0.434*photons_std[low_index:high_index]/photons_mean[low_index:high_index], alpha = 0.5, color = colours[i])	


	ax.set_xlabel(r"z", fontsize = label_size + extra_size)
	ax.set_ylabel(r'Total $\dot{N}_{\gamma,\mathrm{HI}} \: [\mathrm{s}^{-1}]$', fontsize = label_size + extra_size)
    
	ax.set_xlim([5.5,13])
	ax.set_ylim([50, 60])

	
        leg = plt.legend(loc= 1, numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
        	t.set_fontsize('medium')


        outputFile = './' + Output_Tag + Output_Format
        plt.savefig(outputFile)  # Save the figure
        print 'Saved file to', outputFile
        plt.close()


##

def LymanAlphaLF(Simulation, Redshift, LymanAlpha, MySim_Len):

    title = []
    Normalization = []

    print LymanAlpha
    ## Plot Parameters ##
    binwidth = 0.1
    Observations = 0
    Output_Tag = "Test_LymanAlpha"
    Frequency = 0 # 0 for a frequency (count) histogram, 1 for a probaility histogram.
    ##

    ## Normalization and Titles for MySim ##
    for i in xrange(0, MySim_Len): 
        tmp = '512: z = %.2f' %(Redshift[i])
        title.append(tmp)

        AllVars.Set_Params_Mysim()
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        Normalization.append(norm)


    ## Normalization and Titles for MySim2 ## 
    for i in xrange(0, MySim_Len): 
        tmp = '1024: z = %.2f' %(Redshift[i])
        title.append(tmp)

        AllVars.Set_Params_Mysim()
        norm = pow(AllVars.BoxSize,3) / pow(AllVars.Hubble_h, 3) * binwidth 
    
        Normalization.append(norm)

### Plotting ###

    plt.figure()  
    ax = plt.subplot(111)  

### Plots the Histograms ###

    for i in xrange(0, len(LymanAlpha)):

        (counts, Bin_Edges, Bin_Middle) = Calculate_Histogram(LymanAlpha[i], binwidth, Frequency)
        if (i < MySim_Len):
            ls = '-'
        else:
            ls = '--'
        plt.plot(Bin_Middle, counts / Normalization[i], colours[i], linestyle = ls, label = title[i])

##


## 

    plt.yscale('log', nonposy='clip')

#    plt.axis([6, 11.5, 1e-6, 1e-1])

    ax.set_xlabel(r'$\log_{10}\ L_{\mathrm{Ly}\alpha}\ [\mathrm{erg} \: \mathrm{s}^{-1}]$')
    ax.set_ylabel(r'$\Phi\ (\mathrm{Mpc}^{-3}\ \mathrm{dex}^{-1})$')
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

### If we want to put observations on the figure ###

    if Observations == 1:
    	Konno_2014 = np.array([[42.5, 2.5e-4, 5e-4, 1e-4],
                               	[42.65, 9.5e-5, 5e-5, 2e-4],  
                               	[42.85, 2e-5, 7.5e-5, 3e-6]], dtype = np.float32) 

        plt.errorbar(Konno_2014[:,0], Konno_2014[:,1], yerr= (Konno_2014[:,1] - Konno_2014[:,3], 10**Song_z6[:,2] - 10**Song_z6[:,1]), xerr = 0.25, alpha = 1.0, lw=2.0, marker='o', ls='none', label = 'Song 2015, z = 6', color = '#bd0026')
    
    leg = plt.legend(loc=1, numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')

    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()
 
##

def PhotonsVsStellarMass_OLD(Simulation, SnapListZ, Mass, Photons):
   
    print "Plotting a Six panel Ngamma vs Stellar Mass plot."
 
    nrows = 2
    ncols = 3

    binwidth = 0.1
    low_mass = 7
    high_mass = 12
   
    Output_Tag = "PhotonsVsStellarMass"

    Limits = [low_mass, high_mass, 53, 57]

    bins = np.arange(low_mass,high_mass, binwidth)
    bins_mid = bins + binwidth/2.0
    bins_mid = bins_mid[:-1] # len(bins_mid) should be 1 less than len(bins) as the last bin doesn't have a midpoint.


    fig = plt.figure()
    for i in xrange(1, nrows*ncols + 1):
	print "On panel %d" %(i)
        ax = fig.add_subplot(nrows, ncols, i)
    
	Photons_Sum = []

	for j in xrange(0, len(bins)-1):
		w = np.where((Mass[i-1] >= bins[j]) & (Mass[i-1] < bins[j+1]))[0]
		print "w", w	
		if (len(w) != 0):
			Photons_Sum.append(Sum_Log(Photons[i-1][w]))
		else:
			Photons_Sum.append(nan)
		print "Photons_Sum[-1]", Photons_Sum[-1]
	ax.plot(bins_mid, np.log10(Photons_Sum))
		
        ax.set_xticks(np.arange(Limits[0], Limits[1] + 1.0, 1))
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f'))
        ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
        
        ax.set_yticks(np.arange(Limits[2], Limits[3] + 0.5, 0.5))
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
        ax.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))

        if i < 4:
            ax.tick_params(axis='both', labelbottom = 'off')
        else:
            ax.xaxis.tick_bottom()
        if i > 3:
            ax.xaxis.get_major_ticks()[-1].set_visible(False)

        if i != 1 and i != 4:
            ax.tick_params(axis = 'y', labelleft = 'off')
        if i == 1:
            ax.yaxis.get_major_ticks()[0].set_visible(True)
        if i == 4:
            ax.yaxis.get_major_ticks()[-1].set_visible(False)

        if i == 6:
            ax.xaxis.get_major_ticks()[-1].set_visible(True)


        ax.set_xlim(Limits[0:2])
        ax.set_ylim(Limits[2:4])
        
        label = "z = %.2f" %(SnapListZ[i-1])
        ax.text(11.25, 56.0, label, fontsize = label_size)


    plt.tight_layout()
    fig.text(0.5, 0.01, "Mass", ha = 'center')
    fig.text(0.001, 0.5, "Photons", va = 'center', rotation = 'vertical')
    fig.subplots_adjust(hspace = 0.0, wspace = 0.0)
    
    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

##

def PhotonsVsStellarMass(Simulation, SnapListZ, Mass, Photons):
   
    print "Plotting histogram for Ngamma vs M*" 
 
    binwidth = 0.1
    low_mass = 7
    high_mass = 12
   
    Output_Tag = "PhotonsVsStellarMass"

    Limits = [low_mass, high_mass, 53, 57]

    bins = np.arange(low_mass,high_mass, binwidth)
    bins_mid = bins + binwidth/2.0
    bins_mid = bins_mid[:-1] # len(bins_mid) should be 1 less than len(bins) as the last bin doesn't have a midpoint.


    ax = plt.subplot(111) 
    for i in xrange(0, len(SnapListZ)):
    
	Photons_Sum = []

	print "REDSHIFT %.4f" %(SnapListZ[i])
	for j in xrange(0, len(bins)-1):
		w = np.where((Mass[i] >= bins[j]) & (Mass[i] < bins[j+1]))[0]
		if (len(w) != 0):
			Photons_Sum.append(Sum_Log(Photons[i][w]))
			print "Photons_Sum[-1]", Photons_Sum[-1], "Bins %.1f-%.1f" %(bins[j], bins[j+1])
		else:
			Photons_Sum.append(nan)


	label = "z = %.2f" %(SnapListZ[i])
	ax.plot(bins_mid, np.log10(Photons_Sum), label = label)
 
    ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval)) 
    ax.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
            
    ax.set_xlim(Limits[0:2])
    ax.set_ylim(Limits[2:4])

    leg = plt.legend(loc=1, numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')

       
    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()


##

def fesc(simulation, SnapListZ, mass, fesc):

    Output_Tag = "fesc_mass_fiducial"

    dilute = 1e100
    ax = plt.subplot(111)


    for i in xrange(0, len(SnapListZ)):
    	if(len(mass[i]) > dilute): 
		w = random.sample(range(0, len(mass[i])), dilute)
	else:
		w = range(0, len(mass[i])) 
	print len(mass[i])
	print "Max mass is %.4f and min is %.4f" %(max(mass[i][w]), min(mass[i][w]))
	label = "z = %.4f" %(SnapListZ[i])
    	ax.scatter(mass[i][w], fesc[i][w], label = label, alpha = 0.5, color = colours[i])

    leg = plt.legend(loc=1, numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')


    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

##

#################################

Simulation = 1 # Set 0 for Mini-Millennium, 1 for My_Simulation, 2 for both (kinda).


binwidth = 0.1
low_mass = 8
high_mass = 12

bins = np.arange(low_mass,high_mass+binwidth, binwidth)
bins_mid = bins + binwidth/2.0

if (Simulation == 0 or Simulation == 2):
    H_Millennium = ReadScripts.ReadHalos('/lustre/projects/p004_swin/jseiler/millennium_trees/trees_063', 0, 7)

    GG_Millennium, Gal_Desc = ReadScripts.ReadGals_SAGE_Photons('./results/millennium_post_processed/NewCode_z0.000', 0, 7, 64)
    G_Merged_Millennium, Merged_Desc = ReadScripts.ReadGals_SAGE_Photons('./results/millennium_post_processed/NewCode_z0.000', 0, 7, 64)
    G_Millennium = ReadScripts.Join_Arrays(GG_Millennium, G_Merged_Millennium, Gal_Desc)

    GG_Millennium2, Gal_Desc = ReadScripts.ReadGals_SAGE_Photons('./results/millennium_post_processed/noreion_z0.000', 0, 7, 64)
    G_Merged_Millennium2, Merged_Desc = ReadScripts.ReadGals_SAGE_Photons('./results/millennium_post_processed/noreion_z0.000', 0, 7, 64)
    G_Millennium2 = ReadScripts.Join_Arrays(GG_Millennium2, G_Merged_Millennium2, Gal_Desc)


    SnapList_Millennium = [15] 
    print "Snapshots analyzing are", SnapList_Millennium

if (Simulation == 1 or Simulation == 2):

    H_MySim = ReadScripts.ReadHalos('/lustre/projects/p004_swin/jseiler/Rockstar_output/Halos_noLL/Ltrees/lhalotree.bin', 0, 26)
    H_MySim2 = ReadScripts.ReadHalos('/lustre/projects/p004_swin/jseiler/Rockstar_output/1024_Halos_noLL/Ltrees/lhalotree.bin', 0, 124)
    #H_MySim2 = ReadScripts.ReadHalos('/lustre/projects/p004_swin/jseiler/Rockstar_output/Halos_noLL/Ltrees/lhalotree.bin', 0, 26)
    
    # 512 
    GG_MySim, Gal_Desc = ReadScripts.ReadGals_SAGE_Photons('/lustre/projects/p004_swin/jseiler/SAGE_output/512/fiducial_z5.038', 0, 26, 100)
    G_Merged_MySim, Merged_Desc = ReadScripts.ReadGals_SAGE_Photons('/lustre/projects/p004_swin/jseiler/SAGE_output/512/fiducial_MergedGalaxies', 0, 26, 100)
    G_MySim = ReadScripts.Join_Arrays(GG_MySim, G_Merged_MySim, Gal_Desc)

    # 1024
    #GG_MySim, Gal_Desc = ReadScripts.ReadGals_SAGE_Photons('/lustre/projects/p004_swin/jseiler/SAGE_output/1024/noreion_reionmine_DRAGONS_z5.000', 0, 124, 101)
    #G_Merged_MySim, Merged_Desc = ReadScripts.ReadGals_SAGE_Photons('/lustre/projects/p004_swin/jseiler/SAGE_output/1024/noreion_reionmine_DRAGONS_MergedGalaxies', 0, 124, 101)
    #G_MySim = ReadScripts.Join_Arrays(GG_MySim, G_Merged_MySim, Gal_Desc)

    # 512
    #GG_MySim2, Gal_Desc = ReadScripts.ReadGals_SAGE_Photons('/lustre/projects/p004_swin/jseiler/SAGE_output/512/noreion_myreion_512_z5.038', 0, 26, 100)
    #G_Merged_MySim2, Merged_Desc = ReadScripts.ReadGals_SAGE_Photons('/lustre/projects/p004_swin/jseiler/SAGE_output/512/noreion_myreion_512_MergedGalaxies', 0, 26, 100)
    #G_MySim2 = ReadScripts.Join_Arrays(GG_MySim, G_Merged_MySim, Gal_Desc)

    # 1024
    GG_MySim2, Gal_Desc = ReadScripts.ReadGals_SAGE_Photons('/lustre/projects/p004_swin/jseiler/SAGE_output/1024/noreion_z5.000', 0, 124, 101)
    G_Merged_MySim2, Merged_Desc = ReadScripts.ReadGals_SAGE_Photons('/lustre/projects/p004_swin/jseiler/SAGE_output/1024/noreion_MergedGalaxies', 0, 124, 101)
    G_MySim2 = ReadScripts.Join_Arrays(GG_MySim2, G_Merged_MySim2, Gal_Desc)

    #SnapList_MySim = np.arange(22,78) 

    SnapList_MySim = [78]
    # Snapshots for z = [6, 7, 8] are [78, 64, 51] 
    print "Snapshots analyzing are", SnapList_MySim



SnapListZ = []

AllVars.Set_Constants()

HaloPart_Low = 41 # Bounds for where we define the cutoff for a 'Dark Matter Halo'. Below this we can't be sure of the results.
HaloPart_High = 51

'''
GG_model2, Gal_Desc = ReadScripts.ReadGals_STARBURST_Millennium('./results/millennium_post_processed/Photons_NoReion_NewSave', 0, 7, 64)
G_Merged_model2, Merged_Desc = ReadScripts.ReadGals_STARBURST_Millennium('./results/millennium_post_processed/Photons_Merged_NoReion_NewSave', 0, 7, 64)
G_model2 = ReadScripts.Join_Arrays(GG_model2, G_Merged_model2, Gal_Desc)

GG_model3, Gal_Desc = ReadScripts.ReadGals_STARBURST_Millennium('./results/millennium_post_processed/Photons_NoSN_NewSave', 0, 7, 64)
G_Merged_model3, Merged_Desc = ReadScripts.ReadGals_STARBURST_Millennium('./results/millennium_post_processed/Photons_Merged_NoSN_NewSave', 0, 7, 64)
G_model3 = ReadScripts.Join_Arrays(GG_model3, G_Merged_model3, Gal_Desc)
'''
## Halo Initialization ##

w_H_MySim = []
mass_H_MySim = []
Photons_HI_H_MySim = []
Photons_HI_Tot_H_MySim = []
SourceEfficiency_H_MySim = 10

w_H_MySim2 = []
mass_H_MySim2 = []
Photons_HI_H_MySim2 = []
Photons_HI_Tot_H_MySim2 = []
SourceEfficiency_H_MySim2 = 10


## MySim Initialization ## 

w_G_MySim = []
mass_G_MySim = []
Photons_HI_G_MySim = []
Photons_HI_Tot_G_MySim = [] 
SFR_G_MySim = []
sSFR_G_MySim = []
Metallicity_Tremonti_G_MySim = []
HaloPart_MySim = []
fesc_MySim = 0.05

mass_Central_MySim = []
Photons_HI_Central_MySim = []
Photons_HI_Tot_Central_MySim = []
SourceEfficiency_MySim = 1 

SnapListZ_MySim = []

w_Ionized_MySim = []
MfiltGnedin_MySim = []
MfiltSobacchi_MySim = []

LymanAlpha_MySim = []
fesc_LymanAlpha_MySim = fesc_MySim

## MySim Model 2 ##

w_G_MySim2 = []
mass_G_MySim2 = []
Photons_HI_G_MySim2 = []
Photons_HI_Tot_G_MySim2 = []
SFR_G_MySim2 = []
sSFR_G_MySim2 = []
Metallicity_Tremonti_G_MySim2 = []
HaloPart_MySim2 = []
fesc_MySim2 = 0.05

mass_Central_MySim2 = []
Photons_HI_Central_MySim2 = []
Photons_HI_Tot_Central_MySim2 = []
SourceEfficiency_MySim2 = 1

w_Ionized_MySim2 = []
MfiltGnedin_MySim2 = []
MfiltSobacchi_MySim2 = []

LymanAlpha_MySim2 = []
fesc_LymanAlpha_MySim2 = fesc_MySim2

fesc_Kimm_MySim2 = []

## Millennium Initialization ##

w_H_Millennium = []
mass_H_Millennium = []
Photons_HI_H_Millennium = []
Photons_HI_Tot_H_Millennium = []
SourceEfficiency_H_Millennium = 50

w_G_Millennium = []
mass_G_Millennium = []
Photons_HI_G_Millennium = []
Photons_HI_Tot_G_Millennium = []
SFR_G_Millennium = []
sSFR_G_Millennium = []
Metallicity_Tremonti_G_Millennium = []
fesc_Millennium = 0.675

mass_Central_Millennium = []
Photons_HI_Central_Millennium = []
Photons_HI_Tot_Central_Millennium = []
SourceEfficiency_Millennium = 50

SnapListZ_Millennium = []

## Model 2 Initialization ##


w_H_Millennium2 = []
mass_H_Millennium2 = []
Photons_HI_H_Millennium2 = []
Photons_HI_Tot_H_Millennium2 = []
SourceEfficiency_H_Millennium2 = 10

w_G_Millennium2 = []
mass_G_Millennium2 = []
Photons_HI_G_Millennium2 = []
Photons_HI_Tot_G_Millennium2 = []
SFR_G_Millennium2 = []
sSFR_G_Millennium2 = []
Metallicity_Tremonti_G_Millennium2 = []
fesc_Millennium2 = 0.675

mass_Central_Millennium2 = []
Photons_HI_Central_Millennium2 = []
Photons_HI_Tot_Central_Millennium2 = []
SourceEfficiency_Millennium2 = 10

##

if (Simulation == 1 or Simulation == 2):
    for i in xrange(0, len(SnapList_MySim)):
      print i
    ###### MYSIM CALCULATIONS #####


      AllVars.Set_Params_Mysim()

      ## Halo Calculations ##

      w_H_MySim.append(np.where((H_MySim.SnapNum == SnapList_MySim[i]) & (H_MySim.Mvir > 0.0))[0])  
      mass_H_MySim.append(np.log10(H_MySim.Mvir[w_H_MySim[i]] * 1.0e10 / AllVars.Hubble_h))
      Photon_Factor = AllVars.Solar_Mass * SourceEfficiency_H_MySim * AllVars.BaryonFrac / AllVars.Ionized_Mass_H / AllVars.Proton_Mass / (AllVars.Lookback_Time[SnapList_MySim[i]] - AllVars.Lookback_Time[SnapList_MySim[i]+1]) / AllVars.Sec_Per_Megayear / 1e3
      Photons_HI_H_MySim.append(np.log10(H_MySim.Mvir[w_H_MySim[i]]* 1.0e10 / AllVars.Hubble_h * Photon_Factor))
      Photons_HI_Tot_H_MySim.append(np.log10(sum(10**Photons_HI_H_MySim[i])))

      ## Model 1 Calculations ##

      #w_G_MySim.append(np.where((G_MySim.GridHistory[:, SnapList_MySim[i]] != -1) & (G_MySim.GridStellarMass[:, SnapList_MySim[i]] > 0.0) & (G_MySim.GridSFR[:, SnapList_MySim[i]] > 0.0) & (G_MySim.GridCentralGalaxyMass[:,SnapList_MySim[i]] > 0.0) & (G_MySim.GridZ[:, SnapList_MySim[i]] > 0.0))[0])
      w_G_MySim.append(np.where((G_MySim.GridHistory[:, SnapList_MySim[i]] != -1) & (G_MySim.GridStellarMass[:, SnapList_MySim[i]] > 0.0) & (G_MySim.GridSFR[:, SnapList_MySim[i]] > 0.0) & (G_MySim.Photons_HI[:, SnapList_MySim[i]] > 0.0))[0])
      mass_G_MySim.append(np.log10(G_MySim.GridStellarMass[w_G_MySim[i], SnapList_MySim[i]] * 1.0e10 / AllVars.Hubble_h))
      Photons_HI_G_MySim.append(G_MySim.Photons_HI[w_G_MySim[i], SnapList_MySim[i]] + np.log10(fesc_MySim))
    
      Photons_HI_Tot_G_MySim.append(np.log10(Sum_Log(Photons_HI_G_MySim[i])))
 
      SFR_G_MySim.append(np.log10(G_MySim.GridSFR[w_G_MySim[i], SnapList_MySim[i]]))
      sSFR_G_MySim.append(SFR_G_MySim[i] - mass_G_MySim[i])
      HaloPart_MySim.append(G_MySim.Len[w_G_MySim[i]])
      

      Metallicity_Tremonti_G_MySim.append(np.log10(G_MySim.GridZ[w_G_MySim[i], SnapList_MySim[i]] / 0.02) + 9.0)

      print SFR_G_MySim[i]
      print "There were %d galaxies for snapshot %d (Redshift %.4f) for model1." %(len(w_G_MySim[i]), SnapList_MySim[i], AllVars.SnapZ[SnapList_MySim[i]])

      mass_Central_MySim.append(np.log10(G_MySim.GridCentralGalaxyMass[w_G_MySim[i], SnapList_MySim[i]] * 1.0e10 / AllVars.Hubble_h))
      Photon_Factor = AllVars.Solar_Mass * SourceEfficiency_MySim * AllVars.BaryonFrac / AllVars.Ionized_Mass_H / AllVars.Proton_Mass / (AllVars.Lookback_Time[SnapList_MySim[i]] - AllVars.Lookback_Time[SnapList_MySim[i]+1]) / AllVars.Sec_Per_Megayear / 1e3
      Photons_HI_Central_MySim.append(np.log10(10**(mass_Central_MySim[i]) * Photon_Factor))
      Photons_HI_Tot_Central_MySim.append(np.log10(sum(10**Photons_HI_Central_MySim[i])))

      w_Ionized_MySim.append(np.where((G_MySim.GridHistory[:, SnapList_MySim[i]] != -1) & (G_MySim.GridStellarMass[:, SnapList_MySim[i]] > 0.0) & (G_MySim.GridSFR[:, SnapList_MySim[i]] > 0.0) & (G_MySim.Photons_HI[:, SnapList_MySim[i]] > 0.0) & (G_MySim.MfiltSobacchi[:, SnapList_MySim[i]] > 1.0))[0])

      MfiltGnedin_MySim.append(G_MySim.MfiltGnedin[w_G_MySim[i], SnapList_MySim[i]])
      #MfiltSobacchi_MySim.append(G_MySim.MfiltSobacchi[w_Ionized_MySim[i], SnapList_MySim[i]])
      MfiltSobacchi_MySim.append(G_MySim.MfiltSobacchi[w_G_MySim[i], SnapList_MySim[i]])

      LymanAlpha_MySim.append(np.log10(0.68*(1.0-fesc_MySim)*fesc_LymanAlpha_MySim * (AllVars.LymanAlpha_Energy* AllVars.eV_to_erg)) + Photons_HI_Central_MySim[i]) 
 
      ## Model 2 Calculations ##

      ## Halos ##

      w_H_MySim2.append(np.where((H_MySim2.SnapNum == SnapList_MySim[i]) & (H_MySim2.Mvir > 0.0))[0]) 
      mass_H_MySim2.append(np.log10(H_MySim2.Mvir[w_H_MySim2[i]] * 1.0e10 / AllVars.Hubble_h))
      Photon_Factor = AllVars.Solar_Mass * SourceEfficiency_H_MySim * AllVars.BaryonFrac / AllVars.Ionized_Mass_H / AllVars.Proton_Mass / (AllVars.Lookback_Time[SnapList_MySim[i]] - AllVars.Lookback_Time[SnapList_MySim[i]+1]) / AllVars.Sec_Per_Megayear / 1e3
      Photons_HI_H_MySim2.append(np.log10(H_MySim2.Mvir[w_H_MySim2[i]]* 1.0e10 / AllVars.Hubble_h * Photon_Factor))
      Photons_HI_Tot_H_MySim2.append(np.log10(sum(10**Photons_HI_H_MySim2[i])))

      ## Galaxies ##
 
      w_G_MySim2.append(np.where((G_MySim2.GridHistory[:, SnapList_MySim[i]] != -1) & (G_MySim2.GridStellarMass[:, SnapList_MySim[i]] > 0.0) & (G_MySim2.GridSFR[:, SnapList_MySim[i]] > 0.0) & (G_MySim2.Photons_HI[:, SnapList_MySim[i]] > 0.0))[0])
      #w_G_MySim2.append(np.where((G_MySim2.GridHistory[:, SnapList_MySim[i]] != -1) & (G_MySim2.GridStellarMass[:, SnapList_MySim[i]] > 0.0) & (G_MySim2.GridSFR[:, SnapList_MySim[i]] > 0.0) & (G_MySim2.GridCentralGalaxyMass[:,SnapList_MySim[i]] > 0.0) & (G_MySim2.GridZ[:, SnapList_MySim[i]] > 0.0))[0])
      mass_G_MySim2.append(np.log10(G_MySim2.GridStellarMass[w_G_MySim2[i], SnapList_MySim[i]] * 1.0e10 / AllVars.Hubble_h))
      Photons_HI_G_MySim2.append(G_MySim2.Photons_HI[w_G_MySim2[i], SnapList_MySim[i]] + np.log10(fesc_MySim2)) 
      Photons_HI_Tot_G_MySim2.append(np.log10(Sum_Log(Photons_HI_G_MySim2[i])))

      SFR_G_MySim2.append(np.log10(G_MySim2.GridSFR[w_G_MySim2[i], SnapList_MySim[i]])) 
      
      sSFR_G_MySim2.append(SFR_G_MySim2[i] - mass_G_MySim2[i])
      HaloPart_MySim2.append(G_MySim2.Len[w_G_MySim2[i]])
      Metallicity_Tremonti_G_MySim2.append(np.log10(G_MySim2.GridZ[w_G_MySim2[i], SnapList_MySim[i]] / 0.02) + 9.0)

      mass_Central_MySim2.append(np.log10(G_MySim2.GridCentralGalaxyMass[w_G_MySim2[i], SnapList_MySim[i]] * 1.0e10 / AllVars.Hubble_h))
      Photon_Factor = AllVars.Solar_Mass * SourceEfficiency_MySim2 * AllVars.BaryonFrac / AllVars.Ionized_Mass_H / AllVars.Proton_Mass / (AllVars.Lookback_Time[SnapList_MySim[i]] - AllVars.Lookback_Time[SnapList_MySim[i]+1]) / AllVars.Sec_Per_Megayear / 1e3
      Photons_HI_Central_MySim2.append(np.log10(10**(mass_Central_MySim2[i]) * Photon_Factor))
      Photons_HI_Tot_Central_MySim2.append(np.log10(sum(10**Photons_HI_Central_MySim2[i]))) 

      print "There were %d galaxies for snapshot %d (Redshift %.4f) for model2." %(len(w_G_MySim2[i]), SnapList_MySim[i], AllVars.SnapZ[SnapList_MySim[i]])

      w_Ionized_MySim2.append(np.where((G_MySim2.GridHistory[:, SnapList_MySim[i]] != -1) & (G_MySim2.GridStellarMass[:, SnapList_MySim[i]] > 0.0) & (G_MySim2.GridSFR[:, SnapList_MySim[i]] > 0.0) & (G_MySim2.Photons_HI[:, SnapList_MySim[i]] > 0.0) & (G_MySim2.MfiltSobacchi[:, SnapList_MySim[i]] < 1.0))[0])

      MfiltGnedin_MySim2.append(G_MySim2.MfiltGnedin[w_G_MySim2[i], SnapList_MySim[i]])
      MfiltSobacchi_MySim2.append(G_MySim2.MfiltSobacchi[w_Ionized_MySim2[i], SnapList_MySim[i]])

      LymanAlpha_MySim2.append(np.log10(0.68*(1.0-fesc_MySim2)*fesc_LymanAlpha_MySim2 * (AllVars.LymanAlpha_Energy* AllVars.eV_to_erg)) + Photons_HI_Central_MySim2[i]) 

      fesc_Kimm_MySim2.append(10**(1.00 - 0.2 * mass_Central_MySim2[i]))

      SnapListZ_MySim.append(AllVars.SnapZ[SnapList_MySim[i]])
      SnapListZ.append(AllVars.SnapZ[SnapList_MySim[i]])

if (Simulation == 0 or Simulation == 2):
    for i in xrange(0, len(SnapList_Millennium)):

    ##### MILLENNIUM CALCULATIONS ####
    
   
      AllVars.Set_Params_MiniMill()
      ## Halos ##

      w_H_Millennium.append(np.where((H_Millennium.SnapNum == SnapList_Millennium[i]) & (H_Millennium.Mvir > 0.0))[0]) 
      mass_H_Millennium.append(np.log10(H_Millennium.Mvir[w_H_Millennium[i]] * 1.0e10 / AllVars.Hubble_h))
      #Photon_Factor = AllVars.Solar_Mass * SourceEfficiency_H_Millennium * AllVars.BaryonFrac / AllVars.Ionized_Mass_H / AllVars.Proton_Mass / (AllVars.Lookback_Time[SnapList_Millennium[i]] - AllVars.Lookback_Time[SnapList_Millennium[i]+1]) / AllVars.Sec_Per_Megayear / 1e3
      #Photons_H_Millennium.append(np.log10(H_Millennium.Mvir[w_H_Millennium[i]]* 1.0e10 / AllVars.Hubble_h * Photon_Factor))
      #Photons_Tot_H_Millennium.append(np.log10(sum(10**Photons_H_Millennium[i])))

      SnapListZ.append(AllVars.SnapZ[SnapList_Millennium[i]])
      SnapListZ_Millennium.append(AllVars.SnapZ[SnapList_Millennium[i]])

      ## Model 1 Calculations ##

      w_G_Millennium.append(np.where((G_Millennium.GridHistory[:, SnapList_Millennium[i]] != -1) & (G_Millennium.GridStellarMass[:, SnapList_Millennium[i]] > 0.0) & (G_Millennium.GridSFR[:, SnapList_Millennium[i]] > 0.0) & (G_Millennium.GridCentralGalaxyMass[:,SnapList_Millennium[i]] > 0.0))[0])
      mass_G_Millennium.append(np.log10(G_Millennium.GridStellarMass[w_G_Millennium[i], SnapList_Millennium[i]] * 1.0e10 / AllVars.Hubble_h))
      Photons_HI_G_Millennium.append(np.log10(10**(G_Millennium.Photons_HI[w_G_Millennium[i], SnapList_Millennium[i]])*fesc_Millennium))
      Photons_HI_Tot_G_Millennium.append(np.log10(sum(10**Photons_HI_G_Millennium[i])))
      SFR_G_Millennium.append(np.log10(G_Millennium.GridSFR[w_G_Millennium[i], SnapList_Millennium[i]]))
      sSFR_G_Millennium.append(SFR_G_Millennium[i] - mass_G_Millennium[i])

      mass_Central_Millennium.append(np.log10(G_Millennium.GridCentralGalaxyMass[w_G_Millennium[i], SnapList_Millennium[i]] * 1.0e10 / AllVars.Hubble_h))
      Photon_Factor = AllVars.Solar_Mass * SourceEfficiency_Millennium * AllVars.BaryonFrac / AllVars.Ionized_Mass_H / AllVars.Proton_Mass / (AllVars.Lookback_Time[SnapList_Millennium[i]] - AllVars.Lookback_Time[SnapList_Millennium[i]+1]) / AllVars.Sec_Per_Megayear / 1e3
      Photons_HI_Central_Millennium.append(np.log10(10**(mass_Central_Millennium[i]) * Photon_Factor))
      Photons_HI_Tot_Central_Millennium.append(np.log10(sum(10**Photons_HI_Central_Millennium[i])))

      ## Model 2 Calculations ##

      w_G_Millennium2.append(np.where((G_Millennium2.GridHistory[:, SnapList_Millennium[i]] != -1) & (G_Millennium2.GridStellarMass[:, SnapList_Millennium[i]] > 0.0) & (G_Millennium2.GridSFR[:, SnapList_Millennium[i]] > 0.0) & (G_Millennium2.GridCentralGalaxyMass[:,SnapList_Millennium[i]] > 0.0))[0])
      mass_G_Millennium2.append(np.log10(G_Millennium2.GridStellarMass[w_G_Millennium2[i], SnapList_Millennium[i]] * 1.0e10 / AllVars.Hubble_h))
      Photons_HI_G_Millennium2.append(np.log10(10**(G_Millennium2.Photons_HI[w_G_Millennium2[i], SnapList_Millennium[i]])*fesc_Millennium2))
      Photons_HI_Tot_G_Millennium2.append(np.log10(sum(10**Photons_HI_G_Millennium2[i])))
      SFR_G_Millennium2.append(np.log10(G_Millennium2.GridSFR[w_G_Millennium2[i], SnapList_Millennium[i]]))
      sSFR_G_Millennium2.append(SFR_G_Millennium2[i] - mass_G_Millennium2[i])

      mass_Central_Millennium2.append(np.log10(G_Millennium2.GridCentralGalaxyMass[w_G_Millennium2[i], SnapList_Millennium[i]] * 1.0e10 / AllVars.Hubble_h))
      Photon_Factor = AllVars.Solar_Mass * SourceEfficiency_Millennium2 * AllVars.BaryonFrac / AllVars.Ionized_Mass_H / AllVars.Proton_Mass / (AllVars.Lookback_Time[SnapList_Millennium[i]] - AllVars.Lookback_Time[SnapList_Millennium[i]+1]) / AllVars.Sec_Per_Megayear / 1e3
      Photons_HI_Central_Millennium2.append(np.log10(10**(mass_Central_Millennium2[i]) * Photon_Factor))
      Photons_HI_Tot_Central_Millennium2.append(np.log10(sum(10**Photons_HI_Central_Millennium2[i])))


      '''

  ## Model 3 Calculations ##

  w_G_model3.append(np.where((G_model3.GridHistory[:, SnapList[i]] != -1) & (G_model3.GridStellarMass[:, SnapList[i]] > 0.0) & (G_model3.GridSFR[:, SnapList[i]] > 0.0) & (G_model3.GridCentralGalaxyMass[:,SnapList[i]] > 0.0))[0])
  mass_G_model3.append(np.log10(G_model3.GridStellarMass[w_G_model3[i], SnapList[i]] * 1.0e10 / AllVars.Hubble_h))
  Photons_G_model3.append(np.log10(10**(G_model3.Photons[w_G_model3[i], SnapList[i]])*fesc_model3))
  Photons_Tot_G_model3.append(np.log10(sum(10**Photons_G_model3[i])))
  SFR_G_model3.append(np.log10(G_model3.GridSFR[w_G_model3[i], SnapList[i]]))
  sSFR_G_model3.append(SFR_G_model3[i] - mass_G_model3[i])

  mass_Central_model3.append(np.log10(G_model3.GridCentralGalaxyMass[w_G_model3[i], SnapList[i]] * 1.0e10 / AllVars.Hubble_h))
  Photon_Factor = AllVars.Solar_Mass * SourceEfficiency_model3 * AllVars.BaryonFrac / AllVars.Ionized_Mass_H / AllVars.Proton_Mass / (AllVars.Lookback_Time[SnapList[i]] - AllVars.Lookback_Time[SnapList[i]+1]) / AllVars.Sec_Per_Megayear / 1e3
  Photons_Central_model3.append(np.log10(10**(mass_Central_model3[i]) * Photon_Factor))
  Photons_Tot_Central_model3.append(np.log10(sum(10**Photons_Central_model3[i])))
      '''

def Calculate_HaloPartStellarMass(ZZ, HaloPart, StellarMass, BoundLow, BoundHigh):

	for i in xrange(0, len(ZZ)):	
		w = np.where((HaloPart[i] > BoundLow) & (HaloPart[i] < BoundHigh))[0]

		Mass = np.mean(10**(StellarMass[i][w]))
		Mass_std = np.std(10**(StellarMass[i][w]))
		count = len(StellarMass[i][w])
		print "Mass %.4e \t Mass_std = %.4e \t Count = %d" %(Mass, Mass_std, count) 

	return Mass


#print "For 512 model:"
#HaloPartStellarMass_MySim = Calculate_HaloPartStellarMass(SnapListZ, HaloPart_MySim, mass_G_MySim, HaloPart_Low, HaloPart_High)

#print "For 1024 model:"
#HaloPartStellarMass_MySim2 = Calculate_HaloPartStellarMass(SnapListZ, HaloPart_MySim2, mass_G_MySim2, HaloPart_Low, HaloPart_High)


#Metallicity(Simulation, SnapListZ, mass_G_MySim, Metallicity_Tremonti_G_model1)
#Photon_Totals(Simulation, [SnapListZ_MySim, SnapListZ_MySim, SnapListZ_MySim, SnapListZ_MySim], [Photons_Tot_Central_MySim, Photons_Tot_G_MySim, Photons_Tot_Central_MySim2, Photons_Tot_G_MySim2], len(SnapList_MySim))
#StellarMassFunction(Simulation, SnapListZ, (mass_G_MySim + mass_G_Millennium + mass_G_MySim2 + mass_G_Millennium2), HaloPartStellarMass_MySim2, len(SnapList_MySim))
#HaloMassFunction(Simulation, SnapListZ, (mass_H_MySim + mass_H_MySim2 + mass_H_Millennium), len(SnapList_MySim)) 
#CentralGalaxy_Comparison(Simulation, SnapListZ_MySim, (mass_Central_MySim2 + mass_Central_MySim2), (Photons_Central_MySim2 + Photons_G_MySim2))
#CentralGalaxy_Comparison_Difference(Simulation, SnapListZ, (mass_Central_MySim + mass_Central_model1), (Photons_Central_model1 + Photons_G_model1))
#CentralGalaxy_Projection(Simulation, SnapListZ, (mass_Central_MySim + mass_Central_MySim2), (Photons_G_MySim + Photons_G_MySim2))
#FilteringMass(Simulation, SnapListZ_MySim, MfiltGnedin_MySim, MfiltSobacchi_MySim, MfiltGnedin_MySim2, mass_Central_MySim, len(SnapList_MySim))
#PlotScripts.Plot_Scatter_SixPanel((mass_Central_model1 + mass_Central_model1 + mass_Central_model3), (Photons_Central_model1 + Photons_G_model1 + Photons_G_model3), 1, 2, ['Halos', 'Fiducial Gal', 'NoSN Gal'], [min_mass_H, max_mass_H, -2, 1], [r'Log Halo Mass [$M_{\odot}$]', r'Log (Halo/Galaxy Ionizing Photons) [s$^{-1}$]'], 2, SnapListZ, 'CentralStellar_Photons_Difference_Fiducial_NoSN', '.png')
#StarFormationRate(Simulation, SnapListZ_MySim, SFR_G_MySim, SFR_G_MySim2)
#PhotonsStellarMass(Simulation, SnapListZ_MySim, mass_G_MySim2, Photons_HI_G_MySim2)
#LymanAlphaLF(Simulation, SnapListZ_MySim, LymanAlpha_MySim2, len(SnapList_MySim))
#PhotonsVsStellarMass(Simulation, SnapListZ_MySim, mass_G_MySim2, Photons_HI_G_MySim2)
fesc(Simulation, SnapListZ_MySim, mass_Central_MySim2, fesc_Kimm_MySim2)
