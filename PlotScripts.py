#!/usr/bin/env python
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')

import os
import numpy as np
import pylab as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from numpy import *
from random import sample, seed
import math
import random
import csv
from io import StringIO
from collections import Counter
from matplotlib.colors import LogNorm
import time
import matplotlib.ticker as mtick

import AllVars
import ObservationalData as Obs

Chabrier_to_Salpeter = 1.8  

def h_converter_x(hubble_h_observational, hubble_h_simulation):
    return -2.0 * np.log10(hubble_h_simulation / hubble_h_observational)

def h_converter_y(hubble_h_observational, hubble_h_simulation):
    return pow(hubble_h_simulation / hubble_h_observational, 3.0)


def Set_Params_Plot():

    global global_labelsize
    global global_fontsize
    global global_legendsize
    global global_linewidth
    global global_tickinterval
    global global_tickwidth
    global global_ticklength
    global global_errorwidth
    global global_axiswidth

    global colors
    global markers
    global linestyles

    global z_plot
    global time_xlim
    global time_tickinterval 

    global_labelsize = 30
    global_fontsize = 22
    global_legendsize = 20
    global_linewidth = 3
    global_tickinterval = 0.25
    global_tickwidth = 2
    global_ticklength = 6
    global_errorwidth = 3
    global_axiswidth = 2

    matplotlib.rcdefaults()
    plt.rc('axes', color_cycle=['k', 'b', 'r', 'g', 'm', 'c','y', '0.5'], labelsize='x-large')
    plt.rc('lines', linewidth=global_linewidth)
    plt.rc('font', weight='bold') 
    plt.rc('legend', numpoints=1, fontsize='x-large')
    plt.rc('text', usetex=True)

    plt.rc('xtick', labelsize=global_fontsize)
    plt.rc('ytick', labelsize=global_fontsize)
    plt.rc('text', usetex=True)

    np.set_printoptions(formatter={'float': lambda x: "{0:0.10e}".format(x)})

    colors = ["#f03b20", "#c51b8a", "#2b8cbe", "#31a354"]
    #colors = ['#1b9e77','#d95f02','#7570b3', '#f03b20', '#2c7fb8']
    #colors = ['#2c7fb8','#f03b20', '#f768a1', "r", "b", "m"]
    #colors = ['k', 'r', 'b', 'g', 'm', 'c', 'k']
    markers = ['x', 'o', '^', 's', 'D']
    linestyles = ['-', '--', '-.', ':', '-']
    z_plot = np.arange(6, 15, 1)  #Range of redshift we wish to plot.
    time_xlim = [290, 930]
    time_tickinterval = 25

def plot_xy(ax, x_data, y_data, y_std, snapshot_labels, model_labels):
	'''
	Plots a simple x-y line for the given data.
	Accepts multiple snapshots and multiple models in nested numpy arrays.
	Parameters
	----------
	ax : Axis object.
		Axis we are plotting on. 
	x_data, y_data, y_std : Nested `np.darray'.  Assumes the inner array is defined as np.array([....])
		The data that we wish to plot.
		Function accepts a nested array that can plot multiple snapshots over multiple models. 
		Data can be in one of the following forms:
			Single Line : np.array([point0, point1, ..., pointN]). 
			Multiple Lines with different colours : [np.array([snap0_point0, snap0_point1, ... , snap0_pointN]) , ... , np.array([snapN_point0, snapN_point1, ..., snapN_pointN])]
			Multiple Lines of different colours, each with multiple models with different linestypes : [[np.array([model0_snap0_point0, ..., model0_snap0_pointN]), ..., np.array([model0_snapN_point0, ..., model0_snapN_pointN])], ..., [np.array([modelN_snap0_point0, ..., modelN_snap0_pointN]), ..., np.array([modelN_snapN_point0, ..., modelN_snapN_pointN])]].
	snapshot_labels, model_labels : `np.darray' of strings.
		Array that contains the labels for each of the different snapshots (different coloured lines) and models (different linestyles)
				
	Returns
	-------
	ax : Axis object.
		Data with the x-y line plotted on it.
	'''	


	if((len(x_data) != len(y_data)) or (len(x_data) != len(y_std)) or (len(y_data) != len(y_std))):
		print("The length of x_data is %d, the length of y_data is %d and the length of y_std is %d." %(len(x_data), len(y_data), len(y_std)))
		raise ValueError("Each of these need to be equal to each other.")

	
	dimension = AllVars.depth(x_data) # Determines the dimension of the input data.

	## Since I want this script to be able to plot multiple snapshots and even multiple models we need to set up cases. ## 

	## The first case is where we are simply plotting a single snapshot. ##
	## Data is of the form [point0, point1, ..., pointN]. ##
	if dimension == 0:			
		ax.plot(x_data, y_data, color = PlotScripts.colors[0], linestyle = PlotScripts.linestyles[0], rasterized = True, label = snapshot_labels[0], linewidth = PlotScripts.global_linewidth) 

	## The second case is where we have multiple snapshots that we are plotting at; our input data is a 2D array. ##
	## Data is of the form [[snap0_point0, snap0_point1, ... , snap0_pointN] , ... , [snapN_point0, snapN_point1, ..., snapN_pointN]]. ##
	if dimension == 1:	
		if(len(snapshot_labels) != len(x_data)):
			print("The length of the snapshot_labels array is %d but the number of snapshots you're trying to plot is %d" %(len(snapshot_labels), len(x_data)))
			raise ValueError("These should have the same length.")

		for snapshot_idx in xrange(0, len(x_data)):
			ax.plot(x_data[snapshot_idx], y_data[snapshot_idx], color = PlotScripts.colors[snapshot_idx], linestyle = PlotScripts.linestyles[0], rasterized = True, label = snapshot_labels[snapshot_idx], linewidth = PlotScripts.global_linewidth)

	## The third case is we have multiple snapshots over multiple modles that we wish to plot; our input data is a 3D array. ##
	## Data is of the form [[[model0_snap0_point0, ..., model0_snap0_pointN], ..., [model0_snapN_point0, ..., model0_snapN_pointN]], ..., [[modelN_snap0_point0, ..., modelN_snap0_pointN], ..., [modelN_snapN_point0, ..., modelN_snapN_pointN]]]. ##
	if dimension == 2: 
		
		if(len(model_labels) != len(x_data)):
			print("The length of the model_labels array is %d but the number of models you're trying to plot is %d." %(len(model_labels), len(x_data)))
			raise ValueError("These should have the same length.")

		if(len(snapshot_labels) != len(x_data[0])):
			print("The length of the snapshot_labels array is %d but the number of snapshots you're trying to plot is %d" %(len(snapshot_labels), len(x_data[0])))
			raise ValueError("These should have the same length.")

		for model_number in xrange(0, len(x_data)):
			for snapshot_idx in xrange(0, len(x_data[model_number])):
				ax.plot(x_data[model_number][snapshot_idx], y_data[model_number][snapshot_idx], color = PlotScripts.colors[snapshot_idx], linestyle = PlotScripts.linestyles[model_number], rasterized = True, label = snapshot_labels[snapshot_idx], linewidth = PlotScripts.global_linewidth)

		for model_number in xrange(0, len(x_data)):
			ax.plot(np.nan, np.nan, color = 'k', linestyle = PlotScripts.linestyles[model_number], rasterized = True, label = model_labels[model_number], linewidth = PlotScripts.global_linewidth)
		

	return ax


def Plot_SMF_z6(ax, errorwidth = 1.5, capsize = 2.0, linewidth = 1.0, alpha = 1.0):
   
    #Obs.Get_Data_SMF()
      
    ax.errorbar(Obs.Gonzalez_SMF_z6[:,0]+h_converter_x(0.6999999, AllVars.Hubble_h)-np.log10(Chabrier_to_Salpeter),\
                10**(Obs.Gonzalez_SMF_z6[:,1])/h_converter_y(0.6999999, AllVars.Hubble_h),\
                yerr= (10**(Obs.Gonzalez_SMF_z6[:,3])/h_converter_y(0.6999999, AllVars.Hubble_h), 10**(Obs.Gonzalez_SMF_z6[:,2])/h_converter_y(0.6999999, AllVars.Hubble_h)),\
                alpha=alpha, elinewidth = errorwidth, lw=linewidth, marker=Obs.SMF_markers[0], ls='none', 
                label = r'$\mathbf{Gonzalez \: et \: al. \: 2011}$', color = Obs.SMF_colors[0], capsize = capsize) 

    ax.errorbar(Obs.Duncan_SMF_z6[:,0]+h_converter_x(0.7, AllVars.Hubble_h),\
                Obs.Duncan_SMF_z6[:,1]/h_converter_y(0.7, AllVars.Hubble_h),\
                yerr = [Obs.Duncan_SMF_z6[:,2], Obs.Duncan_SMF_z6[:,3]],\
                capsize = capsize, alpha=alpha, elinewidth = errorwidth, lw=1.0, marker=Obs.SMF_markers[2], 
                ls='none', label = r'$\mathbf{Duncan \: et \: al. \: 2014}$', color = Obs.SMF_colors[2], rasterized=True) 

    ax.errorbar(Obs.Song_SMF_z6[:,0]+h_converter_x(0.7, AllVars.Hubble_h)-np.log10(Chabrier_to_Salpeter),\
                10**(Obs.Song_SMF_z6[:,1])/h_converter_y(0.7, AllVars.Hubble_h),\
                yerr= (10**(Obs.Song_SMF_z6[:,1])/h_converter_y(0.7, AllVars.Hubble_h) - 10**(Obs.Song_SMF_z6[:,3])/h_converter_y(0.7, AllVars.Hubble_h), 10**(Obs.Song_SMF_z6[:,2])/h_converter_y(0.7, AllVars.Hubble_h) - 10**(Obs.Song_SMF_z6[:,1])/h_converter_y(0.7, AllVars.Hubble_h)),\
                xerr = 0.25,\
                capsize = capsize, alpha=alpha, elinewidth = errorwidth, lw=1.0, marker=Obs.SMF_markers[1], ls='none', 
                label = r'$\mathbf{Song \: et \: al. \: 2016}$', color = Obs.SMF_colors[1], rasterized=True) 

def Plot_SMF_z7(ax, errorwidth = 1.5, capsize = 2.0, linewidth = 1.0, alpha = 1.0):

    ax.errorbar(Obs.Gonzalez_SMF_z7[:,0]+h_converter_x(0.6999999, AllVars.Hubble_h)-np.log10(Chabrier_to_Salpeter),\
                10**(Obs.Gonzalez_SMF_z7[:,1])/h_converter_y(0.6999999, AllVars.Hubble_h),\
                yerr= (10**(Obs.Gonzalez_SMF_z7[:,3])/h_converter_y(0.6999999, AllVars.Hubble_h), 10**(Obs.Gonzalez_SMF_z7[:,2])/h_converter_y(0.6999999, AllVars.Hubble_h)),\
                alpha=alpha, elinewidth = errorwidth, lw=linewidth, marker=Obs.SMF_markers[0], ls='none', label = 'Gonzalez et al. 2011', color = Obs.SMF_colors[0], capsize = capsize) 

    ax.errorbar(Obs.Duncan_SMF_z7[:,0]+h_converter_x(0.7, AllVars.Hubble_h),\
                Obs.Duncan_SMF_z7[:,1]/(h_converter_y(0.7, AllVars.Hubble_h)),\
                yerr = [Obs.Duncan_SMF_z7[:,2], Obs.Duncan_SMF_z7[:,3]],\
                capsize = capsize, alpha=alpha, elinewidth = errorwidth, lw=linewidth, marker=Obs.SMF_markers[2], ls='none', label = 'Duncan et al. 2014', color = Obs.SMF_colors[2], rasterized=True) 

    ax.errorbar(Obs.Song_SMF_z7[:,0]+h_converter_x(0.7, AllVars.Hubble_h)-np.log10(Chabrier_to_Salpeter),\
                10**(Obs.Song_SMF_z7[:,1])/h_converter_y(0.7, AllVars.Hubble_h),\
                yerr= (10**(Obs.Song_SMF_z7[:,1])/h_converter_y(0.7, AllVars.Hubble_h) - 10**(Obs.Song_SMF_z7[:,3])/h_converter_y(0.7, AllVars.Hubble_h), 10**(Obs.Song_SMF_z7[:,2])/h_converter_y(0.7, AllVars.Hubble_h) - 10**(Obs.Song_SMF_z7[:,1])/h_converter_y(0.7, AllVars.Hubble_h)),\
                xerr = 0.25,\
                capsize = capsize, alpha=alpha, elinewidth = errorwidth, lw=1.0, marker=Obs.SMF_markers[1], ls='none', label = 'Song et al. 2016', color = Obs.SMF_colors[1], rasterized=True) 

def Plot_SMF_z8(ax, errorwidth = 1.5, capsize = 2.0, linewidth = 1.0, alpha = 1.0):

    ax.errorbar(Obs.Song_SMF_z8[:,0]+h_converter_x(0.7, AllVars.Hubble_h)-np.log10(Chabrier_to_Salpeter),\
                10**(Obs.Song_SMF_z8[:,1])/h_converter_y(0.7, AllVars.Hubble_h),\
                yerr= (10**(Obs.Song_SMF_z8[:,1])/h_converter_y(0.7, AllVars.Hubble_h) - 10**(Obs.Song_SMF_z8[:,3])/h_converter_y(0.7, AllVars.Hubble_h), 10**(Obs.Song_SMF_z8[:,2])/h_converter_y(0.7, AllVars.Hubble_h) - 10**(Obs.Song_SMF_z8[:,1]/h_converter_y(0.7, AllVars.Hubble_h))),\
                xerr = 0.25,\
                capsize = capsize, alpha=alpha, elinewidth = errorwidth, lw=linewidth, marker=Obs.SMF_markers[1], ls='none', label = 'Song et al. 2016', color = Obs.SMF_colors[1], rasterized=True) 
  
def plot_SMBH_z8(ax, linewidth = 1.0, alpha = 1.0):
  
    ax.plot(Obs.Mstar, Obs.Huang_z8_BHSM, alpha = alpha, lw = linewidth, ls = '-.', label = "Huang et al. 2018", color = Obs.SMBH_colors[0], rasterized = True)
   
