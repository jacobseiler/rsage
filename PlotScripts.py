#!/usr/bin/env python
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

def Set_Params_Plot():

	global global_labelsize
	global global_fontsize
	global global_legendsize
	global global_linewidth
	global global_tickinterval

	global colors
	global markers
	global linestyles

	global z_plot
	global time_xlim
	global time_tickinterval 

	global_labelsize = 20
	global_fontsize = 20
	global_legendsize = 18
	global_linewidth = 3
	global_tickinterval = 0.25

	matplotlib.rcdefaults()
	plt.rc('axes', color_cycle=['k', 'b', 'r', 'g', 'm', 'c','y', '0.5'], labelsize='x-large')
	plt.rc('lines', linewidth=global_linewidth)
	# plt.rc('font', variant='monospace')
	plt.rc('legend', numpoints=1, fontsize='x-large')
	plt.rc('text', usetex=True)

	plt.rc('xtick', labelsize=global_fontsize)
	plt.rc('ytick', labelsize=global_fontsize)
	plt.rc('text', usetex=True)

	np.set_printoptions(formatter={'float': lambda x: "{0:0.10e}".format(x)})

	colors = ['r', 'b', 'g', 'm', 'c', 'k']
	markers = ['x', 'o', '^', 's', 'D']
	linestyles = ['-', '--', '-.', ':']
	z_plot = np.arange(6, 14)  #Range of redshift we wish to plot.
	time_xlim = [315, 930]
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
		print "The length of x_data is %d, the length of y_data is %d and the length of y_std is %d." %(len(x_data), len(y_data), len(y_std))
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
			print "The length of the snapshot_labels array is %d but the number of snapshots you're trying to plot is %d" %(len(snapshot_labels), len(x_data))
			raise ValueError("These should have the same length.")

		for snapshot_idx in xrange(0, len(x_data)):
			ax.plot(x_data[snapshot_idx], y_data[snapshot_idx], color = PlotScripts.colors[snapshot_idx], linestyle = PlotScripts.linestyles[0], rasterized = True, label = snapshot_labels[snapshot_idx], linewidth = PlotScripts.global_linewidth)

	## The third case is we have multiple snapshots over multiple modles that we wish to plot; our input data is a 3D array. ##
	## Data is of the form [[[model0_snap0_point0, ..., model0_snap0_pointN], ..., [model0_snapN_point0, ..., model0_snapN_pointN]], ..., [[modelN_snap0_point0, ..., modelN_snap0_pointN], ..., [modelN_snapN_point0, ..., modelN_snapN_pointN]]]. ##
	if dimension == 2: 
		
		if(len(model_labels) != len(x_data)):
			print "The length of the model_labels array is %d but the number of models you're trying to plot is %d." %(len(model_labels), len(x_data))
			raise ValueError("These should have the same length.")

		if(len(snapshot_labels) != len(x_data[0])):
			print "The length of the snapshot_labels array is %d but the number of snapshots you're trying to plot is %d" %(len(snapshot_labels), len(x_data[0]))
			raise ValueError("These should have the same length.")

		for model_number in xrange(0, len(x_data)):
			for snapshot_idx in xrange(0, len(x_data[model_number])):
				ax.plot(x_data[model_number][snapshot_idx], y_data[model_number][snapshot_idx], color = PlotScripts.colors[snapshot_idx], linestyle = PlotScripts.linestyles[model_number], rasterized = True, label = snapshot_labels[snapshot_idx], linewidth = PlotScripts.global_linewidth)

		for model_number in xrange(0, len(x_data)):
			ax.plot(np.nan, np.nan, color = 'k', linestyle = PlotScripts.linestyles[model_number], rasterized = True, label = model_labels[model_number], linewidth = PlotScripts.global_linewidth)
		

	return ax

