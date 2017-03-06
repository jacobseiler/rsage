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


#colours = ['r', 'b', 'g', 'm', 'c', 'k']
colours = ['r', 'b', 'g', 'c', 'm', 'k', 'r', 'b', 'g', 'c', 'm', 'k'] 
markers = ['x', 'o', 'x', 'o', 'D', 's']
linestyles = ['-', '--', '-.', ':']

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


def Plot_Histogram(Data, Labels_Graph, Limits, Bin_Width, Normalization, Weights, Log, Labels_Axis, Legend_Location, Observations, Output_Tag, Output_Format):
    
# Plots a 1D histogram.

## Input ##
# Data is an array of arrays each containing the data to be binned.
# Labels_Graph is an array with each element being the label to be placed on the graph to describe the histogram.
# Limits is an array of the form [Lower_x, Upper_x, Lower_y, Upper_y] which describes the graph limits.
# Bin_Width is the width of the bins.
# Normalization is an array of values which the end histograms should be divided by.
# Weights is either 0 or 1.  0: Denotes that the histogram should be a frequency (count) histogram. 1: Denotes that the histogram should be a probability histogram.
# Log is an array of the form [Log_x, Log_y] where each element is either 0 or 1.  0: Denotes that the axis should be in linear scale. 1: Denotes that the axis should be in log scale.
# Labels_Axis is an array of the form [Label_x, Label_y] which describe the x and y axis labels.
# Legend_Location is a number denoting where the legend should be placed.
#### 1 = Upper Right, 2 = Upper Left, 3 = Lower Left, 4 = Lower Right, 5 = Right, 6 = Center Left, 7 = Center Right, 8 = Lower Center, 9 = Upper Center, 10 = Center. ####
# Output_Tag is the name of the output file.
# Output_Format is the format of the saved graph.

## Output ##
# A file in the form '*Output_Tag**Output_Format' saved in the current working directory.

    plt.figure()  # New figure
    ax = plt.subplot(111)  # 1 plot on the figure
    for i in xrange(0, len(Data)):

        (counts, Bin_Edges, Bin_Middle) = Calculate_Histogram(Data[i], Bin_Width, Weights)
        if (i < 4):
            ls = '-'
        else:
            ls = '--'
        plt.plot(Bin_Middle, counts / Normalization[i], colours[i], linestyle = ls, label = Labels_Graph[i])


    if (Log[0] == 1):
        plt.xscale('log', nonposy='clip')
    if (Log[1] == 1):
        plt.yscale('log', nonposy='clip')


    plt.axis(Limits)

    ax.set_xlabel(Labels_Axis[0])
    ax.set_ylabel(Labels_Axis[1])
    ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1))

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

        plt.errorbar(Song_z6[:,0], 10**Song_z6[:,1], yerr= (10**Song_z6[:,3], 10**Song_z6[:,2]) , alpha=0.75, lw=1.0, marker='o', ls='none', label = 'Song 2015, z = 6', color = 'green')

        Song_z7 = np.array([[7.25, -1.63, -1.63 + 0.54, -1.63 - 0.54],
                            [7.75, -2.07, -2.07 + 0.45, -2.07 - 0.41],
                            [8.25, -2.49, -2.49 + 0.38, -2.49 - 0.32],
                            [8.75, -2.96, -2.96 + 0.32, -2.96 - 0.30],
                            [9.25, -3.47, -3.47 + 0.32, -3.47 - 0.35],
                            [9.75, -4.11, -4.11 + 0.41, -4.11 - 0.57],
                            [10.25, -4.61, -4.61 + 0.72, -4.61 - 0.82],
                            [10.75, -5.24, -5.24 + 0.90, -5.25 - 0.57]], dtype = np.float32)

        plt.errorbar(Song_z7[:,0], 10**Song_z7[:,1], yerr= (10**Song_z7[:,3], 10**Song_z7[:,2]) , alpha=0.75, lw=1.0, marker='o', ls='none', label = 'Song 2015, z = 7', color = 'blue')

        Song_z8 = np.array([[7.25, -1.73, -1.73 + 1.01, -1.73 - 0.84],
                            [7.75, -2.28, -2.28 + 0.84, -2.28 - 0.64],
                            [8.25, -2.88, -2.88 + 0.75, -2.88 - 0.57],
                            [8.75, -3.45, -3.45 + 0.57, -3.45 - 0.60],
                            [9.25, -4.21, -4.21 + 0.63, -4.21 - 0.78],
                            [9.75, -5.31, -5.31 + 1.01, -5.31 - 1.64]], dtype = np.float32)

        plt.errorbar(Song_z8[:,0], 10**Song_z8[:,1], yerr= (10**Song_z8[:,3], 10**Song_z8[:,2]) , alpha=0.75, lw=1.0, marker='o', ls='none', label = 'Song 2015, z = 8', color = 'red')

    leg = plt.legend(loc=Legend_Location, numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')

    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

def Plot_Scatter_TwoHists(Data_x, Data_y, Labels_Graph, Limits, Hist_Limits, Bin_Width, Normalization, Weights, Labels_Axis, Observations, Output_Tag, Output_Format):
    
# This creates a scatter plot in the bottom left of the figure with histogram projections of the x and y axis (located in the top left and bottom right respectively).

## Input ##
# Data_x is an array of arrays that contains the x data.
# Data_y is an array of arrays that contains the y data.
# Labels_Graph is an array of strings that contains the labels for the x and y axis.
# Limits is an array of the form [Lower_x, Upper_x, Lower_y, Upper_y] which describes the graph limits.
# Hist_Limits is an array of the form [Lower_y, Upper_y, Lower_x, Upper_x] which describes the y limits of the UPPER LEFT histogram, and the x  limits of the LOWER RIGHT histogram.
# Bin_Width is the width of the bins.
# Normalization is the value by which the end histogram should be divided by.
# Weights is either 0 or 1.  0: Denotes that the histogram should be a frequency (count) histogram. 1: Denotes that the histogram should be a probability histogram.
# Labels_Axis is an array of the form [Label_x, Label_y] which describe the x and y axis labels.
# Output_Tag is the name of the output file.
# Output_Format is the format of the saved graph.

## Output ##
# A file in the form '*Output_Tag**Output_Format' saved in the current working directory.

    dilute = 500

    label_size = 12
    extra_size = 2
    plt.rc('xtick', labelsize=label_size)
    plt.rc('ytick', labelsize=label_size)
    plt.rc('text', usetex=True)
    tick_interval = 0.25
    
    fig = plt.figure()
    ax3 = fig.add_subplot(223) # Scatter plot in the bottom left.
    ax1 = fig.add_subplot(221) # x histogram projection for the top left.
    ax4 = fig.add_subplot(224) # y histogram projection for the bottom right.

    for i in xrange(0, len(Data_x)):
        print Data_x[i]
        if (len(Data_x[i]) > dilute):
            w = sample(range(len(Data_x[i])), dilute)
        else:
            w = np.arange(0, len(Data_x[i]))
        
        ax3.scatter(Data_x[i][w], Data_y[i][w], color = colours[i], marker = markers[i], label = Labels_Graph[i], alpha = 0.5)

        (counts, Bin_Edges, Bin_Middle) = Calculate_Histogram(Data_x[i], Bin_Width, Weights) # Histogram for top left.
        ax1.plot(Bin_Middle, counts / Normalization, color = colours[i])

        (counts, Bin_Edges, Bin_Middle) = Calculate_Histogram(Data_y[i], Bin_Width, Weights) # Histogram for bottom right.
        ax4.plot(counts / Normalization, Bin_Middle, color = colours[i])


    ax3.set_xlabel(Labels_Axis[0], fontsize = label_size + extra_size)
    ax3.set_ylabel(Labels_Axis[1], fontsize = label_size + extra_size)

    ax3.set_xlim(Limits[0:2])
    ax3.set_ylim(Limits[2:4])

    ax1.set_xlim(Limits[0:2])
    ax1.set_ylim(Hist_Limits[0:2])

    ax4.set_ylim(Limits[2:4])
    ax4.set_xlim(Hist_Limits[2:4])
    
    ax1.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    ax1.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    
    ax3.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    ax3.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    
    ax4.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    ax4.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))

    plt.tight_layout()
    fig.subplots_adjust(hspace = 0.0, wspace = 0.0)

    ax1.tick_params(axis='x', labelbottom = 'off')
    ax1.yaxis.get_major_ticks()[0].label.set_visible(False)

    ax4.tick_params(axis='y', labelleft = 'off')
    ax4.xaxis.get_major_ticks()[0].label.set_visible(False)

    if Weights == 0:
        ax1.set_ylabel("Count", fontsize = label_size + extra_size)
        ax4.set_xlabel("Count", fontsize = label_size + extra_size)

    else:
        ax1.set_ylabel("Probability", fontsize = label_size + extra_size)
        ax4.set_xlabel("Probability", fontsize = label_size + extra_size)

    if Observations == 1:  # Tremonti et al. 2003 Metallicity Data
        w = np.arange(6.0, 13.0, 0.1)
        Zobs_Tremonti = -1.492 + 1.847*w - 0.08026*w*w
        ax3.plot(np.log10((10**w * 1.5 / 1.8)), Zobs_Tremonti, 'k-', lw=2.0, label = 'Tremonti et al. 2003')


    leg = ax3.legend(bbox_to_anchor = (1.6, 2), scatterpoints=1, labelspacing=0.1, handlelength = 2, handletextpad = 0.05)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')
    for legobj in leg.legendHandles:
        legobj.set_linewidth(1.0)

    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

def Plot_Line(Data_x, Data_y, Error, Multi_Lines, Labels_Graph, Limits, Log, Labels_Axis, Legend_Location, Output_Tag, Output_Format):
    
# Plots a simple x-y line.

## Input ##
# Data_x is an array of arrays that contains the x data.
# Data_y is an array of arrays that contains the y data.
# Multi_Lines is an array of the form [Number_Lines, Line_Switch, Line_Labels].  Number_Lines is a number denoting how many different linestyles there are.  Line_Switch is a number denoting when to switch to the next linestyle.  Line_Labels is an array containing the labels for each of the different lines.
# Labels_Graph is an array of strings that contains the labels for the x and y axis.
# Limits is an array of the form [Lower_x, Upper_x, Lower_y, Upper_y] which describes the graph limits.
# Log is an array of the form [Log_x, Log_y] where each element is either 0 or 1.  0: Denotes that the axis should be in linear scale. 1: Denotes that the axis should be in log scale.
# Labels_Axis is an array of the form [Label_x, Label_y] which describe the x and y axis labels.
# Output_Tag is the name of the output file.
# Output_Format is the format of the saved graph.

## Output ##
# A file in the form '*Output_Tag**Output_Format' saved in the current working directory.

    ax = plt.subplot(111)

    Count = 0

    for i in xrange(0, len(Data_x)):
        if (i > Multi_Lines[1]-1 and Multi_Lines[1] != 0): # Check to see if we need to switch linestyle.
            if (i == Multi_Lines[1]): # Resets the colours.
                Count = 0
            ax.plot(Data_x[i], Data_y[i], color = colours[Count], linestyle = linestyles[1])
    #ax.fill_between(Data_x[i], Data_y[i]-Error[i], Data_y[i]+Error[i], alpha = 0.5, color = colours[Count])
        else:
            title = "SFR = %s" %(Labels_Graph[i])
            ax.plot(Data_x[i], Data_y[i], color = colours[Count], label = title, linestyle = linestyles[0])
        Count += 1

    for i in xrange(0, Multi_Lines[0]):
        plt.plot(Limits[0]-1e10, Limits[3]+1e10, ls = linestyles[i], label = Multi_Lines[2][i], color = 'k')

    ax.set_xlabel(Labels_Axis[0])
    ax.set_ylabel(Labels_Axis[1])
    
    ax.set_xlim(Limits[0:2])
    ax.set_ylim(Limits[2:4])

    if (Log[0] == 1):
        plt.xscale('log', nonposy='clip')
    if (Log[1] == 1):
        plt.yscale('log', nonposy='clip')

    leg = plt.legend(loc= Legend_Location, numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')

    ax.text(0.1, 1e11, "Large Scales")
    ax.text(3.0, 1e11, "Small Scales")

    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

###

def Plot_Colourmap(Data, Labels_Axis, Dimension, Limits, Limits_Map, Label_Cbar, Output_Tag, Output_Format):

    if Dimension[0] == 3:
        Data_Projected = Data[:, :, Dimension[1]:Dimension[1]+1].mean(axis = -1)
    else:
        Data_Projected = Data

    ax = plt.subplot(111)

    im = ax.imshow(Data_Projected[:,:], interpolation='none', origin='low', extent =[Limits[0], Limits[1], Limits[2], Limits[3]], vmin = Limits_Map[0], vmax = Limits_Map[1], cmap = 'bone')

    cbar = plt.colorbar(im, ax = ax)
    cbar.set_label(Label_Cbar)
				    
    ax.set_xlabel(Labels_Axis[0])
    ax.set_ylabel(Labels_Axis[1])
                            
    ax.set_xlim(Limits[0:2])
    ax.set_ylim(Limits[2:4])

    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

###

def Plot_Scatter(Data_x, Data_y, Labels_Graph, Limits, Log, Labels_Axis, Legend_Location, Observations, Output_Tag, Output_Format):

    dilute = 5000
    
    fig = plt.figure()
    ax = plt.subplot(111)
  
    for i in xrange(0, len(Data_x)):
      
        if (len(Data_x[i]) > dilute):
            w = sample(range(len(Data_x[i])), dilute)
            ax.scatter(Data_x[i][w], Data_y[i][w], color = colours[i], marker = markers[i], label = Labels_Graph[i], alpha = 0.7)
        elif (len(Data_x[i]) == 1):
            ax.scatter(Data_x[i], Data_y[i], color = colours[i], marker = markers[i], label = Labels_Graph[i], alpha = 0.7)
        else:
            ax.scatter(Data_x[i], Data_y[i], color = colours[i], marker = markers[i], label = Labels_Graph[i], alpha = 0.7)

        #ax.errorbar(np.mean(Data_x[i]), np.mean(Data_y[i]), xerr = np.std(Data_x[i]), yerr = np.std(Data_y[i]), alpha = 0.9, fmt = '^', color = colours[i])

    ax.set_xlabel(Labels_Axis[0], fontsize = label_size + extra_size)
    ax.set_ylabel(Labels_Axis[1], fontsize = label_size + extra_size)
   
    if (Log[0] == 1):
        plt.xscale('log', nonposy='clip')
    if (Log[1] == 1):
        plt.yscale('log', nonposy='clip')
 
    ax.set_xlim(Limits[0:2])
    ax.set_ylim(Limits[2:4])

    if Observations == 1:  # Tremonti et al. 2003 Metallicity Data
        w = np.arange(6.0, 13.0, 0.1)
        Zobs_Tremonti = -1.492 + 1.847*w - 0.08026*w*w
        ax.plot(np.log10((10**w * 1.5 / 1.8)), Zobs_Tremonti, 'k-', lw=2.0, label = 'Tremonti et al. 2003')
    
    ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
    ax.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))

    leg = ax.legend(loc = Legend_Location, numpoints=1, labelspacing=0.1, scatterpoints = 1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')
    
    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

###

def Plot_Scatter_Twox(Data_x, Data_y, Labels_Graph, Limits, Labels_Axis, Legend_Location, z, Output_Tag, Output_Format):


    dilute = 100

    Labels_Graph[0] = Labels_Graph[0] + ' (Count = %d)' %(len(Data_x[0]))
    Labels_Graph[1] = Labels_Graph[1] + ' (Count = %d)' %(len(Data_x[1]))

    fig = plt.figure()
    ax = plt.subplot(111)
    
    for i in xrange(0, len(Data_x)):
        
        if i == 0:
            if (len(Data_x[i]) > dilute):
                w = sample(range(len(Data_x[i])), dilute)
                ax.scatter(Data_x[i][w], Data_y[i][w], color = colours[i], marker = markers[i], label = Labels_Graph[i], alpha = 0.5)
            else:
                w = np.arange(0,len(Data_x[i]))
                ax.scatter(Data_x[i], Data_y[i], color = colours[i], marker = markers[i], label = Labels_Graph[i], alpha = 0.5)
            ax.errorbar(np.mean(Data_x[i]), np.mean(Data_y[i]), xerr = np.std(Data_x[i]), yerr = np.std(Data_y[i]), alpha = 0.9, fmt = '^', color = colours[i])
        else:
            ax2 = ax.twiny()
            if (len(Data_x[i]) > dilute):
                w = sample(range(len(Data_x[i])), dilute)
                ax2.scatter(Data_x[i][w], Data_y[i][w], color = colours[i], marker = markers[i], label = Labels_Graph[i], alpha = 0.5)
            else:
                w = np.arange(0,len(Data_x[i]))
                ax2.scatter(Data_x[i], Data_y[i], color = colours[i], marker = markers[i], label = Labels_Graph[i], alpha = 0.5)
            ax2.errorbar(np.mean(Data_x[i]), np.mean(Data_y[i]), xerr = np.std(Data_x[i]), yerr = np.std(Data_y[i]), alpha = 0.9, fmt = '^', color = colours[i])


    ax2.xaxis.tick_top()
    ax.xaxis.tick_bottom()
    
    ax.set_xlabel(Labels_Axis[0], labelsize = label_size + extra_size)
    ax.set_ylabel(Labels_Axis[1], labelsize = label_size + extra_size)
    ax2.set_xlabel(Labels_Axis[2])
    
    ax.set_xlim(Limits[0:2])
    ax.set_ylim(Limits[2:4])
    
    ax2.set_xlim(Limits[4:6])
    ax2.set_ylim(Limits[2:4])


    ax.scatter(np.nan, np.nan, color = colours[1], marker = markers[1], label = Labels_Graph[1])
    leg = ax.legend(loc = Legend_Location, numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')

    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()



def Plot_SFR_Density(z, SFRD, Output_Tag, Output_Format):


    fig = plt.figure()
    ax = plt.subplot(111)
    
    ax.scatter(z, SFRD, color = 'r', marker = 'x', label = 'Model Galaxies', alpha = 0.8)

    Observations_z = [0.045, 0.3, 0.5, 0.7, 1.0, 0.05, 1.25, 0.3, 0.5, 0.7, 0.9, 1.1, 1.45, 2.1, 3.0, 4.0, 1.125, 1.75, 2.225, 2.3, 3.05, 3.8, 4.9, 5.9, 7.0, 7.9, 7.0, 8.0,
    0.03, 0.03, 0.55, 0.85, 1.15, 1.55, 2.05, 0.55, 0.85, 1.15, 1.55, 2.05, 0.15, 0.375, 5.25, 0.7, 0.9, 1.1, 1.45, 1.85, 2.25, 2.75, 3.6]
                      
    Observations_z_Error = [0.055, 0.1, 0.1, 0.1, 0.2, 0, 0.075, 0.1 ,0.1, 0.1, 0.1, 0.1, 0.25, 0.4, 0.5, 0.5, 0.205, 0.13, 0.145, 0.4, 0.35, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0.15, 0.15, 0.15, 0.25, 0.25, 0.15, 0.15, 0.15, 0.25, 0.25, 0.15, 0.075, 0.075, 0.1, 0.1, 0.1, 0.25, 0.15, 0.25, 0.25, 0.6]
                            
    Observations_SFRD = [-1.82, -1.50, -1.39, -1.20, -1.25, -1.77, -1.75, -1.55, -1.44, -1.24, -0.99, -0.94, -0.95, -0.75, -1.04, -1.69, -1.02, -0.75, -0.87, -0.25, -0.97, -1.29, -1.42, -1.65, -1.79, -2.09, -2.00, -2.21,
         -1.72, -1.95, -1.34, -0.96, -0.89, -0.91, -0.89, -1.22, -1.10, -0.96, -0.94, -0.80, -1.64, -1.42, -1.32, -1.14, -0.94, -0.81, -0.84, -0.86, -0.91, -0.86, -1.36]

    Observations_SFRD_Upper = [0.09, 0.05, 0.15, 0.31, 0.31, 0.08, 0.18, 0.12, 0.10, 0.10, 0.09, 0.09, 0.15, 0.49, 0.26, 0.22, 0.08, 0.12, 0.09, 0.09, 0.11, 0.05, 0.06, 0.08, 0.10, 0.11, 0.10, 0.14,
        0.02, 0.20, 0.22, 0.15, 0.27, 0.17, 0.21, 0.08, 0.10, 0.13, 0.13, 0.18, 0.09, 0.03, 0.05, 0.06, 0.05, 0.04, 0.04, 0.02, 0.09, 0.15, 0.23]
        
    Observations_SFRD_Lower = [0.02, 0.05, 0.08, 0.13, 0.13, 0.09, 0.18, 0.12, 0.10, 0.10, 0.08, 0.09, 0.08, 0.09, 0.15, 0.32, 0.08, 0.12, 0.09, 0.11, 0.15, 0.05, 0.06, 0.08, 0.10, 0.11, 0.11, 0.14,
        0.03, 0.20, 0.11, 0.19, 0.21, 0.21, 0.25, 0.11, 0.13, 0.20, 0.18, 0.15, 0.11, 0.04, 0.05, 0.06, -0.06, 0.05, 0.04, 0.03, 0.12, 0.23, 0.50]

    ax.errorbar(Observations_z, Observations_SFRD, yerr = [Observations_SFRD_Lower, Observations_SFRD_Upper], xerr = [Observations_z_Error, Observations_z_Error], fmt='o', color = 'b', alpha = 0.5, label = 'Madau & Dickinson (2014)')

    ax.set_xlabel(r'z')
    ax.set_ylabel(r'SFR [M$_{\odot}$ yr$^{-1}$ Mpc$^{-3}$]')
    
    ax.set_xlim(-0.0, 12.0)
    ax.set_ylim(-3.0, -0.5)
    
    leg = ax.legend(loc = 1, numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')

    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

def Fit_Linear(Data_x, Data_y, Labels_Graph, Limits, Labels_Axis, Legend_Location, Output_Tag, Output_Format):


    dilute = 10000
    
    fig = plt.figure()
    ax = plt.subplot(111)
    
    for i in xrange(0, len(Data_x)):
        
        if (len(Data_x[i]) > dilute):
            w = sample(range(len(Data_x[i])), dilute)
            ax.scatter(Data_x[i][w], Data_y[i][w], color = colours[i], marker = markers[i], label = Labels_Graph[i], alpha = 0.5)
        else:
            w = np.arange(0,len(Data_x[i]))
            ax.scatter(Data_x[i], Data_y[i], color = colours[i], marker = markers[i], label = Labels_Graph[i], alpha = 0.5)

        LinearFit = np.polyfit(Data_x[i], Data_y[i], 1)
        print "The linear fit y = a*x + b had coefficients a = %.4e and b = %.4e" %(LinearFit[0], LinearFit[1])
        ax.plot(Data_x[i], LinearFit[0]*Data_x[i] + LinearFit[1], color = colours[i])
        

    ax.set_xlabel(Labels_Axis[0])
    ax.set_ylabel(Labels_Axis[1])

    ax.set_xlim(Limits[0:2])
    ax.set_ylim(Limits[2:4])

    leg = ax.legend(loc = Legend_Location, numpoints=1, labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize('medium')

    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

def Plot_Scatter_SixPanel_Twox(Data_x, Data_y, Labels_Graph, Limits, Labels_Axis, Legend_Location, SnapList, Output_Tag, Output_Format):

    dilute = 50
    
    label_size = 10
    plt.rc('xtick', labelsize=label_size)
    plt.rc('ytick', labelsize=label_size)
    plt.rc('text', usetex=True)
    
    tick_interval = 0.5
    
    nrows = 2
    ncols = 3
    
    fig = plt.figure()
    for i in xrange(1, nrows*ncols + 1):
        ax = fig.add_subplot(nrows, ncols, i)
    
        if (len(Data_x[i-1]) > dilute):
            w = sample(range(len(Data_x[i-1])), dilute)
            H = ax.scatter(Data_x[i-1][w], Data_y[i-1][w], color = colours[0], marker = markers[0], label = Labels_Graph[0], alpha = 0.5)
        else:
            H = ax.scatter(Data_x[i-1], Data_y[i-1], color = colours[0], marker = markers[0], label = Labels_Graph[0], alpha = 0.5)
        #ax.errorbar(np.mean(Data_x[i-1]), np.mean(Data_y[i-1]), xerr = np.std(Data_x[i-1]), yerr = np.std(Data_y[i-1]), alpha = 0.9, fmt = '^', color = colours[0])
        
        ax2 = ax.twiny()
        if (len(Data_x[(i-1)+nrows*ncols]) > dilute):
            w = sample(range(len(Data_x[(i-1)+nrows*ncols])), dilute)
            G = ax2.scatter(Data_x[(i-1)+nrows*ncols][w], Data_y[(i-1)+nrows*ncols][w], color = colours[1], marker = markers[1], label = Labels_Graph[1], alpha = 0.5)
        else:
            G = ax2.scatter(Data_x[(i-1)+nrows*ncols], Data_y[(i-1)+nrows*ncols], color = colours[1], marker = markers[1], label = Labels_Graph[1], alpha = 0.5)
        #ax2.errorbar(np.mean(Data_x[(i-1)+nrows*ncols]), np.mean(Data_y[(i-1)+nrows*ncols]), xerr = np.std(Data_x[(i-1)+nrows+ncols]), yerr = np.std(Data_y[(i-1)+nrows*ncols]), alpha = 0.9, fmt = '^', color = colours[1])
    
    #LinearFit = np.polyfit(Data_x[(i-1)+ncols*nrows], Data_y[(i-1)+ncols*nrows], 1)
        #print "The linear fit y = a*x + b had coefficients a = %.4e and b = %.4e" %(LinearFit[0], LinearFit[1])
        
        LinearFit = np.polyfit(Data_x[(i-1)], Data_y[(i-1)], 1)
        print "The linear fit y = a*x + b had coefficients a = %.4e and b = %.4e" %(LinearFit[0], LinearFit[1])

    
        ax.set_xticks(np.arange(Limits[0], Limits[1] + 1, 1))
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f'))
        ax.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))
        
        ax.set_yticks(np.arange(Limits[2], Limits[3] + 0.5, 0.5))
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
        ax.yaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))

        ax2.set_xticks(np.arange(Limits[4], Limits[5] + 1, 1))
        ax2.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f'))
        ax2.xaxis.set_minor_locator(mtick.MultipleLocator(tick_interval))

        if i < 4:
            ax.tick_params(axis='both', labelbottom = 'off')
            ax2.xaxis.get_major_ticks()[-1].set_visible(False)
        else:
            ax.xaxis.tick_bottom()
        if i > 3:
            ax2.tick_params(axis='both', labeltop = 'off')
            ax.xaxis.get_major_ticks()[-1].set_visible(False)
        else:
            ax2.xaxis.tick_top()
        if i != 1 and i != 4:
            ax.tick_params(axis = 'y', labelleft = 'off')
        if i == 1:
            ax.yaxis.get_major_ticks()[0].set_visible(True)
            ax2.xaxis.get_major_ticks()[0].set_visible(True)
            leg = ax.legend(handles = [H, G], loc = 2, scatterpoints=1, labelspacing=0.0)
            leg.draw_frame(False)  # Don't want a box frame
            for t in leg.get_texts():  # Reduce the size of the text
                    t.set_fontsize(label_size-3)
        if i == 4:
            ax.yaxis.get_major_ticks()[-1].set_visible(False)
        if i == 3:
            ax2.xaxis.get_major_ticks()[-1].set_visible(True)
        if i == 6:
            ax.xaxis.get_major_ticks()[-1].set_visible(True)

        ax.set_xlim(Limits[0:2])
        ax2.set_xlim(Limits[4:6])
        ax.set_ylim(Limits[2:4])

        label = "z = %.2f" %(SnapList[i-1])
        ax.text(9.25, 54.75, label)

    plt.tight_layout()
    fig.text(0.5, 0.97, Labels_Axis[2], ha = 'center')
    fig.text(0.5, 0.01, Labels_Axis[0], ha = 'center')
    fig.text(0.001, 0.5, Labels_Axis[1], va = 'center', rotation = 'vertical')
    fig.subplots_adjust(hspace = 0.0, wspace = 0.0)
    
    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

#####

def Plot_Scatter_SixPanel(Data_x, Data_y, Difference, Multiple_Scatter, Labels_Graph, Limits, Labels_Axis, Legend_Location, SnapList, Output_Tag, Output_Format):
    
    dilute = 1000

    label_size = 10
    plt.rc('xtick', labelsize=label_size)
    plt.rc('ytick', labelsize=label_size)
    plt.rc('text', usetex=True)
    
    tick_interval = 0.25
    
    nrows = 2
    ncols = 3

    fig = plt.figure()
    for i in xrange(1, nrows*ncols + 1):
        ax = fig.add_subplot(nrows, ncols, i)
    
        if i > 7:
            marker = markers[1]
            colour = colours[1]
            label = Labels_Graph[1]
        else:
            marker = markers[0]
            colour = colours[0]
            label = Labels_Graph[0]
            
        if (len(Data_x[i-1]) > dilute):
            w = sample(range(len(Data_x[i-1])), dilute)
        else:
            w = np.arange(0, len(Data_x[i-1]))

        if Difference == 0:
            Datay = Data_y[i-1][w]
            G1 = ax.scatter(Data_x[(i-1)+nrows*ncols][w], Data_y[(i-1)+nrows*ncols][w], color = colours[1], marker = markers[1], label = Labels_Graph[1], alpha = 0.5)
            label = Labels_Graph[0]
            marker = markers[0]
            colour = colours[0]
        else:
            marker = markers[1]
            colour = colours[1  ]
            Datay = Data_y[(i-1)+nrows*ncols][w] - Data_y[i-1][w]
            if Multiple_Scatter > 1:
                label = Labels_Graph[1]
            else:
                label = ''
        H = ax.scatter(Data_x[i-1][w], Datay, color = colour, marker = marker, label = label, alpha = 0.5)

        if Multiple_Scatter > 1:
            if (len(Data_x[i-1]) > dilute):
                w = sample(range(len(Data_x[(i-1)+nrows*ncols*2])), dilute)
            else:
                w = np.arange(0, len(Data_x[(i-1)+nrows*ncols*2]))
            if Difference == 0:
                G2 = ax.scatter(Data_x[(i-1)+nrows*ncols*2][w], Data_y[(i-1)+nrows*ncols*2][w], color = colours[2], marker = markers[2], label = Labels_Graph[2], alpha = 0.5)
            else:
                G2 = ax.scatter(Data_x[(i-1)+nrows*ncols*2][w], Data_y[(i-1)+nrows*ncols*2][w] - Data_y[i-1][w], color = colours[2], marker = markers[2], label = Labels_Graph[2], alpha = 0.5)

        if Multiple_Scatter > 2:
            if (len(Data_x[i-1]) > dilute):
                w = sample(range(len(Data_x[(i-1)+nrows*ncols*3])), dilute)
            else:
                w = np.arange(0, len(Data_x[(i-1)+nrows*ncols*3]))
            if Difference == 0:
                G3 = ax.scatter(Data_x[(i-1)+nrows*ncols*3][w], Data_y[(i-1)+nrows*ncols*3][w], color = colours[3], marker = markers[3], label = Labels_Graph[3], alpha = 0.5)
            else:
                G3 = ax.scatter(Data_x[(i-1)+nrows*ncols*3][w], Data_y[(i-1)+nrows*ncols*3][w] - Data_y[i-1][w], color = colours[3], marker = markers[3], label = Labels_Graph[3], alpha = 0.5)

        ax.set_xticks(np.arange(Limits[0], Limits[1] + 1, 1))
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
            if (Multiple_Scatter == 1 and Difference == 0):
                leg = ax.legend(handles = [H, G1], loc = 2, scatterpoints=1, labelspacing=0.0)
            elif (Multiple_Scatter == 1 and Difference == 1):
                leg = ax.legend(handles = [H], loc = 2, scatterpoints=1, labelspacing=0.0)
            elif (Multiple_Scatter == 2 and Difference == 0):
                leg = ax.legend(handles = [H, G1, G2], loc = 2, scatterpoints=1, labelspacing=0.0)
            elif (Multiple_Scatter == 2 and Difference == 1):
                leg = ax.legend(handles = [H, G2], loc = 2, scatterpoints=1, labelspacing=0.0)
            elif (Multiple_Scatter == 3 and Difference == 0):
                leg = ax.legend(handles = [H, G1, G2, G3], loc = 2, scatterpoints=1, labelspacing=0.0)
            elif (Multiple_Scatter == 3 and Difference == 1):
                leg = ax.legend(handles = [H, G2, G3], loc = 2, scatterpoints=1, labelspacing=0.0)
            leg.draw_frame(False)  # Don't want a box frame
            for t in leg.get_texts():  # Reduce the size of the text
                t.set_fontsize(label_size-3)
        if i == 4:
            ax.yaxis.get_major_ticks()[-1].set_visible(False)

        if i == 6:
            ax.xaxis.get_major_ticks()[-1].set_visible(True)


        ax.set_xlim(Limits[0:2])
        ax.set_ylim(Limits[2:4])
        
        label = "z = %.2f" %(SnapList[i-1])
        if Difference == 0:
            ax.text(9.25, 54.75, label, fontsize = label_size)
        else:
            ax.text(9.25, 0.75, label, fontsize = label_size)
            ax.axhline(y = 0, ls = '--', color = 'k', linewidth = 1.0)

    plt.tight_layout()
    fig.text(0.5, 0.01, Labels_Axis[0], ha = 'center')
    fig.text(0.001, 0.5, Labels_Axis[1], va = 'center', rotation = 'vertical')
    fig.subplots_adjust(hspace = 0.0, wspace = 0.0)
    
    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

#####

def Plot_LinearRegressions(Data_x, Data_y, SnapList, Output_Tag, Output_Format):

    ax = plt.subplot(111)

    Gradient = []
    Intercept = []
    Mean = []
    
    for i in xrange(0, len(Data_x)):
        LinearFit = np.polyfit(Data_x[i], Data_y[i], 1)
        print "The linear fit y = a*x + b had coefficients a = %.4e and b = %.4e" %(LinearFit[0], LinearFit[1])
        Gradient.append(LinearFit[0])
        Intercept.append(LinearFit[1])
        Mean.append(np.mean(Data_y[i]))


    ax.scatter(SnapList, Gradient)

    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()


######

def Plot_HeatMap_SixPanel(Data_x, Data_y, Difference, Multiple_Scatter, Labels_Graph, Limits, Labels_Axis, Legend_Location, SnapList, Output_Tag, Output_Format):
    
    Nbins = 40
    dilute = 1000
    
    label_size = 10
    plt.rc('xtick', labelsize=label_size)
    plt.rc('ytick', labelsize=label_size)
    plt.rc('text', usetex=True)
    
    tick_interval = 0.25
    
    nrows = 2
    ncols = 3
    
    fig = plt.figure()
    for i in xrange(1, nrows*ncols + 1):
        ax = fig.add_subplot(nrows, ncols, i)
        
        if i > 7:
            marker = markers[1]
            colour = colours[1]
            label = Labels_Graph[1]
        else:
            marker = markers[0]
            colour = colours[0]
            label = Labels_Graph[0]
        
        if (len(Data_x[i-1]) > dilute):
            w = sample(range(len(Data_x[i-1])), dilute)
        else:
            w = np.arange(0, len(Data_x[i-1]))
        
        if Difference == 0:
            Datay = Data_y[i-1][w]
            heatmap, xedges, yedges = np.histogram2d(Data_x[(i-1)+nrows*ncols], Data_y[(i-1)+nrows*ncols], bins = Nbins)
            print "Xedges[0] = %.2f xedges[-1] = %.2f yedges[0] = %.2f yedges[-1] = %.2f" %(xedges[0], xedges[-1], yedges[0], yedges[-1])
            ax.imshow(heatmap.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], origin='lower', cmap = 'Blues')
            label = Labels_Graph[0]
            marker = markers[0]
            colour = colours[0]
        else:
            marker = markers[1]
            colour = colours[1 ]
            Datay = Data_y[(i-1)+nrows*ncols][w] - Data_y[i-1][w]
            if Multiple_Scatter > 1:
                label = Labels_Graph[1]
            else:
                label = ''
        H = ax.scatter(Data_x[i-1][w], Datay, color = colour, marker = marker, label = label, alpha = 0.5)
        
        if Multiple_Scatter > 1:
            if (len(Data_x[i-1]) > dilute):
                w = sample(range(len(Data_x[(i-1)+nrows*ncols*2])), dilute)
            else:
                w = np.arange(0, len(Data_x[(i-1)+nrows*ncols*2]))
            if Difference == 0:
                G2 = ax.scatter(Data_x[(i-1)+nrows*ncols*2][w], Data_y[(i-1)+nrows*ncols*2][w], color = colours[2], marker = markers[2], label = Labels_Graph[2], alpha = 0.5)
            else:
                G2 = ax.scatter(Data_x[(i-1)+nrows*ncols*2][w], Data_y[(i-1)+nrows*ncols*2][w] - Data_y[i-1][w], color = colours[2], marker = markers[2], label = Labels_Graph[2], alpha = 0.5)

        if Multiple_Scatter > 2:
            if (len(Data_x[i-1]) > dilute):
                w = sample(range(len(Data_x[(i-1)+nrows*ncols*3])), dilute)
            else:
                w = np.arange(0, len(Data_x[(i-1)+nrows*ncols*3]))
            if Difference == 0:
                G3 = ax.scatter(Data_x[(i-1)+nrows*ncols*3][w], Data_y[(i-1)+nrows*ncols*3][w], color = colours[3], marker = markers[3], label = Labels_Graph[3], alpha = 0.5)
            else:
                G3 = ax.scatter(Data_x[(i-1)+nrows*ncols*3][w], Data_y[(i-1)+nrows*ncols*3][w] - Data_y[i-1][w], color = colours[3], marker = markers[3], label = Labels_Graph[3], alpha = 0.5)
                    
        ax.set_xticks(np.arange(Limits[0], Limits[1] + 1, 1))
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
        if (Multiple_Scatter == 1 and Difference == 0):
            leg = ax.legend(handles = [H], loc = 2, scatterpoints=1, labelspacing=0.0)
        elif (Multiple_Scatter == 1 and Difference == 1):
            leg = ax.legend(handles = [H], loc = 2, scatterpoints=1, labelspacing=0.0)
        elif (Multiple_Scatter == 2 and Difference == 0):
            leg = ax.legend(handles = [H], loc = 2, scatterpoints=1, labelspacing=0.0)
        elif (Multiple_Scatter == 2 and Difference == 1):
            leg = ax.legend(handles = [H], loc = 2, scatterpoints=1, labelspacing=0.0)
        elif (Multiple_Scatter == 3 and Difference == 0):
            leg = ax.legend(handles = [H], loc = 2, scatterpoints=1, labelspacing=0.0)
        elif (Multiple_Scatter == 3 and Difference == 1):
            leg = ax.legend(handles = [H], loc = 2, scatterpoints=1, labelspacing=0.0)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize(label_size-3)
        if i == 4:
            ax.yaxis.get_major_ticks()[-1].set_visible(False)
        
        if i == 6:
            ax.xaxis.get_major_ticks()[-1].set_visible(True)


        ax.set_xlim(Limits[0:2])
        ax.set_ylim(Limits[2:4])
        
        label = "z = %.2f" %(SnapList[i-1])
        if Difference == 0:
            ax.text(9.25, 54.75, label)
        else:
            ax.text(9.25, 0.75, label)
            ax.axhline(y = 0, ls = '--', color = 'k', linewidth = 1.0)

    plt.tight_layout()
    fig.text(0.5, 0.01, Labels_Axis[0], ha = 'center')
    fig.text(0.001, 0.5, Labels_Axis[1], va = 'center', rotation = 'vertical')
    fig.subplots_adjust(hspace = 0.0, wspace = 0.0)
    
    outputFile = './' + Output_Tag + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

#####
