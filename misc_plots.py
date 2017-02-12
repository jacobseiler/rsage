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

Output_Format = '.png'

def plot_filtervalues(z, Mfilt):


    AllVars.Set_Params_Mysim()
    snaplist = np.arange(24, 78)

    ZZ = np.empty(len(snaplist))

    z_Mfilt = np.empty(len(snaplist))    
    z_Mfilt_sd = np.empty(len(snaplist))
    z_Mfilt_error = np.empty(len(snaplist))

    z_Mvir_reion = np.empty(len(snaplist))    
    z_Mvir_reion_sd = np.empty(len(snaplist))
    z_Mvir_reion_error = np.empty(len(snaplist))

    z_Mvir = np.empty(len(snaplist))    
    z_Mvir_sd = np.empty(len(snaplist))
    z_Mvir_error = np.empty(len(snaplist))

    for i in xrange(0, len(snaplist)):
	ZZ[i] = AllVars.SnapZ[snaplist[i]]
	
	w = np.where((z > (ZZ[i] - 0.01)) & (z < (ZZ[i] + 0.01)))[0]
	
	if (len(Mfilt[w]) > 0):

		z_Mfilt[i] = np.log10(np.mean(Mfilt[w])) 
		z_Mfilt_sd[i] = np.log10(np.std(Mfilt[w]))
		z_Mfilt_error[i] = 0.434*np.std(Mfilt[w])/np.mean(Mfilt[w])

		#z_Mvir_reion[i] = np.log10(np.mean(Mvir[w]))
		#z_Mvir_reion_sd[i] = np.log10(np.std(Mvir[w]))
		#z_Mvir_reion_error[i] = 0.434*np.std(Mvir[w])/np.mean(Mvir[w])

	else:
		z_Mfilt[i] = -1
		z_Mfilt_sd[i] = 0.001
		z_Mfilt_error[i] = 0.001

		#z_Mvir_reion[i] = -1
		#z_Mvir_reion_sd[i] = 0.001
		#z_Mvir_reion_error[i] = 0.001

	#print "z = %.4f \t Total = %.4e \t Count = %.4e \t z_Mfilt = %.4e \t z_Mfilt_sd = %.4e \t z_Mvir = %.4e \t z_Mvir_sd = %.4e" %(ZZ[i], np.sum(Mfilt[w]), len(Mfilt[w]), z_Mfilt[i], z_Mfilt_sd[i], z_Mvir_reion[i], z_Mvir_reion_sd[i])
	print "z = %.4f \t Total = %.4e \t Count = %.4e \t z_Mfilt = %.4e \t z_Mfilt_sd = %.4e" %(ZZ[i], np.sum(Mfilt[w]), len(Mfilt[w]), z_Mfilt[i], z_Mfilt_sd[i])

    fig = plt.figure()
    ax = plt.subplot(111)
          
    #ax.errorbar(ZZ[1:], z_Mfilt[1:], z_Mfilt_error[1:], color = 'r', fmt = '^', alpha = 0.7)
    ax.plot(ZZ[1:], z_Mfilt[1:], color = 'r', linestyle = '-', label = r'$M_\mathrm{filt}$')
    #ax.plot(ZZ[1:], z_Mvir_reion[1:], color = 'b', linestyle = '-', label = r'$M_\mathrm{vir}$ Ionized Regions')

    ax.fill_between(ZZ[1:], z_Mfilt[1:]-z_Mfilt_error[1:], z_Mfilt[1:]+z_Mfilt_error[1:], alpha = 0.5, color = 'r')
    #ax.fill_between(ZZ[1:], z_Mvir_reion[1:]-z_Mvir_reion_error[1:], z_Mvir_reion[1:]+z_Mvir_reion_error[1:], alpha = 0.5, color = 'b')

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


    outputFile = './z_Mfilt_1024_DRAGONS' + Output_Format
    plt.savefig(outputFile)  # Save the figure
    print 'Saved file to', outputFile
    plt.close()

if __name__ == '__main__':

    #fname_filtervalues = "/home/jseiler/SAGE-stuff/post_processed_SAGE/reionization_filter_1024_DRAGONS.dat" 
    #fd = open(fname_filtervalues, 'rb')
    #z, Mfilt, Mvir, Mod= np.loadtxt(fd, dtype = np.float32, unpack=True)  
    #fd.close()


    fname_filtervalues = "/home/jseiler/SAGE-stuff/post_processed_SAGE/whateven.dat"
    fd = open(fname_filtervalues, 'rb')
    z, Mfilt = np.loadtxt(fd, dtype = np.float32, unpack=True)  
    fd.close()

    print Mfilt

    #plot_filtervalues(z, Mfilt, Mvir, Mod)
    plot_filtervalues(z, Mfilt)

