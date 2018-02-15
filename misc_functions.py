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

cut_slice = 40
output_format = ".png"

def create_redshift_grid(SnapList, GridSize, precision, xHII_base, redshift_output_base):

    reionization_redshift_grid = np.full((pow(GridSize, 3)), -1.0)

    for snapshot in SnapList:
        xHII_fname = "{0}_{1:03d}".format(xHII_base, snapshot)
        xHII_grid = ReadScripts.read_binary_grid(xHII_fname, GridSize, precision, False) 
        
        w_ionized_snap = np.where((xHII_grid > 0.9))[0] # Indices of the ionized cells.
        w_not_ionized = np.where((reionization_redshift_grid < -0.5))[0] # Indices of the cells that have yet to be ionized. 
    
        w_to_update = np.intersect1d(w_ionized_snap, w_not_ionized) 

        print("")
        print("For snapshot {0} there were {1} cells already ionized, {2} cells ionized in this snapshot resulting in {3} cells we are about to ionized".format(snapshot, pow(GridSize, 3) - len(w_not_ionized), len(w_ionized_snap), len(w_to_update))) 
        print("{0:.2f}% of cells will be ionized after this snapshot.".format((pow(GridSize, 3) - len(w_not_ionized) + len(w_to_update)) / pow(GridSize, 3) * 100.0))
        print("")
        reionization_redshift_grid[w_to_update] = AllVars.SnapZ[snapshot]
        
    fname_out = "{0}_{1}_{2}".format(redshift_output_base, SnapList[0], SnapList[-1])  
    reionization_redshift_grid.tofile(fname_out)

    print("Reionization redshift grid saved to file {0}".format(fname_out))

def read_redshift_grid(SnapList,GridSize, precision, redshift_input_base, output_tag):

    fname_in = "{0}_{1}_{2}".format(redshift_input_base, SnapList[0], SnapList[-1])      
    reionization_redshift_grid = ReadScripts.read_binary_grid(fname_in, GridSize, precision, True)

    print(reionization_redshift_grid[0,0,0:10])
    ax1 = plt.subplot(111)

    im = ax1.imshow(reionization_redshift_grid[:, :, cut_slice:cut_slice+1].mean(axis = -1), interpolation='bilinear', origin='low', extent = [0, AllVars.BoxSize, 0, AllVars.BoxSize], cmap = 'Dark2', vmin = AllVars.SnapZ[SnapList[0]], vmax = AllVars.SnapZ[SnapList[-1]])

    cbar = plt.colorbar(im, ax = ax1)
    cbar.set_label(r'$z_{\mathrm{reion}}$')

    ax1.set_xlabel(r'$\mathrm{x}  (h^{-1}Mpc)$')
    ax1.set_ylabel(r'$\mathrm{y}  (h^{-1}Mpc)$')

    ax1.set_xlim([0.0, AllVars.BoxSize])
    ax1.set_ylim([0.0, AllVars.BoxSize])

    outputFile = "./" + output_tag + output_format
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))

    plt.close()

if __name__ == '__main__':

    if (len(sys.argv) != 7):
        print("Usage: python3 self_consistent_SAGE.py <SnapLow> <SnapHigh> <Ionization Field Base Name> <GridSize> <Precision> <Redshift Ionization Grid Output Base Name>")
        exit()

    AllVars.Set_Params_Kali()

    SnapLow = int(sys.argv[1])
    SnapHigh = int(sys.argv[2])

    if (SnapLow > SnapHigh):
        raise ValueError("SnapLow should be smaller (or equal) to SnapHigh")
    SnapList = np.arange(SnapLow, SnapHigh + 1)    

    xHII_base = sys.argv[3]
    GridSize = int(sys.argv[4])
    precision = int(sys.argv[5])

    if (precision < 1 or precision > 2):
        raise ValueError("Precision should either be 1 (float) or 2 (double)") 

    redshift_output_base = sys.argv[6]

    #create_redshift_grid(SnapList, GridSize, precision, xHII_base, redshift_output_base)
    read_redshift_grid(SnapList, GridSize, precision, redshift_output_base, "test_reion_redshift") 
