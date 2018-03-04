#!/usr/bin/env python
from __future__ import print_function

import sys

import numpy as np
from numpy import *
np.set_printoptions(threshold = np.nan, linewidth = 1000000)

import matplotlib
matplotlib.use('Agg')

import os
import pylab as plt
sys.path.append('/home/jseiler/SAGE-stuff/output/')

import PlotScripts
import ReadScripts
import AllVars

from optparse import OptionParser

from random import sample, seed
from os.path import getsize as getFileSize
import math
import random
import csv
from matplotlib.colors import LogNorm
import time

from matplotlib.ticker import MultipleLocator
from numpy.fft import fftn, ifftn
import matplotlib.patheffects as PathEffects
from astropy import units as u
from astropy import cosmology
import matplotlib.ticker as mtick
import itertools
from matplotlib import patches

from mpi4py import MPI

cut_slice = 40
output_format = ".png"

def create_redshift_grid(SnapList, SAGE_params, precision):

    GridSize = SAGE_params['GridSize'][0]
    xHII_base = SAGE_params['PhotoionDir'][0] + "XHII_" + SAGE_params['FileNameGalaxies'][0]
    redshift_output_base = SAGE_params['PhotoionDir'][0] + "ReionRedshift"

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

def read_redshift_grid(SnapList, SAGE_params, precision, output_tag, output_dir="./"):

    GridSize = SAGE_params['GridSize'][0]
    redshift_input_base = SAGE_params['PhotoionDir'][0] + "ReionRedshift"

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

    outputFile = "{0}{1}{2}".format(output_dir, output_tag, output_format) 
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))

    plt.close()

def create_SAGE_ini(SAGE_params, SAGE_params_names, increment_HighSnap, outputdir="./"):
    fname = "{0}test.ini".format(outputdir)

    if (increment_HighSnap == 1):
        SAGE_params['HighSnap'] = [SAGE_params['HighSnap'][0] + 1] 
        SAGE_params['ReionSnap'] = [SAGE_params['ReionSnap'][0] + 1] 
    
    with open (fname, "w+") as f:
        for name in SAGE_params_names:
            string = "{0} {1}\n".format(name, SAGE_params[name][0])
            f.write(string)
    print("Successfully wrote to {0}".format(fname))

if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option("-r", "--reionredshift", dest="reionredshift", help="Set to 1 to generate the reionization redshift grid for the snapshot range specified. Default: 0.", default = 0, type = int)
    parser.add_option("-n", "--snap_range", dest="snap_range", nargs = 2, help="Snapshot range of interest.  Range is inclusive.  Default: Range specified by the SAGE ini file (LowSnap to HighSnap inclusive).", default = (0, 0))
    parser.add_option("-f", "--SAGE_fname", dest="SAGE_fname", help="Location of the SAGE ini file. REQUIRED")
    parser.add_option("-p", "--precision", type = int, dest="precision", help="Precision of the grid files. 0 for int, 1 for float, 2 for double. Default: 2", default = 2) 
    parser.add_option("-g", "--galaxy_name", dest="galaxy_name", help="Overwrites the name of the galaxies specified within the SAGE ini file.  Default: Specified by SAGE ini file.") 
    parser.add_option("-d", "--galaxy_dir", dest="galaxy_dir", help="Overwrites the directory containing the galaxies specified within the SAGE ini file.  Default: Specified by SAGE ini file.") 
    parser.add_option("-c", "--create_SAGE", dest="create_SAGE_ini", help="Set to 1 to generate a SAGE ini file with HighSnap incremented by 1.  Default: 0", default = 0, type = int)

    (opt, args) = parser.parse_args()

    if (opt.SAGE_fname == None):
        parser.print_help()
        exit()

    if (opt.precision== None):
        parser.print_help()
        exit()

    AllVars.Set_Params_Kali()

    precision = opt.precision

    SAGE_params, SAGE_params_names = ReadScripts.read_SAGE_ini(opt.SAGE_fname)

    if (opt.galaxy_name != None):
        print("Overwriting the name of galaxies from the SAGE ini file.")
        SAGE_params['FileNameGalaxies'] = opt.galaxy_name     

    if (opt.galaxy_dir != None):
        print("Overwriting the directory contain the galaxies from the SAGE ini file.") 
        SAGE_params['OutputDir'] = opt.galaxy_name     

    if (opt.snap_range[0] == 0 and opt.snap_range[1] == 0):
        print("The snapshot range was not specified.  We are using the SnapLow and SnapHigh values specified in the SAGE ini file of ({0}, {1}).".format(SAGE_params['LowSnap'][0], SAGE_params['HighSnap'][0]))
        LowSnap = SAGE_params['LowSnap'][0]
        HighSnap = SAGE_params['HighSnap'][0]
    else:
        LowSnap = opt.snap_range[0]
        HighSnap = opt.snap_range[1]

    if (LowSnap > HighSnap):
        raise ValueError("SnapLow should be smaller (or equal) to SnapHigh")
    SnapList = np.arange(LowSnap, HighSnap) 
           
    if opt.reionredshift:     
        create_redshift_grid(SnapList, SAGE_params, precision)

    if opt.create_SAGE_ini:
        create_SAGE_ini(SAGE_params, SAGE_params_names, 1)
    #read_redshift_grid(SnapList, SAGE_params, precision, "test_reion_redshift") 
