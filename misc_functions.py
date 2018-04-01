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
sys.path.append('/home/jseiler/self_consistent_SAGE/output/')

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

cut_slice = 40
output_format = ".png"

def create_redshift_grid(SnapList, SAGE_params, precision):

    AllVars.Set_Params_Kali()

    GridSize = SAGE_params["GridSize"][0]
    xHII_base = SAGE_params["PhotoionDir"][0] + "/XHII" 
    redshift_output_base = SAGE_params["PhotoionDir"][0] + "/" + SAGE_params["ReionRedshiftName"][0] 

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
        
    fname_out = "{0}".format(redshift_output_base) 
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

def increment_ini(SAGE_params, SAGE_params_names, cifog_params,
                  cifog_params_names, cifog_headers, increment_snaps, opt):

    if (increment_snaps == 1):
        SAGE_params['HighSnap'] = [SAGE_params['HighSnap'][0] + 1] 
        SAGE_params['ReionSnap'] = [SAGE_params['ReionSnap'][0] + 1] 

        cifog_params["stopSnapshot"] = [cifog_params["stopSnapshot"][0] + 1]
 
    fname_SAGE = "{0}/SAGE_snap{1}.ini".format(opt["ini_dir"],
                                               SAGE_params["HighSnap"][0])
    fname_cifog = "{0}/cifog_snap{1}.ini".format(opt["ini_dir"],
                                                 cifog_params["stopSnapshot"][0])
   
    with open (fname_SAGE, "w+") as f:
        for name in SAGE_params_names:
            string = "{0} {1}\n".format(name, SAGE_params[name][0])
            f.write(string)

    with open (fname_cifog, "w+") as f:
        for name in cifog_params_names:
            if name in cifog_headers:
                header_string = "{0}".format(cifog_headers[name])
                f.write(header_string)
            string = "{0} = {1}\n".format(name, cifog_params[name][0])
            f.write(string)

    print("Successfully wrote to {0}".format(fname_SAGE))
    print("Successfully wrote to {0}".format(fname_cifog))

def parse_input_arguments():

    parser = OptionParser()

    parser.add_option("-f", "--SAGE_ini", dest="SAGE_fname", help="Location of the SAGE ini file.  Required.")
    parser.add_option("-c", "--cifog_ini", dest="cifog_fname", help="Location of the cifog ini file.")
    parser.add_option("-n", "--snap_range", dest="snap_range", nargs = 2, help="Snapshot range of interest.  Range is inclusive.  Default: Range specified by the SAGE ini file (LowSnap to HighSnap inclusive).", default = (0, 0)) 
    parser.add_option("-p", "--precision", type = int, dest="precision", help="Precision of the grid files. 0 for int, 1 for float, 2 for double. Default: 2", default =2) 
    parser.add_option("-r", "--reionredshift", dest="reionredshift", help="Set to 1 to generate the reionization redshift grid for the snapshot range specified. Default: 0.", default = 0, type = int)
    parser.add_option("-g", "--galaxy_name", dest="galaxy_name", help="Overwrites the name of the galaxies specified within the SAGE ini file.  Default: Specified by SAGE ini file.") 
    parser.add_option("-d", "--galaxy_dir", dest="galaxy_dir", help="Overwrites the directory containing the galaxies specified within the SAGE ini file.  Default: Specified by SAGE ini file.") 
    parser.add_option("-i", "--ini_dir", dest="ini_dir", 
                      help="Specifies the directory to create two new .ini "
                      "files that are identical to the specified .ini files "
                      "but with value of HighSnap and stopSnapshot "
                      "incremented by 1.  Default: None", 
                      default = None)

    (opt, args) = parser.parse_args()

    if (opt.SAGE_fname == None):
        parser.print_help()
        exit()

    if (opt.precision == None):
        parser.print_help()
        exit()

    if (opt.ini_dir != None and opt.cifog_fname == None):
        print("To increment the ini files both a SAGE and cifog ini file must be specified.")
        parser.print_help() 
        exit()

    return vars(opt) 

if __name__ == '__main__':


    opt = parse_input_arguments()
    
    precision = opt["precision"]

    SAGE_params, SAGE_params_names = ReadScripts.read_SAGE_ini(opt["SAGE_fname"])

    if (opt["galaxy_name"] != None):
        print("Overwriting the name of galaxies from the SAGE ini file.")
        SAGE_params['FileNameGalaxies'] = opt["galaxy_name"] 

    if (opt["galaxy_dir"] != None):
        print("Overwriting the directory contain the galaxies from the SAGE ini file.") 
        SAGE_params['OutputDir'] = opt["galaxy_name"] 

    if (opt["snap_range"][0] == 0 and opt["snap_range"][1] == 0):
        print("The snapshot range was not specified.  We are using the SnapLow and SnapHigh values specified in the SAGE ini file of ({0}, {1}).".format(SAGE_params['LowSnap'][0], SAGE_params['HighSnap'][0]))
        LowSnap = SAGE_params['LowSnap'][0] + 1 
        HighSnap = SAGE_params['HighSnap'][0] + 1
    else:
        LowSnap = opt["snap_range"][0] 
        HighSnap = opt["snap_range"][1] 

    if (LowSnap > HighSnap):
        raise ValueError("SnapLow should be smaller (or equal) to SnapHigh")
    SnapList = np.arange(LowSnap, HighSnap + 1) 
           
    if opt["reionredshift"]:
        create_redshift_grid(SnapList, SAGE_params, precision)

    if opt["ini_dir"] != None:
        cifog_params, cifog_params_names, cifog_headers = ReadScripts.read_cifog_ini(opt["cifog_fname"])
        increment_ini(SAGE_params, SAGE_params_names, cifog_params,
                      cifog_params_names, cifog_headers, 1, opt)
    #read_redshift_grid(SnapList, SAGE_params, precision, "test_reion_redshift") 
