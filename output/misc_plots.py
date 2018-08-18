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
plt.rc('text', usetex=True)

output_format = '.png'


if __name__ == '__main__':

    gal_filepath="/home/jseiler/self_consistent_SAGE/tests/test_output/galaxies/kali_test_z5.782"
    merged_gal_filepath="/home/jseiler/self_consistent_SAGE/tests/test_output/galaxies/kali_test_MergedGalaxies"
    snap = 77

    GG, Gal_Desc = ReadScripts.ReadGals_SAGE(gal_filepath, 0, 99) # Read galaxies 
    G_Merged, _ = ReadScripts.ReadGals_SAGE(merged_gal_filepath, 0, 99) 
    G = ReadScripts.Join_Arrays(GG, G_Merged, Gal_Desc) # Then join them together for all galaxies. 


    w = np.where((G.GridHistory[:,snap] != -1) & \
                 (G.GridStellarMass[:,snap] > 0.0))[0]

    w_wrong = w[np.where(G.GridNgamma_HI[w,snap] == 0)[0]]
    w_right = w[np.where(G.GridNgamma_HI[w,snap] > 0)[0]]

    print("There were {0} galaxies at snapshot {1}.  Of these, {2} had an "
          "Ngamma value of 0.".format(len(w), snap, len(w_wrong))) 

    no_sat = np.zeros(len(w_wrong))
    for i in range(99):
        w_sat = np.where(G.GridType[w_wrong,i] > 0)[0] 
        no_sat[w_sat] = 1
    w_nosat = np.where(no_sat == 0)[0]
   
    for my_idx in w_nosat: 
        idx = w_nosat[my_idx] 
        #print(max(G.GridStellarMass[w_wrong,snap]))
        #print(np.argmax(G.GridStellarMass[w_wrong,snap]))
                   
        fig = plt.figure(figsize = (8,8))
        ax = fig.add_subplot(111)
        print(G.TreeNr[w_wrong[idx]])
        check_snaps = np.arange(20, 99) 
        for snap in check_snaps: 
            print("Snap {0}\tColdGas {1:.4e}\tStellarMass {2:.4e}\t"
                  "Nion {3:.4e}\tColdCrit {4:.4e}\tGridType {5}\t"
                  "LenMergerGal {6}\tFoFMass {7:.4e}".format(snap,
                  G.GridColdGas[w_wrong[idx],snap],
                  G.GridStellarMass[w_wrong[idx],snap],
                  G.GridNgamma_HI[w_wrong[idx],snap],
                  G.ColdCrit[w_wrong[idx],snap],
                  G.GridType[w_wrong[idx],snap],
                  G.LenMergerGal[w_wrong[idx],snap],
                  G.GridFoFMass[w_wrong[idx],snap]))

            if (G.LenMergerGal[w_wrong[idx],snap] != -1):
                ax.axvline(snap, lw = 1, ls = '--', color = 'k')

        AllVars.Set_Params_Kali()
        PlotScripts.Set_Params_Plot()
        Mass = (G.GridStellarMass[w_wrong[idx], :] * 1.0e10 / AllVars.Hubble_h)
        ColdGas = (G.GridColdGas[w_wrong[idx], :] * 1.0e10 / AllVars.Hubble_h)
        HotGas = (G.GridHotGas[w_wrong[idx], :] * 1.0e10 / AllVars.Hubble_h)
        ColdCrit= (G.ColdCrit[w_wrong[idx], :] * 1.0e10 / AllVars.Hubble_h)
        FoFMass = G.GridFoFMass[w_wrong[idx], :] * 1.0e10 / AllVars.Hubble_h

        ax.plot(np.arange(0,99), Mass, color = 'k',
                label = "Stellar Mass")
        ax.plot(np.arange(0,99), HotGas, color = 'r',
                label = "Hot Gas")
        ax.plot(np.arange(0,99), ColdGas, color = 'b',
                label = "Cold Gas")
        ax.plot(np.arange(0,99), ColdCrit, color = 'b',
                ls = '--', label = "Cold Crit")
        ax.plot(np.arange(0,99), FoFMass, color = 'g',
                ls = '-', label = "FoF Mass")


        ax.set_xlabel("Snapshot Number", size = PlotScripts.global_labelsize)
        ax.set_ylabel("Mass [Msun]", size = PlotScripts.global_labelsize)
        ax.set_yscale('log')

        leg = ax.legend(loc='upper left', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize(PlotScripts.global_legendsize - 2)
        output_tag = "Masses"
        outputFile = "Mass/{0}_{2}{1}".format(output_tag, output_format, idx)
        plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
        print('Saved file to {0}'.format(outputFile))
        plt.close()


    no_sat = np.zeros(len(w_right))
    for i in range(99):
        w_sat = np.where(G.GridType[w_right,i] > 0)[0] 
        no_sat[w_sat] = 1
    w_nosat = np.where(no_sat == 0)[0]
   
    for my_idx in w_nosat: 
        idx = w_nosat[my_idx] 
        #print(max(G.GridStellarMass[w_right,snap]))
        #print(np.argmax(G.GridStellarMass[w_right,snap]))

        '''
        print("Stellar Mass")
        print(G.GridStellarMass[w_right[idx],:])

        print("Cold Gas")
        print(G.GridColdGas[w_right[idx],:])
        w_nogas = np.where(G.GridColdGas[w_right[idx],:] == 0)[0]

        print("Nion")
        print(G.GridNgamma_HI[w_right[idx],:])

        print("Snapshot where no cold gas is {0}".format(w_nogas))
        print("Stellar Mass at this snap")
        print(G.GridStellarMass[w_right[idx], w_nogas])
        print("Nion at this snap")
        print(G.GridNgamma_HI[w_right[idx], w_nogas])
        '''
      
        fig = plt.figure(figsize = (8,8))
        ax = fig.add_subplot(111)
        print(G.TreeNr[w_right[idx]])
        check_snaps = np.arange(20, 99) 
        for snap in check_snaps: 
            print("Snap {0}\tColdGas {1:.4e}\tStellarMass {2:.4e}\t"
                  "Nion {3:.4e}\tColdCrit {4:.4e}\tGridType {5}\t"
                  "LenMergerGal {6}\tFoFMass {7:.4e}".format(snap,
                  G.GridColdGas[w_right[idx],snap],
                  G.GridStellarMass[w_right[idx],snap],
                  G.GridNgamma_HI[w_right[idx],snap],
                  G.ColdCrit[w_right[idx],snap],
                  G.GridType[w_right[idx],snap],
                  G.LenMergerGal[w_right[idx],snap],
                  G.GridFoFMass[w_right[idx],snap]))

            if (G.LenMergerGal[w_right[idx],snap] != -1):
                ax.axvline(snap, lw = 1, ls = '--', color = 'k')

        AllVars.Set_Params_Kali()
        PlotScripts.Set_Params_Plot()
        Mass = (G.GridStellarMass[w_right[idx], :] * 1.0e10 / AllVars.Hubble_h)
        ColdGas = (G.GridColdGas[w_right[idx], :] * 1.0e10 / AllVars.Hubble_h)
        HotGas = (G.GridHotGas[w_right[idx], :] * 1.0e10 / AllVars.Hubble_h)
        ColdCrit= (G.ColdCrit[w_right[idx], :] * 1.0e10 / AllVars.Hubble_h)
        FoFMass = G.GridFoFMass[w_right[idx], :] * 1.0e10 / AllVars.Hubble_h

        ax.plot(np.arange(0,99), Mass, color = 'k',
                label = "Stellar Mass")
        ax.plot(np.arange(0,99), HotGas, color = 'r',
                label = "Hot Gas")
        ax.plot(np.arange(0,99), ColdGas, color = 'b',
                label = "Cold Gas")
        ax.plot(np.arange(0,99), ColdCrit, color = 'b',
                ls = '--', label = "Cold Crit")
        ax.plot(np.arange(0,99), FoFMass, color = 'g',
                ls = '-', label = "FoF Mass")


        ax.set_xlabel("Snapshot Number", size = PlotScripts.global_labelsize)
        ax.set_ylabel("Mass [Msun]", size = PlotScripts.global_labelsize)
        ax.set_yscale('log')

        leg = ax.legend(loc='upper left', numpoints=1, labelspacing=0.1)
        leg.draw_frame(False)  # Don't want a box frame
        for t in leg.get_texts():  # Reduce the size of the text
            t.set_fontsize(PlotScripts.global_legendsize - 2)
        output_tag = "Masses"
        outputFile = "Mass_correct/{0}_{2}{1}".format(output_tag, output_format, idx)
        plt.savefig(outputFile, bbox_inches='tight')  # Save the figure
        print('Saved file to {0}'.format(outputFile))
        plt.close()
