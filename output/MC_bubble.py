#!/usr/bin/env python
from __future__ import print_function

import sys
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import pylab as plt
import time
import random
import itertools
import matplotlib.ticker as mtick

from mpi4py import MPI

import PlotScripts
import ReadScripts
import AllVars

comm= MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

output_format = ".png"
matplotlib.rcdefaults()

plt.rc('text', usetex=True)


def calculate_bubble_MC(z, ionized_cells, Ncell, OutputDir, output_tag):
# Calculate the size of ionized/neutral bubbles by choosing an ionized/neutral cell and measuring the distance to a cell of the opposite phase.

## Input ##
# z is the redshift we are doing the MC walk at.
# ionized_cells is the array that contains the ionization state of a cell.

## Output ##
# The output file will be of the form '*OutputDir*/MC/*output_tag*_*z*.dat'

    def MC_walk(cells, indices, phase, N, Ncell, outfile):
        # Function for the MC walk to calculate bubble size.
    
        ## Input ##
        # cells is the array that contains the ionization field (0 or 1).
        # indices is a length 3 array containing the x,y,z indices of the cells that are neutral (phase = 0) or ionized (phase = 1)
        # phase is 0 if we are starting with neutral or 1 if we are starting with ionized cell.
        # N is the number of random points we select.
    
        ## Output ##
        # A .dat file containing the results of the MC walk (number of steps until it encounters a phase transition).

        radii = []

        for j in range(0,int(N)):
                        
            if (j%20000 == 0):
                print("On walk number {0}".format(j))
    
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
                phase_transition = cells[(walk_x + R*x) % Ncell, (walk_y + R*y) % Ncell, (walk_z + R*z) % Ncell]
                if (phase_transition > 0.8):
                    phase_transition = 1
                else:
                    phase_transition = 0	
                if (R >= Ncell): # If the radius has gone beyond the number of available cells, 
                    phase_transition = (phase + 1) % 2 # We force the change.

            radii.append(R)
			
        np.savetxt(outfile, radii, delimiter = ',')
        print("MC file saved as {0}".format(outfile))

    print("")
    print("Calculating bubble size using MC walk.")	

    N = 1e5

    start_time = time.time()

    ionized_indices= np.where(ionized_cells > 0.8) # Calculate the array indices corresponding to neutral/ionized cells.
    unionized_indices = np.where(ionized_cells < 0.2)

    MCDir = OutputDir + 'MC/' 	
    if not os.path.exists(MCDir):
        os.makedirs(MCDir)
    outfile = MCDir + output_tag + '_z_%.3f.dat' %(z)
    MC_walk(ionized_cells, ionized_indices, 1, N, Ncell, outfile)

    print("MC took {0} seconds.".format(time.time() - start_time))


##########

def plot_bubble_MC(ZZ, fractions_HI, model_tags, file_tags, Ncell, OutputDir, output_tag):
# Plot the results from the MC walk.  Note that input_tag and labels are arrays that contains the input file tag and label for different realizations.

## Input ##
# ZZ is the redshift range we are plotting the results over.
# upperZ is the upper bound we are plotting to.
# input_tag is an array with the file tags for the input MC file for different realizations #### NOTE: This should be the same as *output_tag* in calculate_bubble_MC
# labels is an array that contains the graph labels for different realizations.

## Output ##
# The output file will be of the form '*OutputDir*/*output_tag*.*output_format*'
       
    def adjust_plot(ax, HI_labels, model_tags):
        """
        Adds text and adjusts the ranges for the MC plot. 

        Parameters
        ---------
        ax: Matplotlib axis handle. Required.
           Axis we're adjusting. 
    
        Returns
        -------
        No returns.
        The axis handle is opened by the outer function and passed in here.
        """

        for i in range(5):

            ax[i].set_xlabel(r'$\mathbf{R \: [h^{-1} Mpc]}$', 
                             size = PlotScripts.global_labelsize)
            ax[i].set_xscale('log')
            ax[i].set_xlim([2, 101]) 

            HI_string = "{0:.2f}".format(HI_labels[i])

            label = r"$\mathbf{\langle \chi_{HI}\rangle = " +HI_string + r"}$"
            ax[i].text(0.05, 0.9, label, transform = ax[i].transAxes,
                       fontsize = PlotScripts.global_fontsize)

            ax[i].tick_params(which = 'both', direction='in',
                              width = PlotScripts.global_tickwidth)

            ax[i].tick_params(which = 'major',
                              length = PlotScripts.global_ticklength)

            ax[i].tick_params(which = 'minor',
                              length = PlotScripts.global_ticklength - 2.5)

            for axis in ['top','bottom','left','right']:
                ax[i].spines[axis].set_linewidth(PlotScripts.global_axiswidth)


            tick_locs = np.arange(-1, 4, 1)
            ax[i].set_xticklabels([r"$\mathbf{10^{%d}}$" % x for x in tick_locs],
                                  fontsize = PlotScripts.global_fontsize)

        ax[0].set_ylabel(r'$\mathbf{R \: \: dP/dR}$', 
                         size = PlotScripts.global_labelsize)

        #ax[0].set_yscale('log', nonposy='clip')
        ax[0].set_ylim([0.0, 0.8])

        ax[0].yaxis.set_minor_locator(mtick.MultipleLocator(0.1))	
        tick_locs = np.arange(0.0, 0.9, 0.2)
        ax[0].set_yticklabels([r"$\mathbf{%.1f}$" % x for x in tick_locs],
                              fontsize = PlotScripts.global_fontsize)


        '''
        labels = ax[0].xaxis.get_ticklabels()
        locs = ax[0].xaxis.get_ticklocs()
        for label, loc in zip(labels, locs):
            print("{0} {1}".format(label, loc)) 

        labels = ax[0].yaxis.get_ticklabels()
        locs = ax[0].yaxis.get_ticklocs()
        for label, loc in zip(labels, locs):
            print("{0} {1}".format(label, loc)) 
        '''

    print("Plotting results from MC walk.")

    MCDir = OutputDir + 'MC/'


    fig, ax = plt.subplots(nrows=1, ncols=5, sharey='row', figsize=(16, 6))

    for model_number in range(len(ZZ)):
        for fraction in range(len(ZZ[model_number])):


            infile = "{0}{1}_z_{2:.3f}.dat".format(MCDir,
                                                   file_tags[model_number],
                                                   ZZ[model_number][fraction])

            if (os.path.exists(infile) == False):
                print("Could not find file {0}.  Skipping and moving on".format(infile))
                exit() 
            fd = open(infile, 'rb')

            print("Plotting Bubble size of file {0}".format(infile))
            R = np.loadtxt(fd)
            print("Maximum radius before scaling is {0} cells.".format(max(R)))
            #print("The ratio is: BoxSize = {0:.3f} Mpc/h, Ncell = {1} => 1cell = {2:.3f} Mpc/h".format(AllVars.BoxSize, Ncell[q], AllVars.BoxSize/float(Ncell[q])))		
            R *= AllVars.BoxSize/float(Ncell[model_number])
            print("Maximum radius after scaling is {0:.4f} Mpc/h.".format(max(R)))

            binwidth = 6*AllVars.BoxSize/float(Ncell[model_number])

            R_low = binwidth/2 
            if max(R) > AllVars.BoxSize/2:
                R_high = AllVars.BoxSize + binwidth 
            else:
                R_high = max(R) + binwidth 

            NB = int(np.ceil((R_high - R_low) / binwidth))

            (counts, binedges) = np.histogram(R, range=(R_low, R_high), bins=NB, density = True)

            # Set the x-axis values to be the centre of the bins
            xaxeshisto = binedges[:-1] + 0.5 * binwidth
            counts *= xaxeshisto
                       
            label = model_tags[model_number]

            ax[fraction].plot(xaxeshisto, counts, 
                              label = model_tags[model_number], 
                              color = PlotScripts.colors[model_number], 
                              ls = PlotScripts.linestyles[model_number]) 

    adjust_plot(ax, HI_fraction_target, model_tags)

    leg = ax[1].legend(loc='lower right', numpoints=1,
                       labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize-2)

    plt.tight_layout()
    plt.subplots_adjust(wspace = 0.0, hspace = 0.0)
    
    outputFile = OutputDir + output_tag + output_format 
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))				
    plt.close()


if __name__ == '__main__':

    ###########################   

    PlotScripts.Set_Params_Plot()      
 
    model_tags = [r"$\mathbf{f_{esc} = 0.35}$",
                  r"$\mathbf{f_{esc} \: \propto \: M_H^{-1}}$",
                  r"$\mathbf{f_{esc} \: \propto \: M_H}$",
                  r"$\mathbf{f_{esc} \: \propto \: f_ej}$",
                  r"$\mathbf{f_\mathrm{esc} = 0.35, photHI2}$"] 
                  #r"$\mathbf{f_{esc} \: \propto \: SFR}$"]

    output_tags = ["Constant",
                   "MH_Neg",
                   "MH_Pos",
                   "Ejected",
                   "SFR"]

    number_models = 5 

    model = 'new_paper_MC'

    OutputDir = "./ionization_plots/" + model + '/'
    if not os.path.exists(OutputDir):
        os.makedirs(OutputDir)

    GridSize_model1 = 256
    precision_model1 = 2

    filepath_model1="/fred/oz004/jseiler/kali/self_consistent_output/constant/grids/cifog/newphoton_SF0.03_XHII"
    filepath_model2="/fred/oz004/jseiler/kali/self_consistent_output/anne/grids/cifog/1e8_1e12_0.99_0.10_XHII"
    filepath_model3="/fred/oz004/jseiler/kali/self_consistent_output/anne/grids/cifog/1e8_1e12_0.01_0.50_XHII"
    filepath_model4="/fred/oz004/jseiler/kali/self_consistent_output/fej/grids/cifog/newphoton_SF0.03_fej0.7_XHII"
    #filepath_model5="/fred/oz004/jseiler/kali/self_consistent_output/SFR/grids/cifog/SFR_0.20_0.30_XHII"
    filepath_model5="/fred/oz004/jseiler/kali/self_consistent_output/constant/grids/cifog/const0.35_photHI2_XHII"

    precision = [precision_model1, 
                 precision_model1, 
                 precision_model1, 
                 precision_model1, 
                 precision_model1]

    GridSize = [GridSize_model1, 
                GridSize_model1, 
                GridSize_model1, 
                GridSize_model1, 
                GridSize_model1]

    fname_ionized = [filepath_model1,
                     filepath_model2,
                     filepath_model3,
                     filepath_model4,
                     filepath_model5]

    cosmo = AllVars.Set_Params_Kali() #  Let's just assume we're using Kali.


    SnapList = [[19, 26, 34, 42, 48],
                [18, 25, 34, 43, 49],
                [22, 29, 37, 45, 52],
                [17, 24, 33, 41, 48],
                #[21, 28, 36, 44, 50]]
                [19, 26, 34, 42, 48]]

    Redshift = np.zeros_like(SnapList, dtype=np.float32)
    HI_fraction_target = [0.90, 0.75, 0.50, 0.25, 0.10] 

    for snap in range(len(SnapList)):
        for inner_snap in range(len(SnapList[snap])):
            SnapList[snap][inner_snap] += 28    
            Redshift[snap][inner_snap] = AllVars.SnapZ[SnapList[snap][inner_snap]]

    have_data = 1

    if have_data == 1:
        plot_bubble_MC(Redshift, HI_fraction_target, model_tags, output_tags,
                       GridSize, OutputDir, "Bubbles_photHI2")
        exit()   
 

    for model_number in range(rank, len(fname_ionized), size):
        for snapnum in range(len(SnapList[model_number])):
            this_snapnum = SnapList[model_number][snapnum]
            XHII_fname = "{0}_{1:03d}".format(fname_ionized[model_number], 
                                              this_snapnum) 

            XHII = ReadScripts.read_binary_grid(XHII_fname,
                                                GridSize[model_number],
                                                precision[model_number])

            calculate_bubble_MC(AllVars.SnapZ[this_snapnum], XHII, 
                                GridSize[model_number], OutputDir, 
                                output_tags[model_number])



