#!/usr/bin/env python
from __future__ import print_function

import sys
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import pylab as plt

from mpi4py import MPI

import PlotScripts
import ReadScripts
import AllVars

comm= MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

output_format = ".png"
matplotlib.rcdefaults()

def calculate_volume_frac(ionized_cells):

    assert(len(np.argwhere(np.isnan(ionized_cells))) == 0)
           
    HII = np.mean(ionized_cells) 
    print("There is {0:.4e} fraction of HII.".format(HII))

    return HII

def hoshen_kopelman(ionized_cells, GridSize):

	print("Running the Hoshen-Kopelman Algorithm")

	## Just a quick function that replicates the !! behaviour in C.
	## If the input value (a) is != 0 it returns 1 otherwise it returns 0.
	def double_not(a):  
		if(a != 0):
			return 1
		else:
			return 0

	l = len(ionized_cells)
	m = len(ionized_cells)
	n = len(ionized_cells)
	

	max_labels = int(l*m*n / 2) # The maximum number of discrete ionized regions.
	labels = np.zeros((max_labels), dtype = np.int)

	test = np.zeros((l, m, n), dtype = np.int32)
	for i in range(0, l):
		for j in range(0, m):
			for k in range(0, n):
				if (ionized_cells[i][j][k] > 0.8):
					test[i][j][k] = 1
				else:
					test[i][j][k] = 0
	
	def make_set(label_number):
		label_number += 1	
		assert(label_number < max_labels), "The current_label value is greater than the maximum label limit."
		labels[label_number] = label_number 
		return label_number 

	def find(x):
		y = x
		while (labels[y] != y):
			y = labels[y]

		while(labels[x] != x):
			z = labels[x]
			labels[x] = y
			x = z

		return y

	def union(x, y):
		labels[find(x)] = find(y)
		return find(y)

	def union_3D(x, y, z):
		labels[find(x)] = find(y)
		labels[find(z)] = find(y)
		return find(y) 

	current_label = 0

	for i in range(0, l):
		for j in range(0, m):
			for k in range(0, l):

				if(test[i][j][k] == 1):

					if(i == 0):
						up = 0
					else:
						up = test[i-1][j][k]

					if(j == 0):
						left = 0
					else:
						left = test[i][j-1][k]

					if(k == 0):
						back = 0
					else:
						back = test[i][j][k-1]

					tmp = double_not(left) + double_not(up) + double_not(back) # Since the labels can be greater than 1, double_not returns 1 if the value is >= otherwise it returns 0.
					# So it will return 1 if there is a cluster at the position.

					if(tmp == 0): # We have a new cluster
						test[i][j][k] = make_set(current_label)
						current_label += 1
	
					elif(tmp == 1): # Part of an existing cluster
						test[i][j][k] = max(up, left, back)

					elif(tmp == 2): # Joins two existing clusters together
						if(up == 0):
							test[i][j][k] = union(left, back)
						elif(left == 0):
							test[i][j][k] = union(up, back)
						elif(back == 0):
							test[i][j][k] = union(up, left)
					elif(tmp == 3): # Joins three existing clusters together
						test[i][j][k] = union_3D(up, left, back)

	new_labels = np.zeros((max_labels), dtype = np.int)
	size_bubbles = np.zeros((max_labels), dtype = np.int32)
	new_labels_len = 0

	print("Finished labelling all the cells. Starting the relabelling.")

	for i in range(0, l):
		for j in range(0, m):
			for k in range(0, n):
				if(test[i][j][k] > 0):
					x = find(test[i][j][k])
					if(new_labels[x] == 0):
						new_labels_len += 1
						new_labels[x] = new_labels_len

					test[i][j][k] = new_labels[x]
					size_bubbles[test[i][j][k]] += 1

	total_clusters = new_labels_len		

	print("There was {0} clusters found".format(total_clusters))

	largest_bubble = max(size_bubbles)
	ionization_fraction = calculate_volume_frac(ionized_cells)
	num_cells_ionized = GridSize**3 * ionization_fraction 

	print("The largest bubble contains {0} cells. This is {1:.4f} of the entire box.  Ionization Fraction = {2:.4f}. The number of cells ionized is {3}. Num Cells Inside Bubble/Num Cells Ionized = {4:.4f}".format(largest_bubble, float(largest_bubble)/float(l*m*n), ionization_fraction, num_cells_ionized, float(largest_bubble)/float(num_cells_ionized))) 

	return float(largest_bubble)/float(num_cells_ionized)

##


def plot_hoshen(order_parameter, HI_frac_model, model_tags, OutputDir, output_tag):

    ax1 = plt.subplot(111)

    print("Plotting the Hoshen-Koppelman results.")

    print(HI_frac_model.shape)
    print(order_parameter.shape)
    for p in range(0, len(HI_frac_model)):
        print(HI_frac_model[p])
        print(order_parameter[p])
        ax1.plot(HI_frac_model[p], order_parameter[p], lw = 3, 
                 color = PlotScripts.colors[p], label = model_tags[p])

    ax1.set_xlabel(r"$\mathbf{\langle \chi_{HII}\rangle}$", 
                   size = PlotScripts.global_labelsize) 
    ax1.set_ylabel(r"$P_{inf}/\chi_\mathrm{HII}$", 
                   size = PlotScripts.global_labelsize)      

    leg = ax1.legend(loc='bottom right', numpoints=1,
             labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize-5)
    plt.tight_layout()

    outputFile = OutputDir + output_tag + output_format 
    plt.savefig(outputFile)  # Save the figure
    print('Saved file to {0}'.format(outputFile))		
    plt.close()

##

if __name__ == '__main__':

    ###########################   

    PlotScripts.Set_Params_Plot()      
 
    model_tags = [r"$\mathbf{f_{esc} = 0.35}$",
                  r"$\mathbf{f_{esc} \: \propto \: M_H^{-1}}$",
                  r"$\mathbf{f_{esc} \: \propto \: M_H}$",
                  r"$\mathbf{f_{esc} \: \propto \: f_ej}$",
                  r"$\mathbf{f_{esc} \: \propto \: SFR}$"]

    number_models = 5 

    model = 'new_paper'

    GridSize_model1 = 256
    precision_model1 = 2

    filepath_model1="/fred/oz004/jseiler/kali/self_consistent_output/constant/grids/cifog/newphoton_SF0.03_XHII"
    filepath_model2="/fred/oz004/jseiler/kali/self_consistent_output/anne/grids/cifog/1e8_1e12_0.99_0.10_XHII"
    filepath_model3="/fred/oz004/jseiler/kali/self_consistent_output/anne/grids/cifog/1e8_1e12_0.01_0.50_XHII"
    filepath_model4="/fred/oz004/jseiler/kali/self_consistent_output/fej/grids/cifog/newphoton_SF0.03_fej0.7_XHII"
    filepath_model5="/fred/oz004/jseiler/kali/self_consistent_output/SFR/grids/cifog/SFR_0.20_0.30_XHII"

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

    SnapList = np.arange(28, 98)

    XHII_array = np.zeros((number_models, len(SnapList)), dtype = np.float64)
    hoshen_array = np.zeros((number_models, len(SnapList)), dtype = np.float64)

    have_data = 1

    for snapnum in range(rank, len(SnapList), size):
        for model_number in range(len(fname_ionized)):
            XHII_fname = "{0}_{1:03d}".format(fname_ionized[model_number], 
                                              SnapList[snapnum])

            XHII = ReadScripts.read_binary_grid(XHII_fname,
                                                GridSize[model_number],
                                                precision[model_number])

            XHII_array[model_number][snapnum] = calculate_volume_frac(XHII)

            if have_data == 0:
                hoshen_array[model_number][snapnum] = hoshen_kopelman(XHII,
                                                                      GridSize[model_number])


    if have_data == 1:
        hoshen_array = np.loadtxt("./slurm_scripts/hoshen_data.txt")
        plot_hoshen(hoshen_array, XHII_array, model_tags, "./", "hoshen")
        exit()

    if rank == 0:
        hoshen_array_master = np.zeros_like(hoshen_array)
    else:
        hoshen_array_master = None

    comm.Reduce([hoshen_array, MPI.DOUBLE], 
                [hoshen_array_master, MPI.DOUBLE], 
                op = MPI.SUM, root = 0) # Sum the arrays across processors.

    if rank == 0:
        np.savetxt("./hoshen_data.txt", hoshen_array_master) 
