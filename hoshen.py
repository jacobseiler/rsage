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


from astropy import units as u
from astropy import cosmology
import matplotlib.ticker as mtick
import PlotScripts
import ReadScripts
import AllVars

label_size = 20 # 12 for paper.


def hoshen_kopelman(ionized_cells, ionization_fraction):

	## Just a quick function that replicates the !! behaviour in C.
	## If the input value (a) is != 0 it returns 1 otherwise it returns 0.

	print "Running Hoshen-Kopelman Algorithm for an ionization fraction %.3f." %(ionization_fraction)

	def double_not(a):  
		if(a != 0):
			return 1
		else:
			return 0


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

	l = len(ionized_cells)
	m = len(ionized_cells)
	n = len(ionized_cells)
	

	max_labels = l*m*n / 2 # The maximum number of discrete ionized regions.
	labels = np.zeros((max_labels), dtype = np.int32)

	test = np.zeros((l, m, n), dtype = np.int32)
	for i in xrange(0, l):
		for j in xrange(0, m):
			for k in xrange(0, n):
				if (ionized_cells[i][j][k] < ionization_fraction):
					test[i][j][k] = 1
				else:
					test[i][j][k] = 0

	cells_ionized = sum(test) # An ionized cell has value of 1 so the number of ionized cells is just the sum of test.	

	current_label = 0

	for i in xrange(0, l):
		for j in xrange(0, m):
			for k in xrange(0, l):

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

	print "Finished labelling all the cells, doing the second sweep through."

	for i in xrange(0, l):
		for j in xrange(0, m):
			for k in xrange(0, n):
				if(test[i][j][k] > 0):
					x = find(test[i][j][k])
					if(new_labels[x] == 0):
						new_labels_len += 1
						new_labels[x] = new_labels_len

					test[i][j][k] = new_labels[x]
					size_bubbles[test[i][j][k]] += 1

	total_clusters = new_labels_len		

	print "There was %d clusters found" %(total_clusters)

	largest_bubble = max(size_bubbles)

	print "The largest bubble contains %d cells. This is %.4f of the entire box.  Ionization Fraction = %.4f. The number of cells ionized is %d. Num Cells Inside Bubble/Num Cells Ionized = %.4f" %(largest_bubble, float(largest_bubble)/float(l*m*n), ionization_fraction, cells_ionized, float(largest_bubble)/float(cells_ionized)) 

	return float(largest_bubble)/float(cells_ionized)

##

def run_algorithm(ionization_fraction):
	matrix = np.zeros((128, 128, 128), dtype = np.float32)

	for i in xrange(0, len(matrix)):
		for j in xrange(0, len(matrix)):
			for k in xrange(0, len(matrix)):
				matrix[i][j][k] = random.uniform(0, 1) 

	order_parameter = []	
	
	for xi in ionization_fraction: 
		order_parameter.append(hoshen_kopelman(matrix, xi))

	return order_parameter

##

def plot_order_parameter(ionization_fraction, order_parameter):

	ax1 = plt.subplot(111)
	
	ax1.plot(ionization_fraction, order_parameter, lw = 3, ls = '-')

	ax1.set_xlabel(r"$x_\mathrm{HII}$", size = label_size)
 	ax1.set_ylabel(r"$P_\mathrm{inf}/x_\mathrm{HII}$", size = label_size)      
 
	plt.tight_layout()

	outputFile = './hosh_test.png' 
	plt.savefig(outputFile)  # Save the figure
	print 'Saved file to', outputFile
	plt.close()

##

if __name__ == '__main__':

	ionization_fraction = np.arange(0.025, 1, 0.025)
	order_parameter = run_algorithm(ionization_fraction)
	#order_parameter = [8.617387973956339e-05, 8.63052713341836e-05, 0.00010859254480240948, 0.00016267747993445054, 0.0004856949453306359, 0.00576613083064141, 0.7136131125149394, 0.8894560383494995, 0.9503515610072987, 0.9763831118285013, 0.9887944314981405, 0.9947896898242197, 0.997792191391382, 0.9991435810856882, 0.9997078420539404, 0.9999194561680467, 0.9999831585660305, 0.9999978795218534, 1.0]

	plot_order_parameter(ionization_fraction, order_parameter)


