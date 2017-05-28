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

colors = ['r', 'b', 'g', 'c', 'm', 'k']

matplotlib.rcdefaults()
plt.rc('axes', color_cycle=['k', 'b', 'r', 'g', 'm', 'c','y', '0.5'], labelsize='x-large')
plt.rc('lines', linewidth='2.0')
# plt.rc('font', variant='monospace')
plt.rc('legend', numpoints=1, fontsize='x-large')
plt.rc('text', usetex=True)

label_size = 20 # 12 for paper.
extra_size = 2
legend_size = 16 # 20 for paper.
plt.rc('xtick', labelsize=label_size)
plt.rc('ytick', labelsize=label_size)


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

def plot_order_parameter(ionization_fraction, order_parameter, model_tags):

	ax1 = plt.subplot(111)

	for i in xrange(0, len(ionization_fraction)):
	
		ax1.plot(np.subtract(1, ionization_fraction[i]), order_parameter[i], lw = 3, ls = '-', label = model_tags[i], color = colors[i])


	ax1.set_xlabel(r"$x_\mathrm{HII}$", size = label_size)
 	ax1.set_ylabel(r"$P_\mathrm{inf}/x_\mathrm{HII}$", size = label_size)      

	leg = ax1.legend(loc='upper right', numpoints=1,
			 labelspacing=0.1)
	leg.draw_frame(False)  # Don't want a box frame
	for t in leg.get_texts():  # Reduce the size of the text
		   t.set_fontsize(legend_size)
 
	plt.tight_layout()

	outputFile = './hosh_comparison.png' 
	plt.savefig(outputFile)  # Save the figure
	print 'Saved file to', outputFile
	plt.close()

##

if __name__ == '__main__':

	ionization_fraction = np.arange(0.025, 1, 0.025)
	#order_parameter = run_algorithm(ionization_fraction)
	#order_parameter = [8.617387973956339e-05, 8.63052713341836e-05, 0.00010859254480240948, 0.00016267747993445054, 0.0004856949453306359, 0.00576613083064141, 0.7136131125149394, 0.8894560383494995, 0.9503515610072987, 0.9763831118285013, 0.9887944314981405, 0.9947896898242197, 0.997792191391382, 0.9991435810856882, 0.9997078420539404, 0.9999194561680467, 0.9999831585660305, 0.9999978795218534, 1.0]

	order_parameter = [[0.0, 0.0002600717855728328, 0.0001524281017914023, 0.00030435464307665015, 0.0002310397864098811, 0.00018833491725983578, 0.00015799961943269454, 0.0001337031333757895, 0.0001525441528027305, 0.00013211035421659787, 0.00011539402771016814, 0.00010182217867472144, 0.00011339241864335947, 0.00010186042537244762, 9.191521888849925e-05, 9.949067640830607e-05, 8.987988946446047e-05, 8.151525039876147e-05, 7.39405818408903e-05, 6.66173733307235e-05, 7.030944764634071e-05, 7.258501149520218e-05, 0.0003114293323236657, 0.00038118038000661017, 0.0005544993179623913, 0.0005840516032551386, 0.0006038736563358097, 0.0007389455342275873, 0.0009985895288745327, 0.0011543526797691158, 0.0021008603168896135, 0.005183920494658369, 0.026295820736016263, 0.2893099261153975, 0.6790186095274281, 0.8146890100710948, 0.8814153538218651, 0.9161085538987875, 0.9402805947479239, 0.9587870118029799, 0.9695256132262272, 0.9778378077551334, 0.983452187956986, 0.9875154363910544, 0.990540786168907, 0.9930090279275877, 0.9940295213043593, 0.9963880636757952, 0.9975231431213227], [0.003449834569863321, 0.00722269504690302, 0.008811504055806172, 0.0093149943844322, 0.009367512339044098, 0.009633008846130293, 0.009354662382135514, 0.008962898766188035, 0.009086893414209759, 0.009859272081505574, 0.009912142875591463, 0.01231311824873962, 0.011985480997543672, 0.012309076652455543, 0.012284714065883233, 0.0172470246141306, 0.022272939660760025, 0.025795645199905005, 0.06629693959084412, 0.06772066380310063, 0.07116502077469639, 0.07435712191929296, 0.0840162329284189, 0.12043528510475125, 0.15175416686107712, 0.21607625325864632, 0.2626273094252207, 0.37502332097445523, 0.3983689681132447, 0.4915966526906067, 0.5530684964126101, 0.6435670064218341, 0.6787071838944451, 0.7205779395461404, 0.7478145804247244, 0.7828786709161016, 0.798599895233373, 0.8098167376367655, 0.8566573437687678, 0.8776220993663492, 0.8868625941609881, 0.8996873956930244, 0.9162180699338778, 0.9246971892029022, 0.93469329906986, 0.9464485253645982, 0.9545566513661694, 0.9651809199752057, 0.9689353709057305], [0.0, 0.0002444507143383135, 0.0002879828506113652, 0.00030115309343995566, 0.00023274018223679838, 0.00018903348823184123, 0.0001574637643011322, 0.00013302116922008362, 0.0001138620994344768, 0.00013165908901626984, 0.0001149887925068919, 0.0001268879180847866, 0.000135597825619836, 0.0001215334837799708, 0.00010972384615424451, 9.876622566698821e-05, 8.907077028526146e-05, 8.005698307930536e-05, 7.15336762753797e-05, 6.340413445861318e-05, 0.0006144888446897098, 0.001129066338702871, 0.0012872351808379502, 0.0012980270211040508, 0.0013095120313362958, 0.0017940217345819229, 0.0025295725793437703, 0.0031526876435801483, 0.0071256150101336934, 0.018745652444882265, 0.0646566932084082, 0.4467642364530629, 0.6588657545205925, 0.7711209767570868, 0.8421693847767213, 0.8855540965835584, 0.9152386021447962, 0.9349072892747065, 0.95115834873681, 0.9625329892604764, 0.9713202572566492, 0.9775890873738998, 0.9826634695718339, 0.9865665320597946, 0.9891875229531089, 0.9915759923074391, 0.9934800979413023, 0.9949923755388796, 0.9961905174600497]] 

	ionization_fraction = [[0.99917077997523962, 0.9981665171515901, 0.99687172409418523, 0.99529985328908188, 0.99380837605142369, 0.99240442773213544, 0.99094610809984418, 0.98930083806945091, 0.98749641597027271, 0.98556245917193319, 0.9834709935110062, 0.98126784696970915, 0.97897402825038649, 0.97659360068154077, 0.97406103341920081, 0.97124330588047048, 0.96816837485820417, 0.9649019915264565, 0.96130645880802956, 0.95705290067479643, 0.95252615090632708, 0.94744511040165436, 0.94181725954802809, 0.93495065977390412, 0.91142599151356751, 0.90447770847845066, 0.89261022925040057, 0.87610354196037032, 0.86199941543947667, 0.84096514945075385, 0.8238694730007281, 0.80591399793933383, 0.786549724172659, 0.76584192607723833, 0.7433345109872207, 0.72026892896736672, 0.69443190689022005, 0.66922717415374378, 0.64130030123690174, 0.61021664791819108, 0.57953270854059546, 0.54580066531955784, 0.50910973505947987, 0.47058817144581366, 0.42872879954815024, 0.3841615254452434, 0.33551981099352562, 0.28278506040017048, 0.22714831420630349], [0.99930889851593374, 0.99867961431264474, 0.99772715752978236, 0.99667263188746491, 0.99536780098856936, 0.99386195854673798, 0.99209915264685244, 0.99005137167002633, 0.98772078751963277, 0.98510378418280398, 0.98196011331046296, 0.97846837376161189, 0.97441852414705377, 0.96993879839948771, 0.96522132371509273, 0.959717588682935, 0.95388542085537698, 0.94722480960953304, 0.93989992119492072, 0.93212952871691312, 0.92360825049310158, 0.91453023623709828, 0.90475882965826249, 0.89426739598167737, 0.88301706788615597, 0.87089114365899201, 0.85830548717582666, 0.84492185529613029, 0.83078742588067567, 0.81594794546164406, 0.80038171620326171, 0.78405934265906596, 0.7666719335018084, 0.7484153690499612, 0.7292302575113645, 0.70897401109218627, 0.68797195739554895, 0.666405189084053, 0.64364680766045523, 0.61963414294643671, 0.59449550315090882, 0.56849899138831761, 0.53988515778263357, 0.51008094667625892, 0.47826743631505042, 0.44590171304902854, 0.41483988736426752, 0.38501117737040863, 0.35589569887669203], [0.99912129880617739, 0.9980493525678833, 0.99668843365366488, 0.99524988616829668, 0.9938536119510557, 0.99243249707768755, 0.99091529736407369, 0.98924598631182836, 0.98743645618942255, 0.98551297409799943, 0.98341274317931304, 0.98121030097268591, 0.97890067236594225, 0.97645897360764822, 0.97392524005040015, 0.97103237539049725, 0.96787921626751483, 0.96426266842475705, 0.96000453075828496, 0.95487639767267429, 0.93947278512659615, 0.93073808751117837, 0.92072687774968842, 0.90963058625836901, 0.89767849531069499, 0.88438035068360121, 0.86993152841437305, 0.85419709388469056, 0.83745438735769429, 0.81949753498533795, 0.80010312108855375, 0.779558295444577, 0.75720714238835118, 0.73335442360917014, 0.70787538210709344, 0.6807066205192559, 0.65199024369796577, 0.62216930029876083, 0.5899732743336632, 0.55613703937792569, 0.52083062371311795, 0.48456771445332719, 0.44615297105768242, 0.40709092677389735, 0.36723508716561898, 0.3259876305142288, 0.28449859036385133, 0.24211208862736555, 0.20058030364074028]]

	
	plot_order_parameter(ionization_fraction, order_parameter, ["Homo + Random", "Homo + Correct", "Correct + Random"])


