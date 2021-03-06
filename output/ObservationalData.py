#!/usr/bin/env python

from __future__ import print_function
import numpy as np
from astropy import units as u
from astropy import cosmology
from scipy import stats
import os

def Get_Data_SMF():

    ## Stellar Mass Function Data.##

    print("Getting SMF data.")

    global Gonzalez_SMF_z6
    global Gonzalez_SMF_z7

    global Song_SMF_z6
    global Song_SMF_z7
    global Song_SMF_z8

    global Duncan_SMF_z6
    global Duncan_SMF_z7

    global Baldry_SMF_z0

    global SMF_colors
    global SMF_markers

    #SMF_colors = ["b", "g", "c", "m"]
    SMF_colors = ["#31a354", "#2b8cbe", "#c51b8a", "m"]
    SMF_markers = ["o", "s", "D"]

    ## For Gonzalez and Song, Columns are redshift, mean value, lower bound, upper bound.##
    ## Units of log10 Phi [dex^-1 Mpc^-3].
    ## Gonzales Data uses h = 0.699999 and a Salpeter IMF ## 
    ## Song Data uses h = 0.7 and a Salpeter IMF ## 
 
    Gonzalez_SMF_z6 = np.array([[7.77, -2.0956, -1.8596, -2.3539],
                [8.27, -2.1742, -1.9494, -2.4101],
                [8.77, -2.5674, -2.3876, -2.7921],
                [9.27, -2.8483, -2.6573, -3.0843],
                [9.77, -3.5787, -3.3764, -3.8258],
                [10.27, -4.3202, -4.0281, -4.5674]], dtype = np.float32)

    Gonzalez_SMF_z7 = np.array([[7.75, -2.1828, -1.7463, -2.5858],
                [8.26, -2.25, -1.8694, -2.2631],
                [8.77, -2.7425, -2.3731, -3.1231],
                [9.27, -3.0672, -2.6753, -3.4142],
                [9.76, -3.8731, -3.4831, -4.2537]], dtype = np.float32)

    Song_SMF_z6 = np.array([[7.25, -1.47, -1.47 + 0.35, -1.47 - 0.23],
                [7.75, -1.81, -1.81 + 0.23, -1.81 - 0.28],
                [8.25, -2.26, -2.26 + 0.21, -2.26 - 0.16],
                [8.75, -2.65, -2.65 + 0.15, -2.65 - 0.15],
                [9.25, -3.14, -3.14 + 0.12, -3.14 - 0.11],
                [9.75, -3.69, -3.69 + 0.12, -3.69 - 0.13],
                [10.25, -4.27, -4.27 + 0.38, -4.27 - 0.86]], dtype = np.float32)

    Song_SMF_z7 = np.array([[7.25, -1.63, -1.63 + 0.54, -1.63 - 0.54],
                [7.75, -2.07, -2.07 + 0.45, -2.07 - 0.41],
                [8.25, -2.49, -2.49 + 0.38, -2.49 - 0.32],
                [8.75, -2.96, -2.96 + 0.32, -2.96 - 0.30],
                [9.25, -3.47, -3.47 + 0.32, -3.47 - 0.35],
                [9.75, -4.11, -4.11 + 0.41, -4.11 - 0.57],
                [10.25, -4.61, -4.61 + 0.72, -4.61 - 0.82],
                [10.75, -5.24, -5.24 + 0.90, -5.25 - 0.57]], dtype = np.float32)

    Song_SMF_z8 = np.array([[7.25, -1.73, -1.73 + 1.01, -1.73 - 0.84],
                [7.75, -2.28, -2.28 + 0.84, -2.28 - 0.64],
                [8.25, -2.88, -2.88 + 0.75, -2.88 - 0.57],
                [8.75, -3.45, -3.45 + 0.57, -3.45 - 0.60],
                [9.25, -4.21, -4.21 + 0.63, -4.21 - 0.78],
                [9.75, -5.31, -5.31 + 1.01, -5.31 - 1.64]], dtype = np.float32)


    ## For Baldry, Columns are Stellar Mass, lower bound, upper bound.##
    ## Units of Phi [dex^-1 Mpc^-3].

    Baldry_SMF_z0 = np.array([
        [7.05, 1.3531e-01, 6.0741e-02],
        [7.15, 1.3474e-01, 6.0109e-02],
        [7.25, 2.0971e-01, 7.7965e-02],
        [7.35, 1.7161e-01, 3.1841e-02],
        [7.45, 2.1648e-01, 5.7832e-02],
        [7.55, 2.1645e-01, 3.9988e-02],
        [7.65, 2.0837e-01, 4.8713e-02],
        [7.75, 2.0402e-01, 7.0061e-02],
        [7.85, 1.5536e-01, 3.9182e-02],
        [7.95, 1.5232e-01, 2.6824e-02],
        [8.05, 1.5067e-01, 4.8824e-02],
        [8.15, 1.3032e-01, 2.1892e-02],
        [8.25, 1.2545e-01, 3.5526e-02],
        [8.35, 9.8472e-02, 2.7181e-02],
        [8.45, 8.7194e-02, 2.8345e-02],
        [8.55, 7.0758e-02, 2.0808e-02],
        [8.65, 5.8190e-02, 1.3359e-02],
        [8.75, 5.6057e-02, 1.3512e-02],
        [8.85, 5.1380e-02, 1.2815e-02],
        [8.95, 4.4206e-02, 9.6866e-03],
        [9.05, 4.1149e-02, 1.0169e-02],
        [9.15, 3.4959e-02, 6.7898e-03],
        [9.25, 3.3111e-02, 8.3704e-03],
        [9.35, 3.0138e-02, 4.7741e-03],
        [9.45, 2.6692e-02, 5.5029e-03],
        [9.55, 2.4656e-02, 4.4359e-03],
        [9.65, 2.2885e-02, 3.7915e-03],
        [9.75, 2.1849e-02, 3.9812e-03],
        [9.85, 2.0383e-02, 3.2930e-03],
        [9.95, 1.9929e-02, 2.9370e-03],
        [10.05, 1.8865e-02, 2.4624e-03],
        [10.15, 1.8136e-02, 2.5208e-03],
        [10.25, 1.7657e-02, 2.4217e-03],
        [10.35, 1.6616e-02, 2.2784e-03],
        [10.45, 1.6114e-02, 2.1783e-03],
        [10.55, 1.4366e-02, 1.8819e-03],
        [10.65, 1.2588e-02, 1.8249e-03],
        [10.75, 1.1372e-02, 1.4436e-03],
        [10.85, 9.1213e-03, 1.5816e-03],
        [10.95, 6.1125e-03, 9.6735e-04],
        [11.05, 4.3923e-03, 9.6254e-04],
        [11.15, 2.5463e-03, 5.0038e-04],
        [11.25, 1.4298e-03, 4.2816e-04],
        [11.35, 6.4867e-04, 1.6439e-04],
        [11.45, 2.8294e-04, 9.9799e-05],
        [11.55, 1.0617e-04, 4.9085e-05],
        [11.65, 3.2702e-05, 2.4546e-05],
        [11.75, 1.2571e-05, 1.2571e-05],
        [11.85, 8.4589e-06, 8.4589e-06],
        [11.95, 7.4764e-06, 7.4764e-06],
        ], dtype=np.float32)
    
    ## For Duncan, Columns are log Stellar Mass, Phi, Error in Phi Lower, Error in Phi Upper ## 
    ## Units of Phi [dex^-1 Mpc^-3].
    ## Duncan data uses h = 0.7, Chabrier IMF. ##

    Duncan_SMF_z6 = np.loadtxt("/home/jseiler/ObservationalData/Duncan14_MF_z6.txt")    
    Duncan_SMF_z7 = np.loadtxt("/home/jseiler/ObservationalData/Duncan14_MF_z7.txt")    

def Get_Data_SMBH():

    ## This is data for the StellarMass - BlackHole relationship.
    global Mstar
    global Huang_z8_BHSM
    global SMBH_colors

    SMBH_colors = ["b", "g", "c", "m"]

    Mstar = np.arange(1.0, 12.0, 0.01) # Array of stellar masses.

    Huang_z8_BHSM = 1.10*Mstar - 3.85 # Relationship is log10(M_BH) = 1.10*log10(M_*/10^11) + 8.25.  Units in Msun.


def get_data_UVLF():

    global Bouwens2015_UVLF_z6
    global Qin2017_UVLF_z6

    global Bouwens2015_UVLF_z7
    global Atek2015_UVLF_z7

    global Bouwens2015_UVLF_z8

    Bouwens2015_UVLF_z6 = np.genfromtxt("/home/jseiler/ObservationalData/GLF_UV/Bouwens2015/z6pt0.dat")
    Qin2017_UVLF_z6 = np.genfromtxt("/home/jseiler/ObservationalData/GLF_UV/Qin2017_Tiamat/z6pt0.dat")

    Bouwens2015_UVLF_z7 = np.genfromtxt("/home/jseiler/ObservationalData/GLF_UV/Bouwens2015/z7pt0.dat")
    Atek2015_UVLF_z7 = np.genfromtxt("/home/jseiler/ObservationalData/GLF_UV/Atek2015/z7pt0.dat")

    Bouwens2015_UVLF_z8 = np.genfromtxt("/home/jseiler/ObservationalData/GLF_UV/Bouwens2015/z8pt0.dat")

    global UVLF_colors
    global UVLF_markers
    global UVLF_linestyles
    global UVLF_linewidths

    UVLF_colors = ["#31a354", "#2b8cbe", "#c51b8a", "m", "r", "b", "g"]
    UVLF_markers = ["o", "s", "D", "x", "O", "p", "8"]
    UVLF_linestyles = ["-", "--", "-."]
    UVLF_linewidth = 3
