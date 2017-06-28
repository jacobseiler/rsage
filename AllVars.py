#!/usr/bin/env python

import numpy as np
from astropy import units as u
from astropy import cosmology

def set_cosmology(Hubble_h, Omega_m):
	
    cosmo = cosmology.FlatLambdaCDM(H0 = Hubble_h*100, Om0 = Omega_m) 
    t_BigBang = cosmo.lookback_time(100000).value # Lookback time to the Big Bang in Gyr.

    return cosmo, t_BigBang

def Set_Params_MiniMill():
    
    print "Setting parameters to Mini Millennium."
    
    global Hubble_h
    global Omega_m
    global Omega_L
    global BoxSize
    global Volume
    global SnapZ
    global BaryonFrac
    
    global SnapZ
    global Lookback_Time

    global cosmo
    global t_BigBang
    
    Hubble_h = 0.73
    Omega_m = 0.25
    Omega_L = 0.75
    BoxSize = 62.5 # Mpc/h
    Volume = BoxSize**3
    BaryonFrac = 0.17
    
    SnapZ = [ 1.27000000e+02, 7.99978940e+01, 4.99995900e+01, 3.00000630e+01,
             1.99156900e+01,  1.82437230e+01, 1.67245250e+01, 1.53430740e+01,
             1.40859140e+01,  1.29407800e+01, 1.18965700e+01, 1.09438640e+01,
             1.00734610e+01,  9.27791500e+00, 8.54991200e+00, 7.88320350e+00,
             7.27218800e+00,  6.71158650e+00, 6.19683360e+00, 5.72386400e+00,
             5.28883360e+00,  4.88844900e+00, 4.51955560e+00, 4.17946860e+00,
             3.86568280e+00,  3.57590500e+00, 3.30809780e+00, 3.06041900e+00,
             2.83118270e+00,  2.61886140e+00, 2.42204400e+00, 2.23948550e+00,
             2.07002740e+00,  1.91263270e+00, 1.76633580e+00, 1.63027070e+00,
             1.50363650e+00,  1.38571810e+00, 1.27584620e+00, 1.17341690e+00,
             1.07787450e+00,  9.88708140e-01, 9.05462400e-01, 8.27699100e-01,
             7.55035640e-01,  6.87108800e-01, 6.23590100e-01, 5.64176600e-01,
             5.08591400e-01,  4.56577240e-01, 4.07899440e-01, 3.62340270e-01,
             3.19703430e-01,  2.79801800e-01, 2.42469090e-01, 2.07548630e-01,
             1.74897610e-01,  1.44383420e-01, 1.15883370e-01, 8.92878300e-02,
             6.44933950e-02,  4.14030630e-02, 1.99325420e-02, 0.00000000e+00]

    Lookback_Time = [13.5672,    13.5551,    13.5305,    13.4761,    13.3929,    13.3680,    13.3403,
                 13.3093,    13.2748,    13.2365,    13.1940,    13.1470,    13.0951,    13.0378,
                 12.9748,    12.9055,    12.8296,    12.7465,    12.6558,    12.5569,    12.4494,
                 12.3327,    12.2064,    12.0699,    11.9227,    11.7644,    11.5945,    11.4127,
                 11.2186,    11.0119,    10.7924,    10.5598,    10.3142,    10.0557,    9.78421,
                 9.50011,    9.20376,    8.89563,    8.57635,    8.24662,    7.90726,    7.55927,
                 7.20365,    6.84147,    6.47396,    6.10229,    5.72772,    5.35152,    4.97492,
                 4.59916,    4.22544,    3.85492,    3.48873,    3.12788,    2.77339,    2.42617,
                 2.08709,    1.75692,    1.43640,    1.12622,    0.82696,    0.53917,
                 0.263375,   0.      ] # In Gyr. 

    cosmo, t_BigBang = set_cosmology(Hubble_h, Omega_m)

    print "######################"
    print "BoxSize = %.3f (Mpc/h)" %(BoxSize)
    print "Hubble_h = %.3f" %(Hubble_h)
    print "Omega_m = %.3f" %(Omega_m)
    print "Omega_L = %.3f" %(Omega_L)
    print "BaryonFrac = %.3f" %(BaryonFrac)
    print "t_BigBang = %.3f Gyr" %(t_BigBang)
    print "######################"

def Set_Params_Mysim():

    print "Setting parameters to my simulation."
    
    global Hubble_h
    global Omega_m
    global Omega_L
    global Omega_b
    global BoxSize
    global Volume
    global SnapZ
    global BaryonFrac
    global Y
    
    global SnapZ
    global Lookback_Time

    global cosmo
    global t_BigBang
    
    Hubble_h = 0.678
    Omega_m = 0.308
    Omega_L = 0.692
    Omega_b = 0.0484
    BoxSize = 100 # Mpc/h
    Volume = BoxSize**3
    BaryonFrac = 0.17
    Y = 0.24

 
    SnapZ = [49.000000,     35.000001,      32.025452,      29.597752,      27.571988,      25.851283,      24.368212,
             23.074276,     21.933618,      20.919111,      20.009818,      19.189299,      18.444455,      17.764701,
             17.141388,     16.567374,      16.036700,      15.544359,      15.086102,      14.658304,      14.257849,
             13.882045,     13.528548,      13.195313,      12.880543,      12.582654,      12.300241,      12.032056,
             11.776987,     11.534035,      11.302305,      11.080989,      10.869359,      10.666753,      10.472570,
             10.286264,     10.107335,      9.935327,       9.769820,       9.610431,       9.456804,       9.308614,
             9.165559,      9.027360,       8.893758,       8.764514,       8.639403,       8.518218,       8.400767,
             8.286867,      8.176350,       8.069058,       7.964843,       7.863565,       7.765095,       7.669310,
             7.576095,      7.485340,       7.396944,       7.310809,       7.226845,       7.144965,       7.065088,
             6.987135,      6.911035,       6.836717,       6.764116,       6.693168,       6.623815,       6.555999,
             6.489667,      6.424768,       6.361251,       6.299072,       6.238185,       6.178547,       6.120119,
             6.062861,      6.006736,       5.951709,       5.897746,       5.844813,       5.792881,       5.741918,
             5.691896,      5.642788,       5.594568,       5.547208,       5.500687,       5.454979,       5.410062,
             5.365914,      5.322515,       5.279845,       5.237883,       5.196612,       5.156013,       5.116069,
             5.076763,      5.038078,       5.000000,       5.000000]

    Lookback_Time = [13.7561, 13.7249, 13.7138, 13.7028, 13.6917, 13.6806, 13.6695,
                 13.6584, 13.6474, 13.6363, 13.6252, 13.6141, 13.6030, 13.5920,
                 13.5809, 13.5698, 13.5587, 13.5477, 13.5366, 13.5255, 13.5144,
                 13.5033, 13.4923, 13.4812, 13.4701, 13.4590, 13.4479, 13.4369,
                 13.4258, 13.4147, 13.4036, 13.3926, 13.3815, 13.3704, 13.3593,
                 13.3482, 13.3372, 13.3261, 13.3150, 13.3039, 13.2928, 13.2818,
                 13.2707, 13.2596, 13.2485, 13.2375, 13.2264, 13.2153, 13.2042,
                 13.1931, 13.1821, 13.1710, 13.1599, 13.1488, 13.1377, 13.1267,
                 13.1156, 13.1045, 13.0934, 13.0824, 13.0713, 13.0602, 13.0491,
                 13.0380, 13.0270, 13.0159, 13.0048, 12.9937, 12.9826, 12.9716,
                 12.9605, 12.9494, 12.9383, 12.9273, 12.9162, 12.9051, 12.8940,
                 12.8829, 12.8719, 12.8608, 12.8497, 12.8386, 12.8275, 12.8165,
                 12.8054, 12.7943, 12.7832, 12.7722, 12.7611, 12.7500, 12.7389,
                 12.7278, 12.7168, 12.7057, 12.6946, 12.6835, 12.6724, 12.6614,
                 12.6503, 12.6392, 12.6281, 12.6281] # In Gyr.


    cosmo, t_BigBang = set_cosmology(Hubble_h, Omega_m)

    print "######################"
    print "BoxSize = %.3f (Mpc/h)" %(BoxSize)
    print "Hubble_h = %.3f" %(Hubble_h)
    print "Omega_m = %.3f" %(Omega_m)
    print "Omega_L = %.3f" %(Omega_L)
    print "BaryonFrac = %.3f" %(BaryonFrac)
    print "t_BigBang = %.3f Gyr" %(t_BigBang)
    print "######################"

    return cosmo

def Set_Constants():

    print "Setting constants (in cgs units)."

    global Proton_Mass
    global Solar_Mass
    global Sec_Per_Year
    global Sec_Per_Megayear
    global Ionized_Mass_H
    global LymanAlpha_Energy
    global eV_to_erg
    global M_Bol_Sun
    global L_Sun
    global W_to_ergs
    global A_to_m
    global c_in_ms
    global pc_to_m
    global m_to_cm
    global Sigmat
    global H0

    Proton_Mass = 1.6726219e-24 # Grams.
    Solar_Mass = 1.98855e33 # Grams.
    Sec_Per_Year = 3.155e7
    Sec_Per_Megayear = 3.155e13
    Ionized_Mass_H = 0.53 # Molecular mass of ionized H.
    LymanAlpha_Energy = 10.2 # Energy of Lyman Alpha photon in eV.
    eV_to_erg = 1.6202e-12 # Conversion from eV to erg. 
    M_Bol_Sun = 4.74 # Bolometric magnitude of the sun.
    L_Sun = 3.828e26 # Luminosity of the Sun in W.
    W_to_ergs = 1.0e7 # Watts to erg s^-1.
    A_to_m = 1.0e-10 # Angstroms to meters.
    c_in_ms = 3.0e8 # Speed of light in m s^-1 
    pc_to_m = 3.086e16 # Parsec to metres. 
    m_to_cm = 1.0e2 # Meters to centimetres.
    Sigmat = 6.652e-29 # Thomson cross-section for an electron in m^2.
    H0 = 3.24078e-18 # In h/sec

def spectralflux_wavelength_to_frequency(Flux, Wavelength):

    # For a given spectral flux at a specific wavelength (f_lamba), this function converts it to the spectral flux density at a specific frequency (f_nu).

    ### INPUT ###
    # Flux: The spectral flux density in units of erg s^-1 A^-1 cm^-2
    # Wavelength: The wavelength we are determining the spectral flux density at in units of Angstroms (A).

    ### OUTPUT ###
    # f_nu: The spectral flux density at a specific frequency in units of Janksy (W Hz^-1 m^-2).

    f_nu = 3.34e4 * pow(Wavelength, 2) * Flux
  
    return f_nu

def Luminosity_to_Flux(Luminosity, Distance):

    # Converts a luminosity to an observed flux at some distance.
    ## NOTE THAT INPUT AND OUTPUT ARE IN LOG UNITS.

    ### INPUT ###
    # Luminosity: The intrinisic luminosity being converted in units of erg s^-1 A^-1. <<<NOTE MUST BE IN LOG UNITS>>>.
    # Distance: The distance at which the flux is being observed in units of parsec.

    ### OUTPUT ###
    # Flux: The observed flux in units of erg s^-1 A^-1 cm^-2. 

    F = Luminosity - np.log10(4*np.pi*pow(Distance * pc_to_m * m_to_cm, 2.0))

    return F

def Luminosity_to_ABMag(Luminosity, Wavelength):
     
    # Converts an intrinsic luminosity into absolute AB magnitude at a specified wavelength.
    ## NOTE THAT INPUT IS IN LOG.

    ### INPUT ###
    # Luminosity: The intrinsic luminosity of the star in units of erg s^-1 A^-1.  <<NOTE MUST BE IN LOG UNITS>>>.
    # Wavelength: The wavelength we want the magnitude at.

    ### OUTPUT ###
    # M: The absolute magnitude in the AB system.

    Flux = Luminosity_to_Flux(Luminosity, 10.0) # Calculate the flux from a distance of 10 parsec, units of erg s^-1 A^-1 cm^-2.  Log Units. 
    f_nu = spectralflux_wavelength_to_frequency(10**Flux, 1600) # Spectral flux density in Janksy.
    M = -2.5 * np.log10(f_nu) + 8.90 # AB Magnitude from http://www.astro.ljmu.ac.uk/~ikb/convert-units/node2.html

    print "Flux from AllVars.py", Flux
    print "M from AllVars.py", M
    return M


def Set_Params_Tiamat():

    print "Setting parameters to Tiamat" 
    
    global Hubble_h
    global Omega_m
    global Omega_L
    global Omega_b
    global BoxSize
    global Volume
    global SnapZ
    global BaryonFrac
    global Y

    global SnapZ
    global Lookback_Time

    global cosmo
    global t_BigBang
    
    Hubble_h = 0.678
    Omega_m = 0.308
    Omega_L = 0.692
    Omega_b = 0.0484
    BoxSize = 100*Hubble_h # Mpc/h
    Volume = BoxSize**3
    BaryonFrac = 0.17
    Y = 0.24   
 
    SnapZ = [49.000000,     35.000001,      32.025452,      29.597752,      27.571988,      25.851283,      24.368212,
             23.074276,     21.933618,      20.919111,      20.009818,      19.189299,      18.444455,      17.764701,
             17.141388,     16.567374,      16.036700,      15.544359,      15.086102,      14.658304,      14.257849,
             13.882045,     13.528548,      13.195313,      12.880543,      12.582654,      12.300241,      12.032056,
             11.776987,     11.534035,      11.302305,      11.080989,      10.869359,      10.666753,      10.472570,
             10.286264,     10.107335,      9.935327,       9.769820,       9.610431,       9.456804,       9.308614,
             9.165559,      9.027360,       8.893758,       8.764514,       8.639403,       8.518218,       8.400767,
             8.286867,      8.176350,       8.069058,       7.964843,       7.863565,       7.765095,       7.669310,
             7.576095,      7.485340,       7.396944,       7.310809,       7.226845,       7.144965,       7.065088,
             6.987135,      6.911035,       6.836717,       6.764116,       6.693168,       6.623815,       6.555999,
             6.489667,      6.424768,       6.361251,       6.299072,       6.238185,       6.178547,       6.120119,
             6.062861,      6.006736,       5.951709,       5.897746,       5.844813,       5.792881,       5.741918,
             5.691896,      5.642788,       5.594568,       5.547208,       5.500687,       5.454979,       5.410062,
             5.365914,      5.322515,       5.279845,       5.237883,       5.196612,       5.156013,       5.116069,
             5.076763,      5.038078,       5.000000,       5.000000]

    Lookback_Time = [13.7561, 13.7249, 13.7138, 13.7028, 13.6917, 13.6806, 13.6695,
                 13.6584, 13.6474, 13.6363, 13.6252, 13.6141, 13.6030, 13.5920,
                 13.5809, 13.5698, 13.5587, 13.5477, 13.5366, 13.5255, 13.5144,
                 13.5033, 13.4923, 13.4812, 13.4701, 13.4590, 13.4479, 13.4369,
                 13.4258, 13.4147, 13.4036, 13.3926, 13.3815, 13.3704, 13.3593,
                 13.3482, 13.3372, 13.3261, 13.3150, 13.3039, 13.2928, 13.2818,
                 13.2707, 13.2596, 13.2485, 13.2375, 13.2264, 13.2153, 13.2042,
                 13.1931, 13.1821, 13.1710, 13.1599, 13.1488, 13.1377, 13.1267,
                 13.1156, 13.1045, 13.0934, 13.0824, 13.0713, 13.0602, 13.0491,
                 13.0380, 13.0270, 13.0159, 13.0048, 12.9937, 12.9826, 12.9716,
                 12.9605, 12.9494, 12.9383, 12.9273, 12.9162, 12.9051, 12.8940,
                 12.8829, 12.8719, 12.8608, 12.8497, 12.8386, 12.8275, 12.8165,
                 12.8054, 12.7943, 12.7832, 12.7722, 12.7611, 12.7500, 12.7389,
                 12.7278, 12.7168, 12.7057, 12.6946, 12.6835, 12.6724, 12.6614,
                 12.6503, 12.6392, 12.6281, 12.6281] # In Gyr.


    cosmo, t_BigBang = set_cosmology(Hubble_h, Omega_m)

    print "######################"
    print "BoxSize = %.3f (Mpc/h)" %(BoxSize)
    print "Hubble_h = %.3f" %(Hubble_h)
    print "Omega_m = %.3f" %(Omega_m)
    print "Omega_L = %.3f" %(Omega_L)
    print "BaryonFrac = %.3f" %(BaryonFrac)
    print "t_BigBang = %.3f Gyr" %(t_BigBang)
    print "######################"

    return cosmo


def Set_Params_Tiamat_extended():

    print "Setting parameters to the extended Tiamat simulation (that ran down to z = 1.8ish)." 
    
    global Hubble_h
    global Omega_m
    global Omega_L
    global Omega_b
    global BoxSize
    global Volume
    global SnapZ
    global BaryonFrac
    global Y

    global SnapZ
    global Lookback_Time

    global cosmo
    global t_BigBang
    
    Hubble_h = 0.678
    Omega_m = 0.308
    Omega_L = 0.692
    Omega_b = 0.0484
    BoxSize = 100*Hubble_h # Mpc/h
    Volume = BoxSize**3
    BaryonFrac = 0.17
    Y = 0.24   
 
    SnapZ = [49.000000,     35.000001,      32.025452,      29.597752,      27.571988,      25.851283,      24.368212,
             23.074276,     21.933618,      20.919111,      20.009818,      19.189299,      18.444455,      17.764701,
             17.141388,     16.567374,      16.036700,      15.544359,      15.086102,      14.658304,      14.257849,
             13.882045,     13.528548,      13.195313,      12.880543,      12.582654,      12.300241,      12.032056,
             11.776987,     11.534035,      11.302305,      11.080989,      10.869359,      10.666753,      10.472570,
             10.286264,     10.107335,      9.935327,       9.769820,       9.610431,       9.456804,       9.308614,
             9.165559,      9.027360,       8.893758,       8.764514,       8.639403,       8.518218,       8.400767,
             8.286867,      8.176350,       8.069058,       7.964843,       7.863565,       7.765095,       7.669310,
             7.576095,      7.485340,       7.396944,       7.310809,       7.226845,       7.144965,       7.065088,
             6.987135,      6.911035,       6.836717,       6.764116,       6.693168,       6.623815,       6.555999,
             6.489667,      6.424768,       6.361251,       6.299072,       6.238185,       6.178547,       6.120119,
             6.062861,      6.006736,       5.951709,       5.897746,       5.844813,       5.792881,       5.741918,
             5.691896,      5.642788,       5.594568,       5.547208,       5.500687,       5.454979,       5.410062,
             5.365914,      5.322515,       5.279845,       5.237883,       5.196612,       5.156013,       5.116069,
             5.076763,      5.038078,       5.000000,       5.000000,       4.928760,       4.858360,	    4.788800,
             4.720060,      4.652140,       4.585030,       4.518710,       4.453180,       4.388430,       4.324450,
             4.261230,      4.198750,       4.137020,       4.076030,       4.015750,       3.956190,       3.897350,
             3.839190,      3.781730,       3.724960,       3.668850,       3.613410,       3.558630,       3.504500,
	     3.451020,      3.398170,       3.345940,       3.294340,       3.243350,       3.192960,       3.143170,
	     3.093980,      3.045370,       2.997330,       2.949870,       2.902970,       2.856620,       2.810830,
	     2.765580,      2.720870,       2.676690,       2.633030,       2.589890,       2.547260,       2.505140,
	     2.463520,      2.422400,       2.381760,       2.341610,       2.301930,       2.262720,       2.223980,
	     2.185700,      2.147870,       2.110490,       2.073560,       2.037060,       2.001000,       1.965370,
	     1.930160,      1.895360,       1.860980,       1.82701]

    cosmo, t_BigBang = set_cosmology(Hubble_h, Omega_m)

    Lookback_Time = cosmo.lookback_time(SnapZ).value # In Gyr 
                                 
    print "######################"
    print "BoxSize = %.3f (Mpc/h)" %(BoxSize)
    print "Hubble_h = %.3f" %(Hubble_h)
    print "Omega_m = %.3f" %(Omega_m)
    print "Omega_L = %.3f" %(Omega_L)
    print "BaryonFrac = %.3f" %(BaryonFrac)
    print "t_BigBang = %.3f Gyr" %(t_BigBang)
    print "######################"

    return cosmo

