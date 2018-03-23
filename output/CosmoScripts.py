#!/usr/bin/env python

import numpy as np
import AllVars

# Hubble parameter as a function of redshift.
## Input ##
# Cosmology: array containing [Hubble_h, Omega_m, Omega_l]
# z: Redshift calculating the Hubble parameter at.

## Output ## 
# Hubble parameter at given redshift in Units of (km/s/Mpc)
def Hz(cosmology, z):

  return cosmology['Hubble_h'] * AllVars.H0 * np.sqrt(cosmology['Omega_m'] * pow(1.0 + z, 3) + cosmology['Omega_l'])


def integralGrowth(cosmology, z):

	dz = 0.001
	integral = 0.0
	
	Hz0 = Hz(cosmology, 0.0)

	for z_dummy in np.arange(z+dz/2.0, z + 100.0, dz):
		Hubz = Hz(cosmology, z_dummy)/Hz0
		integral = (1.0+z_dummy)/pow(Hubz, 3)
#		print "integral = ", integral
	integral *= dz	

	return integral


def getGrowth(cosmology, z):

	growth = Hz(cosmology, z)/Hz(cosmology, 0.0) * integralGrowth(cosmology, z) / integralGrowth(cosmology, 0.0) 

	return growth 

def dzdt(cosmology, z):

	return (-Hz(cosmology, z) * (1.0 + z))

def dgrowthdt(cosmology, z):
	
	dz = 1.0e-3 
	tmp1 = getGrowth(cosmology, z+dz/2.0)
	tmp2 = getGrowth(cosmology, z-dz/2.0)
	tmp3 = dzdt(cosmology, z)
	print "growth(z+dz/2.0) = ", tmp1
	print "growth(z-dz/2.0) = ", tmp2
	print "dzdt(z) = ", tmp3 
	print "growth(z+dz/2.0) - growth(z-dz/2.0) = ", tmp1 - tmp2	
	return ((getGrowth(cosmology, z + dz/2.0) - getGrowth(cosmology, z - dz/2.0)) / dz * dzdt(cosmology,z))
