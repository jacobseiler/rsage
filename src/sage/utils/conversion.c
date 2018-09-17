#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>

#include "conversion.h"

#define pc_to_m 3.086e16 // Parsec to metres.
#define m_to_cm 1.0e2 // Meters to centimetres.

/*
   Converts a given spectral flux at a specific wavelength (f_lambda) to the
  spectral flux density at the equivalent frequency (f_nu).

  Parameters
  ----------

  flux : Float.
      The spectral flux density at the wavelength `wavelength`.

  wavelength : Float.
      The wavelength we're converting at.

  Returns
  ---------

  f_nu : Float.
      The spectral flux density at the `wavelength` equivalent frequency.

  Units
  ---------

  `flux` units are erg s^-1 A^-1 cm^-2 (where A is Angstroms).
  `wavelength` units are A.
  `f_nu` units are Janksy (W Hz^-1 m^-2).
*/

float spectralflux_wavelength_to_frequency(float flux, float wavelength)
{

    float f_nu = 3.34e4 * pow(wavelength, 2) * flux;

    return f_nu; 
}

/*
  Converts an intrinsic luminosity to an observed flux.

  Parameters
  ----------

  luminosity : Float.
      The instrinsic luminosity we're converting.

  distance : Float.
      The distance the flux is being measured/observed at.

  Returns
  ---------

  flux : Float.
      The observed flux.

  Units
  ---------

  `luminosity` units are log10(erg s^-1 A^-1) (where A is Angstroms).
  `distance` units are pc.
  `flux` units are log10(erg s^-1 A^-1 cm^-2).
*/

float luminosity_to_flux(float luminosity, float distance)
{

    float flux = luminosity - log10(4*3.14159*pow(distance * pc_to_m * m_to_cm, 2.0));

    return flux;
}

/*
  Converts an intrinsic luminosity to an absolute AB magnitude at a specified
  wavelength.

  Parameters
  ----------

  luminosity : Float.
      The instrinsic luminosity we're converting.

  wavelength : Float.
      The wavelength we're calculating the magnitude at.

  Returns
  ---------

  M : Float.
      The AB absolute magnitude.

  Units
  ---------

  `luminosity` units are log10(erg s^-1 A^-1) (where A is Angstroms).
  `wavelength` units are A.
  `M` units are unitless.
*/

float luminosity_to_ABMag(float luminosity, float wavelength)
{
    float flux = luminosity_to_flux(luminosity, 10.0); // Calculate the flux from a distance of 10 parsec, units of erg s^-1 A^-1 cm^-2.  Log Units.
    float f_nu = spectralflux_wavelength_to_frequency(pow(10, flux), wavelength); // Spectral flux density in Janksy.
    float M = -2.5 * log10(f_nu) + 8.90; //AB Magnitude from http://www.astro.ljmu.ac.uk/~ikb/convert-units/node2.html

    return M;
}
