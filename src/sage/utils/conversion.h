#ifndef CONVERSION_H 
#define CONVERSION_H 

// Proto-Types //
float spectralflux_wavelength_to_frequency(float flux, float wavelength);
float luminosity_to_flux(float luminosity, float distance);
float luminosity_to_ABMag(float luminosity, float wavelength);

#endif
