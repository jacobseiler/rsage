#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



void init_galaxy(int p, int halonr)
{
  int j, step;

	assert(halonr == Halo[halonr].FirstHaloInFOFgroup);

  Gal[p].Type = 0;

  Gal[p].GalaxyNr = GalaxyCounter;
  GalaxyCounter++;
  
  Gal[p].HaloNr = halonr;
  Gal[p].MostBoundID = Halo[halonr].MostBoundID;
  Gal[p].SnapNum = Halo[halonr].SnapNum - 1;

  Gal[p].mergeType = 0;
  Gal[p].mergeIntoID = -1;
  Gal[p].mergeIntoSnapNum = -1;
  Gal[p].dT = -1.0;

  for(j = 0; j < 3; j++)
  {
    Gal[p].Pos[j] = Halo[halonr].Pos[j];
    Gal[p].Vel[j] = Halo[halonr].Vel[j];
  }

  Gal[p].Len = Halo[halonr].Len;
  Gal[p].Vmax = Halo[halonr].Vmax;
  Gal[p].Vvir = get_virial_velocity(halonr);
  Gal[p].Mvir = get_virial_mass(halonr);
  Gal[p].Rvir = get_virial_radius(halonr);

  Gal[p].deltaMvir = 0.0;

  Gal[p].ColdGas = 0.0;
  Gal[p].StellarMass = 0.0;
  Gal[p].BulgeMass = 0.0;
  Gal[p].HotGas = 0.0;
  Gal[p].EjectedMass = 0.0;
  Gal[p].BlackHoleMass = 0.0;
  Gal[p].ICS = 0.0;

  Gal[p].MetalsColdGas = 0.0;
  Gal[p].MetalsStellarMass = 0.0;
  Gal[p].MetalsBulgeMass = 0.0;
  Gal[p].MetalsHotGas = 0.0;
  Gal[p].MetalsEjectedMass = 0.0;
  Gal[p].MetalsICS = 0.0;
  
  for(step = 0; step < STEPS; step++)
  {
    Gal[p].SfrDisk[step] = 0.0;
    Gal[p].SfrBulge[step] = 0.0;
    Gal[p].SfrDiskColdGas[step] = 0.0;
    Gal[p].SfrDiskColdGasMetals[step] = 0.0;
    Gal[p].SfrBulgeColdGas[step] = 0.0;
    Gal[p].SfrBulgeColdGasMetals[step] = 0.0;
  }

  Gal[p].DiskScaleRadius = get_disk_radius(halonr, p);
  Gal[p].MergTime = 999.9;
  Gal[p].Cooling = 0.0;
  Gal[p].Heating = 0.0;
  Gal[p].r_heat = 0.0;
  Gal[p].QuasarModeBHaccretionMass = 0.0;
  Gal[p].TimeOfLastMajorMerger = -1.0;
  Gal[p].TimeOfLastMinorMerger = -1.0;
  Gal[p].OutflowRate = 0.0;
	Gal[p].TotalSatelliteBaryons = 0.0;
	// infall properties
  Gal[p].infallMvir = -1.0;  
  Gal[p].infallVvir = -1.0;
  Gal[p].infallVmax = -1.0;
 
  if(NULL == (Gal[p].GridHistory = malloc(sizeof(int) * MAXSNAPS)))
  {
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridHistory.", sizeof(int)*MAXSNAPS); 
    exit(EXIT_FAILURE);
  }

  if(NULL == (Gal[p].GridStellarMass = malloc(sizeof(float) * MAXSNAPS)))
  {
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridStellarMass.", sizeof(float)*MAXSNAPS); 
    exit(EXIT_FAILURE);
  } 

  if(NULL == (Gal[p].GridSFR = malloc(sizeof(float) * MAXSNAPS)))
  {
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridSFR.", sizeof(float)*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  if (NULL == (Gal[p].GridZ = malloc(sizeof(float) * MAXSNAPS)))
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridSFR.", sizeof(float)*MAXSNAPS);
    exit(EXIT_FAILURE);
  }
 
  if (NULL == (Gal[p].GridCentralGalaxyMass = malloc(sizeof(float) * MAXSNAPS)))
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridCentralGalaxyMass.", sizeof(float)*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  if (NULL == (Gal[p].GridPhotons_HI = malloc(sizeof(float) * MAXSNAPS))) 
  {
    fprintf(stderr, "Out of memoery allocating %ld bytes, could not allocate GridPhotons_HI.", sizeof(float)*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  if (NULL == (Gal[p].GridPhotons_HeI = malloc(sizeof(float) * MAXSNAPS)))
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridPhotons_HeI.", sizeof(float)*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  if (NULL == (Gal[p].GridPhotons_HeII = malloc(sizeof(float) * MAXSNAPS)))
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridPhotons_HeII.", sizeof(float)*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  if (NULL == (Gal[p].MfiltGnedin = malloc(sizeof(double) * MAXSNAPS)))
  {
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate MfiltGnedin.", sizeof(float)*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  if (NULL == (Gal[p].MfiltSobacchi = malloc(sizeof(double) * MAXSNAPS)))
  {   
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate MfiltSobacchi.", sizeof(float)*MAXSNAPS);
    exit(EXIT_FAILURE);
  }


  for (j = 0; j < MAXSNAPS; ++j)
  {
    Gal[p].GridHistory[j] = -1;
    Gal[p].GridStellarMass[j] = 0.0;
    Gal[p].GridSFR[j] = 0.0;
    Gal[p].GridZ[j] = -1;
    Gal[p].GridCentralGalaxyMass[j] = -1.0;
    Gal[p].GridPhotons_HI[j] = 0.0;
    Gal[p].GridPhotons_HeI[j] = 0.0;
    Gal[p].GridPhotons_HeII[j] = 0.0;
    Gal[p].MfiltGnedin[j] = 1.0;
    Gal[p].MfiltSobacchi[j] = 1.0;
  }
 
}



double get_disk_radius(int halonr, int p)
{
  double SpinMagnitude, SpinParameter;
  
	if(Gal[p].Vvir > 0.0 && Gal[p].Rvir > 0.0)
	{
		// See Mo, Shude & White (1998) eq12, and using a Bullock style lambda.
		SpinMagnitude = sqrt(Halo[halonr].Spin[0] * Halo[halonr].Spin[0] + 
			Halo[halonr].Spin[1] * Halo[halonr].Spin[1] + Halo[halonr].Spin[2] * Halo[halonr].Spin[2]);
  
		SpinParameter = SpinMagnitude / (1.414 * Gal[p].Vvir * Gal[p].Rvir);
		return (SpinParameter / 1.414) * Gal[p].Rvir;		
	}
	else
		return 0.1 * Gal[p].Rvir;

}



double get_metallicity(double gas, double metals)
{
  double metallicity;

  if(gas > 0.0 && metals > 0.0)
  {
    metallicity = metals / gas;
    if(metallicity < 1.0)
      return metallicity;
    else
      return 1.0;
  }
  else
    return 0.0;

}



double dmax(double x, double y)
{
  if(x > y)
    return x;
  else
    return y;
}



double get_virial_mass(int halonr)
{
  if(halonr == Halo[halonr].FirstHaloInFOFgroup && Halo[halonr].Mvir >= 0.0)
  {

    return Halo[halonr].Mvir;   /* take spherical overdensity mass estimate */
  } 
  else
  {

    return Halo[halonr].Len * PartMass;
  }
}



double get_virial_velocity(int halonr)
{
	double Rvir;
	
	Rvir = get_virial_radius(halonr);
	
  if(Rvir > 0.0)
		return sqrt(G * get_virial_mass(halonr) / Rvir);
	else
		return 0.0;
}



double get_virial_radius(int halonr)
{
  // return Halo[halonr].Rvir;  // Used for Bolshoi

  double zplus1, hubble_of_z_sq, rhocrit, fac;
  
  zplus1 = 1 + ZZ[Halo[halonr].SnapNum];
  hubble_of_z_sq =
    Hubble * Hubble *(Omega * zplus1 * zplus1 * zplus1 + (1 - Omega - OmegaLambda) * zplus1 * zplus1 +
    OmegaLambda);
  
  rhocrit = 3 * hubble_of_z_sq / (8 * M_PI * G);
  fac = 1 / (200 * 4 * M_PI / 3.0 * rhocrit);
  
  return cbrt(get_virial_mass(halonr) * fac);
}

// INPUT: 
// Galaxy index (p).
// Halo index (halonr).
// Number of steps completed (for merged galaxies this is <= STEPS otherwise = STEPS) (steps_completed)
//
// OUTPUT/USE:
// Tracks a number of properties that will be used by the gridding code.  
// These properties are in the form of a length SNAPNUM array so only the values at the Galaxy redshift will be altered. 

void update_grid_array(int p, int halonr, int steps_completed, int centralgal)
{
    int x_grid, y_grid, z_grid, grid_position, step;
    int SnapCurr = Halo[halonr].SnapNum;
    x_grid = Gal[p].Pos[0]*GridSize/BoxSize; // Convert the (x,y,z) position to a grid (x,y,z).
    y_grid = Gal[p].Pos[1]*GridSize/BoxSize;
    z_grid = Gal[p].Pos[2]*GridSize/BoxSize; 

    grid_position = (x_grid*GridSize+y_grid)*GridSize+z_grid; // Convert the grid (x,y,z) to a 1D value.

    if(grid_position > CUBE(GridSize) || grid_position < 0) // Sanity check to ensure that no Grid Positions are outside the box.
    {
	printf("Found a Grid Position outside the bounds of the box or negative; grid_position = %d, Galaxy Index = %d, halonr = %d\n", grid_position, p, halonr);
	exit(0);
    }

    Gal[p].GridPos = grid_position; 

    // NOTE: We use the Snapshot number of the FOF-Halo (i.e. the main halo the galaxy belongs to) because the snapshot number of the galaxy has been shifted by -1. //
    // This is self-consistent with the end of the 'evolve_galaxies' function which shifts Gal[p].SnapNum by +1. // 
  
    Gal[p].GridHistory[SnapCurr] = grid_position; // Remember the grid history of the galaxy over the Snapshot range.
    Gal[p].GridStellarMass[SnapCurr] = Gal[p].StellarMass; // Stellar mass at this snapshot.

    for(step = 0; step < steps_completed; step++)
    {
      Gal[p].GridSFR[SnapCurr] += Gal[p].SfrBulge[step] + Gal[p].SfrDisk[step]; // Star formation rate at this snapshot.
    }

    Gal[p].GridZ[SnapCurr] = get_metallicity(Gal[p].ColdGas, Gal[p].MetalsColdGas); // Metallicity at this snapshot.
    Gal[p].GridCentralGalaxyMass[SnapCurr] = get_virial_mass(Halo[Gal[p].HaloNr].FirstHaloInFOFgroup); // Virial mass of the central galaxy (i.e. virial mass of the host halo).  
    //printf("Gal[p].Halonr = %d \t %.4e \t get_virial = %.4e \t Halo[halonr].Mvir = %.4e\n", halonr, Gal[p].GridCentralGalaxyMass[SnapCurr], get_virial_mass(halonr), Halo[halonr].Mvir);
//    printf("Gal[p].Halonr = %d \t %.4e \t get_virial = %.4e \t Halo[halonr].Mvir = %.4e\n", halonr, Gal[p].GridCentralGalaxyMass[SnapCurr], get_virial_mass(halonr), Halo[halonr].Mvir);
    float SFR_conversion = UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS / STEPS;
    float Ngamma_HI, Ngamma_HeI, Ngamma_HeII; 
    calculate_photons(Gal[p].GridSFR[SnapCurr]*SFR_conversion, Gal[p].GridZ[SnapCurr], &Ngamma_HI, &Ngamma_HeI, &Ngamma_HeII);
    Gal[p].GridPhotons_HI[SnapCurr] = Ngamma_HI; 
    Gal[p].GridPhotons_HeI[SnapCurr] = Ngamma_HeI; 
    Gal[p].GridPhotons_HeII[SnapCurr] = Ngamma_HeII; 
    Gal[p].MfiltGnedin[SnapCurr] = do_reionization(centralgal, ZZ[SnapCurr], 1);
    if (ReionizationOn == 2)
    {  
      Gal[p].MfiltSobacchi[SnapCurr] = do_myreionization(centralgal, ZZ[SnapCurr], 1); 
    }
}

// INPUT:
// Star formation rate of the galaxy (in units of Msun/yr) (SFR).
// Metallicity (NOT Solar Units) (Z).
//
// OUTPUT/USE:
// Returns the number of HI ionizing photons for the galaxy.  
//
// NOTE: These relationships have been fit manually from STARBURST99 results.  
// DOUBLE NOTE: These relationships assume a constant starformation scenario; a Starburst scenario is completely different.
void calculate_photons(float SFR, float Z, float *Ngamma_HI, float *Ngamma_HeI, float *Ngamma_HeII)
{

  assert(Ngamma_HI);
  assert(Ngamma_HeI);
  assert(Ngamma_HeII);

  if (SFR == 0)
  {
    *Ngamma_HI = 0;
    *Ngamma_HeI = 0;
    *Ngamma_HeII = 0; 
  }
 
  
  if (Z < 0.0025) // 11
  { 
    *Ngamma_HI = log10(SFR) + 53.354;
    *Ngamma_HeI = log10(SFR) + 52.727;
    *Ngamma_HeII = log10(SFR) + 48.941;
  }
  else if (Z >= 0.0025 && Z < 0.006) // 12
  {
    *Ngamma_HI = log10(SFR) + 53.290;
    *Ngamma_HeI = log10(SFR) + 52.583;
    *Ngamma_HeII = log10(SFR) + 49.411;
  }
  else if (Z>= 0.006 && Z < 0.014) // 13
  {
    *Ngamma_HI = log10(SFR) + 53.248;
    *Ngamma_HeI = log10(SFR) + 52.481;
    *Ngamma_HeII = log10(SFR) + 49.254;
  }
  else if (Z >= 0.014 && Z < 0.030) // 14
  {
    *Ngamma_HI = log10(SFR) + 53.166;
    *Ngamma_HeI = log10(SFR) + 52.319;
    *Ngamma_HeII = log10(SFR) + 48.596;
  }
  else // 15
  {
    *Ngamma_HI = log10(SFR) + 53.041;
    *Ngamma_HeI = log10(SFR) + 52.052;
    *Ngamma_HeII = log10(SFR) + 47.939;
  }

  if (SFR != 0)
  {
    assert(*Ngamma_HI > 0.0);
    assert(*Ngamma_HeI > 0.0);
    assert(*Ngamma_HeII > 0.0);
  }
}
