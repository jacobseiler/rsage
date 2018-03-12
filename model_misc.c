#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



void init_galaxy(int p, int halonr, int treenr)
{
  int j, step;
  
  ++count_gal;
	assert(halonr == Halo[halonr].FirstHaloInFOFgroup);

  Gal[p].Type = 0;
  Gal[p].TreeNr = treenr;

  Gal[p].GalaxyNr = GalaxyCounter;
  GalaxyCounter++;
  
  Gal[p].HaloNr = halonr;
  Gal[p].MostBoundID = Halo[halonr].MostBoundID;
  //Gal[p].MostBoundID = -1; 
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
 
  Gal[p].IsMerged = -1;

//  if(Gal[p].IsMalloced != 1)
//  {
  	malloc_grid_arrays(&Gal[p]);
  ++gal_mallocs;	
//  }
  
  for (j = 0; j < MAXSNAPS; ++j)
  {
    Gal[p].GridHistory[j] = -1;
    Gal[p].GridStellarMass[j] = 0.0;
    Gal[p].GridSFR[j] = 0.0;
    Gal[p].GridZ[j] = -1;
    Gal[p].GridCentralGalaxyMass[j] = -1.0;
    Gal[p].MfiltGnedin[j] = 1.0;
    Gal[p].MfiltSobacchi[j] = 1.0;
    Gal[p].EjectedFraction[j] = -1.0;
    Gal[p].LenHistory[j] = -1;
    Gal[p].GridOutflowRate[j] = 0.0;
    Gal[p].GridInfallRate[j] = 0.0;
    Gal[p].GridEjectedMass[j] = 0.0;
    Gal[p].QuasarActivity[j] = 0;
    Gal[p].DynamicalTime[j] = 0.0;
    Gal[p].QuasarSubstep[j] = -1;
    Gal[p].GridColdGas[j] = 0.0;
    Gal[p].LenMergerGal[j] = -1;
    Gal[p].GridBHMass[j] = 0.0;
  }

  for (j = 0; j < SN_Array_Len; ++j)
  {
    Gal[p].Stars[j] = 0.0;
  }

  Gal[p].Total_SF_Time = 0.0;
  Gal[p].Total_Stars = 0.0;

  Gal[p].GrandSum = 0.0;
 
  Gal[p].StellarAge_Numerator = 0.0;
  Gal[p].StellarAge_Denominator = 0.0;

  Gal[p].reheated_mass = 0.0;
  Gal[p].ejected_mass = 0.0;
  Gal[p].mass_stars_recycled = 0.0;
  Gal[p].mass_metals_new = 0.0; 


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
    ++count_Mvir;
    return Halo[halonr].Mvir;   /* take spherical overdensity mass estimate */
  } 
  else
  {    
    ++count_Len;
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

int32_t determine_1D_idx(float pos_x, float pos_y, float pos_z, int32_t *grid_1D)
{

  int32_t x_grid, y_grid, z_grid;

  x_grid = pos_x * GridSize/BoxSize;
  y_grid = pos_y * GridSize/BoxSize;
  z_grid = pos_z * GridSize/BoxSize;

  *grid_1D = (z_grid*GridSize+y_grid)*GridSize+x_grid; // Convert the grid (x,y,z) to a 1D value.

  if(*grid_1D > CUBE(GridSize) || *grid_1D < 0) // Sanity check to ensure that no Grid Positions are outside the box.
  {
    fprintf(stderr, "Found a Grid Position outside the bounds of the box or negative; grid_position = %d\nPos[0] = %.4f\t Pos[1] = %.4f\tPos[2] = %.4f", *grid_1D, pos_x, pos_y, pos_z); 
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
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
    int32_t grid_position, step, status;
    int32_t SnapCurr = Halo[halonr].SnapNum;
    double MfiltSobacchi, reionization_modifier; 

    status = determine_1D_idx(Gal[p].Pos[0], Gal[p].Pos[1], Gal[p].Pos[2], &grid_position); 
    if (status == EXIT_FAILURE)
    {
      exit(EXIT_FAILURE);
    }
    
    Gal[p].GridPos = grid_position; 

    // NOTE: We use the Snapshot number of the FOF-Halo (i.e. the main halo the galaxy belongs to) because the snapshot number of the galaxy has been shifted by -1. //
    // This is self-consistent with the end of the 'evolve_galaxies' function which shifts Gal[p].SnapNum by +1. //
    Gal[p].GridHistory[SnapCurr] = grid_position; // Remember the grid history of the galaxy over the Snapshot range.
    Gal[p].GridStellarMass[SnapCurr] = Gal[p].StellarMass; // Stellar mass at this snapshot.

    for(step = 0; step < steps_completed; step++) // We loop over the number of steps completed to allow merged galaxies to be updated. 
    {
      Gal[p].GridSFR[SnapCurr] += Gal[p].SfrBulge[step] + Gal[p].SfrDisk[step]; // Star formation rate at this snapshot.
    }

    Gal[p].GridZ[SnapCurr] = get_metallicity(Gal[p].ColdGas, Gal[p].MetalsColdGas); // Metallicity at this snapshot.
    Gal[p].GridCentralGalaxyMass[SnapCurr] = get_virial_mass(Halo[Gal[p].HaloNr].FirstHaloInFOFgroup); // Virial mass of the central galaxy (i.e. virial mass of the host halo).  
    //printf("Gal[p].Halonr = %d \t %.4e \t get_virial = %.4e \t Halo[halonr].Mvir = %.4e\n", halonr, Gal[p].GridCentralGalaxyMass[SnapCurr], get_virial_mass(halonr), Halo[halonr].Mvir);
//    printf("Gal[p].Halonr = %d \t %.4e \t get_virial = %.4e \t Halo[halonr].Mvir = %.4e\n", halonr, Gal[p].GridCentralGalaxyMass[SnapCurr], get_virial_mass(halonr), Halo[halonr].Mvir);
//    float SFR_conversion = UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS / STEPS;
//    float Ngamma_HI, Ngamma_HeI, Ngamma_HeII; 
//    calculate_photons(Gal[p].GridSFR[SnapCurr]*SFR_conversion, Gal[p].GridZ[SnapCurr], &Ngamma_HI, &Ngamma_HeI, &Ngamma_HeII);
//    Gal[p].GridPhotons_HI[SnapCurr] = Ngamma_HI; 
//    Gal[p].GridPhotons_HeI[SnapCurr] = Ngamma_HeI; 
//    Gal[p].GridPhotons_HeII[SnapCurr] = Ngamma_HeII; 
//    Gal[p].MfiltGnedin[SnapCurr] = do_reionization(centralgal, ZZ[SnapCurr], 1);
    if (ReionizationOn == 2) 
    {  
      reionization_modifier = do_myreionization(centralgal, ZZ[SnapCurr], &MfiltSobacchi);
      Gal[p].MfiltSobacchi[SnapCurr] = MfiltSobacchi; 
    }
    else
    {
      reionization_modifier = 0.0;
      Gal[p].MfiltSobacchi[SnapCurr] = 1.0;
    }

 
    if((Gal[p].EjectedMass < 0.0) || ((Gal[p].HotGas + Gal[p].ColdGas + Gal[p].EjectedMass) == 0.0))
    {
      Gal[p].EjectedFraction[SnapCurr] = 0.0; // Check divide by 0 case.
    }
    else
    { 
        Gal[p].EjectedFraction[SnapCurr] = Gal[p].EjectedMass/(Gal[p].HotGas + Gal[p].ColdGas + Gal[p].EjectedMass);
        // EjectedFraction is the fraction of baryons in the ejected reservoir.
    }
    if (Gal[p].EjectedFraction[SnapCurr] < 0.0 || Gal[p].EjectedFraction[SnapCurr] > 1.0)
    {
      fprintf(stderr, "Found ejected fraction = %.4e \t p = %d \t Gal[p].EjectedMass = %.4e \t Gal[p].HotGas = %.4e \t Gal[p].ColdGas = %.4e\n\n", Gal[p].EjectedFraction[SnapCurr], p, Gal[p].EjectedMass, Gal[p].HotGas, Gal[p].ColdGas); 
      return(EXIT_FAILURE);
    }

    Gal[p].LenHistory[SnapCurr] = Gal[p].Len;
    if (Gal[p].LenHistory[SnapCurr] < 0)
    {  
      fprintf(stderr, "Have a galaxy with Len < 0.  Galaxy number %d with Len %d.\n", p, Gal[p].Len);
      exit(EXIT_FAILURE);
    }
    if (Gal[p].EjectedMass > 0.0)
    {
        Gal[p].GridEjectedMass[SnapCurr] += Gal[p].EjectedMass; 
    }

    Gal[p].DynamicalTime[SnapCurr] = Gal[p].Rvir / Gal[p].Vvir; 
    Gal[p].GridColdGas[SnapCurr] = Gal[p].ColdGas;
    Gal[p].GridBHMass[SnapCurr] = Gal[p].BlackHoleMass;
}

