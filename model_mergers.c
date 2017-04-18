#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdbool.h>

#include "core_allvars.h"
#include "core_proto.h"



double estimate_merging_time(int sat_halo, int mother_halo, int ngal)
{
  double coulomb, mergtime, SatelliteMass, SatelliteRadius;

  if(sat_halo == mother_halo) 
  {
    printf("\t\tSnapNum, Type, IDs, sat radius:\t%i\t%i\t%i\t%i\t--- sat/cent have the same ID\n", 
      Gal[ngal].SnapNum, Gal[ngal].Type, sat_halo, mother_halo);
    return -1.0;
  }
  
  coulomb = log(Halo[mother_halo].Len / ((double) Halo[sat_halo].Len) + 1);

  SatelliteMass = get_virial_mass(sat_halo) + Gal[ngal].StellarMass + Gal[ngal].ColdGas;
  SatelliteRadius = get_virial_radius(mother_halo);

  if(SatelliteMass > 0.0 && coulomb > 0.0)
    mergtime = 2.0 *
    1.17 * SatelliteRadius * SatelliteRadius * get_virial_velocity(mother_halo) / (coulomb * G * SatelliteMass);
  else
    mergtime = -1.0;
  
  return mergtime;

}



void deal_with_galaxy_merger(int p, int merger_centralgal, int centralgal, double time, double dt, int halonr, int step)
{
  double mi, ma, mass_ratio;

  // calculate mass ratio of merging galaxies 
  if(Gal[p].StellarMass + Gal[p].ColdGas <
    Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas)
  {
    mi = Gal[p].StellarMass + Gal[p].ColdGas;
    ma = Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas;
  }
  else
  {
    mi = Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas;
    ma = Gal[p].StellarMass + Gal[p].ColdGas;
  }

  if(ma > 0)
    mass_ratio = mi / ma;
  else
    mass_ratio = 1.0;
  add_galaxies_together(merger_centralgal, p);

  // grow black hole through accretion from cold disk during mergers, a la Kauffmann & Haehnelt (2000) 
  if(AGNrecipeOn)
  {
    grow_black_hole(merger_centralgal, mass_ratio);
  }
  // starburst recipe similar to Somerville et al. 2001

  collisional_starburst_recipe(mass_ratio, merger_centralgal, centralgal, time, dt, halonr, 0, step);

  if(mass_ratio > 0.1)
		Gal[merger_centralgal].TimeOfLastMinorMerger = time;

  if(mass_ratio > ThreshMajorMerger)
  {
    make_bulge_from_burst(merger_centralgal);
    Gal[merger_centralgal].TimeOfLastMajorMerger = time;
    Gal[p].mergeType = 2;  // mark as major merger
  }
  else
  {
    Gal[p].mergeType = 1;  // mark as minor merger
  }

}



void grow_black_hole(int merger_centralgal, double mass_ratio)
{
  double BHaccrete, metallicity;

  if(Gal[merger_centralgal].ColdGas > 0.0)
  {
    BHaccrete = BlackHoleGrowthRate * mass_ratio / 
      (1.0 + pow(280.0 / Gal[merger_centralgal].Vvir, 2.0)) * Gal[merger_centralgal].ColdGas;

    // cannot accrete more gas than is available! 
    if(BHaccrete > Gal[merger_centralgal].ColdGas)
      BHaccrete = Gal[merger_centralgal].ColdGas;

    metallicity = get_metallicity(Gal[merger_centralgal].ColdGas, Gal[merger_centralgal].MetalsColdGas);
    Gal[merger_centralgal].BlackHoleMass += BHaccrete;
    Gal[merger_centralgal].ColdGas -= BHaccrete;
    Gal[merger_centralgal].MetalsColdGas -= metallicity * BHaccrete;

    Gal[merger_centralgal].QuasarModeBHaccretionMass += BHaccrete;

    quasar_mode_wind(merger_centralgal, BHaccrete);
  }
}



void quasar_mode_wind(int gal, float BHaccrete)
{
  float quasar_energy, cold_gas_energy, hot_gas_energy;

  // work out total energies in quasar wind (eta*m*c^2), cold and hot gas (1/2*m*Vvir^2)
  quasar_energy = QuasarModeEfficiency * 0.1 * BHaccrete * (C / UnitVelocity_in_cm_per_s) * (C / UnitVelocity_in_cm_per_s);
  cold_gas_energy = 0.5 * Gal[gal].ColdGas * Gal[gal].Vvir * Gal[gal].Vvir;
  hot_gas_energy = 0.5 * Gal[gal].HotGas * Gal[gal].Vvir * Gal[gal].Vvir;
   
  // compare quasar wind and cold gas energies and eject cold
  if(quasar_energy > cold_gas_energy)
  {
    Gal[gal].EjectedMass += Gal[gal].ColdGas;
    Gal[gal].MetalsEjectedMass += Gal[gal].MetalsColdGas;

    Gal[gal].ColdGas = 0.0;
    Gal[gal].MetalsColdGas = 0.0;
  }
  
  // compare quasar wind and cold+hot gas energies and eject hot
  if(quasar_energy > cold_gas_energy + hot_gas_energy)
  {
    Gal[gal].EjectedMass += Gal[gal].HotGas;
    Gal[gal].MetalsEjectedMass += Gal[gal].MetalsHotGas; 

    Gal[gal].HotGas = 0.0;
    Gal[gal].MetalsHotGas = 0.0;
  }
}



void add_galaxies_together(int t, int p)
{
   // t is the index of the central galaxy for the merger.
   // p is the index of the merging galaxy.

  int step, i;

  Gal[t].ColdGas += Gal[p].ColdGas;
  Gal[t].MetalsColdGas += Gal[p].MetalsColdGas;
  
  Gal[t].StellarMass += Gal[p].StellarMass;
  Gal[t].MetalsStellarMass += Gal[p].MetalsStellarMass;

  Gal[t].HotGas += Gal[p].HotGas;
  Gal[t].MetalsHotGas += Gal[p].MetalsHotGas;
  
  Gal[t].EjectedMass += Gal[p].EjectedMass;
  Gal[t].MetalsEjectedMass += Gal[p].MetalsEjectedMass;
  
  Gal[t].ICS += Gal[p].ICS;
  Gal[t].MetalsICS += Gal[p].MetalsICS;

  Gal[t].BlackHoleMass += Gal[p].BlackHoleMass;

  // add merger to bulge
	Gal[t].BulgeMass += Gal[p].StellarMass;
	Gal[t].MetalsBulgeMass += Gal[p].MetalsStellarMass;		

  for(step = 0; step < STEPS; step++)
  {
    Gal[t].SfrBulge[step] += Gal[p].SfrDisk[step] + Gal[p].SfrBulge[step];
    Gal[t].SfrBulgeColdGas[step] += Gal[p].SfrDiskColdGas[step] + Gal[p].SfrBulgeColdGas[step];
    Gal[t].SfrBulgeColdGasMetals[step] += Gal[p].SfrDiskColdGasMetals[step] + Gal[p].SfrBulgeColdGasMetals[step];
  }



  // Our delayed SN scheme requires the stars formed by current galaxy and also its progenitors; so need to go back through the central galaxy of the merger and add all the stars from the merging galaxy.
  //fprintf(stderr, "Adding merger stars\n");
  for(i = 0; i < SN_Array_Len; ++i) // Careful that we add up to the current snapshot number (inclusive) as we need to account for the stars just formed.
  {
  //  fprintf(stderr, "i in add together mergers = %d\n", i); 
    Gal[t].Stars[i] += Gal[p].Stars[i];

  }
  //fprintf(stderr, "Finished adding merger stars\n");  

  Gal[t].PreviousReheatedMass[Gal[t].SnapNum] += Gal[p].PreviousReheatedMass[Gal[p].SnapNum];
  
}



void make_bulge_from_burst(int p)
{
  int step;
  
  // generate bulge 
  Gal[p].BulgeMass = Gal[p].StellarMass;
  Gal[p].MetalsBulgeMass = Gal[p].MetalsStellarMass;

  // update the star formation rate 
  for(step = 0; step < STEPS; step++)
  {
    Gal[p].SfrBulge[step] += Gal[p].SfrDisk[step];
    Gal[p].SfrBulgeColdGas[step] += Gal[p].SfrDiskColdGas[step];
    Gal[p].SfrBulgeColdGasMetals[step] += Gal[p].SfrDiskColdGasMetals[step];
    Gal[p].SfrDisk[step] = 0.0;
    Gal[p].SfrDiskColdGas[step] = 0.0;
    Gal[p].SfrDiskColdGasMetals[step] = 0.0;
  }
}



void collisional_starburst_recipe(double mass_ratio, int merger_centralgal, int centralgal, double time, double dt, int halonr, int mode, int step)
{
  double stars, eburst;
  double FracZleaveDiskVal; 
  double reheated_mass = 0.0; 
  double mass_metals_new = 0.0; 
  double mass_stars_recycled = 0.0; 
  double ejected_mass = 0.0;

  // This is the major and minor merger starburst recipe of Somerville et al. 2001. 
  // The coefficients in eburst are taken from TJ Cox's PhD thesis and should be more accurate then previous. 

  // the bursting fraction 
  if(mode == 1)
    eburst = mass_ratio;
  else
    eburst = 0.56 * pow(mass_ratio, 0.7);


  stars = eburst * Gal[merger_centralgal].ColdGas;
  if(stars < 0.0)
    stars = 0.0;
 
  if(SupernovaRecipeOn == 1)
  {    
    if(IRA == 1)
      do_current_SN(merger_centralgal, merger_centralgal, halonr, &stars, &reheated_mass, &mass_metals_new, &mass_stars_recycled, &ejected_mass);   
  } 
  
  if(stars > Gal[merger_centralgal].ColdGas) // we do this check in 'do_current_sn()' but if supernovarecipeon == 0 then we don't do the check.
  {
    double factor = Gal[merger_centralgal].ColdGas / stars; 
    stars *= factor; 
  }

  update_from_star_formation(merger_centralgal, stars, dt, step, true); 
  update_from_SN_feedback(merger_centralgal, merger_centralgal, reheated_mass, ejected_mass, mass_stars_recycled, mass_metals_new);

  // check for disk instability
  if(DiskInstabilityOn && mode == 0)
    if(mass_ratio < ThreshMajorMerger)
    check_disk_instability(merger_centralgal, centralgal, halonr, time, dt, step);

}



void disrupt_satellite_to_ICS(int centralgal, int gal)
{  
  Gal[centralgal].HotGas += Gal[gal].ColdGas + Gal[gal].HotGas;
  Gal[centralgal].MetalsHotGas += Gal[gal].MetalsColdGas + Gal[gal].MetalsHotGas;
  
  Gal[centralgal].EjectedMass += Gal[gal].EjectedMass;
  Gal[centralgal].MetalsEjectedMass += Gal[gal].MetalsEjectedMass;
  
  Gal[centralgal].ICS += Gal[gal].ICS;
  Gal[centralgal].MetalsICS += Gal[gal].MetalsICS;

  Gal[centralgal].ICS += Gal[gal].StellarMass;
  Gal[centralgal].MetalsICS += Gal[gal].MetalsStellarMass;
  
  // what should we do with the disrupted satellite BH?
  
  Gal[gal].mergeType = 4;  // mark as disruption to the ICS
  
}

void add_galaxy_to_merger_list(int p)
{
  int j, step;

  MergedGal[MergedNr].Type = Gal[p].Type;

  MergedGal[MergedNr].GalaxyNr = Gal[p].GalaxyNr; 
  
  MergedGal[MergedNr].HaloNr = Gal[p].HaloNr;
  MergedGal[MergedNr].MostBoundID = Gal[p].MostBoundID; 
  MergedGal[MergedNr].SnapNum = Gal[p].SnapNum; 

  MergedGal[MergedNr].mergeType = Gal[p].mergeType;
  MergedGal[MergedNr].mergeIntoID = Gal[p].mergeIntoID; 
  MergedGal[MergedNr].mergeIntoSnapNum = Gal[p].mergeIntoSnapNum;
  MergedGal[MergedNr].dT = Gal[p].dT; 

  for(j = 0; j < 3; j++)
  {
    MergedGal[MergedNr].Pos[j] = Gal[p].Pos[j];
    MergedGal[MergedNr].Vel[j] = Gal[p].Vel[j];
  }

  MergedGal[MergedNr].Len = Gal[p].Len;
  MergedGal[MergedNr].Vmax = Gal[p].Vmax;
  MergedGal[MergedNr].Vvir = Gal[p].Vvir; 
  MergedGal[MergedNr].Mvir = Gal[p].Mvir;
  MergedGal[MergedNr].Rvir = Gal[p].Rvir;

  MergedGal[MergedNr].deltaMvir = Gal[p].deltaMvir;

  MergedGal[MergedNr].ColdGas = Gal[p].ColdGas;
  MergedGal[MergedNr].StellarMass = Gal[p].StellarMass;
  MergedGal[MergedNr].BulgeMass = Gal[p].BulgeMass;
  MergedGal[MergedNr].HotGas = Gal[p].HotGas;
  MergedGal[MergedNr].EjectedMass = Gal[p].EjectedMass;
  MergedGal[MergedNr].BlackHoleMass = Gal[p].BlackHoleMass;
  MergedGal[MergedNr].ICS = Gal[p].ICS;

  MergedGal[MergedNr].MetalsColdGas = Gal[p].MetalsColdGas;
  MergedGal[MergedNr].MetalsStellarMass = Gal[p].MetalsStellarMass;
  MergedGal[MergedNr].MetalsBulgeMass = Gal[p].MetalsBulgeMass;
  MergedGal[MergedNr].MetalsHotGas = Gal[p].MetalsHotGas;
  MergedGal[MergedNr].MetalsEjectedMass = Gal[p].MetalsEjectedMass;
  MergedGal[MergedNr].MetalsICS = Gal[p].MetalsICS;
  
  for(step = 0; step < STEPS; step++)
  {
    MergedGal[MergedNr].SfrDisk[step] = Gal[p].SfrDisk[step];
    MergedGal[MergedNr].SfrBulge[step] = Gal[p].SfrBulge[step];
    MergedGal[MergedNr].SfrDiskColdGas[step] = Gal[p].SfrDiskColdGas[step];
    MergedGal[MergedNr].SfrDiskColdGasMetals[step] = Gal[p].SfrDiskColdGasMetals[step];
    MergedGal[MergedNr].SfrBulgeColdGas[step] = Gal[p].SfrBulgeColdGas[step];

  }

  MergedGal[MergedNr].DiskScaleRadius = Gal[p].DiskScaleRadius; 
  MergedGal[MergedNr].MergTime = Gal[p].MergTime;
  MergedGal[MergedNr].Cooling = Gal[p].Cooling;
  MergedGal[MergedNr].Heating = Gal[p].Heating;
  MergedGal[MergedNr].r_heat = Gal[p].r_heat;
  MergedGal[MergedNr].QuasarModeBHaccretionMass = Gal[p].QuasarModeBHaccretionMass;
  MergedGal[MergedNr].TimeOfLastMajorMerger = Gal[p].TimeOfLastMajorMerger;
  MergedGal[MergedNr].TimeOfLastMinorMerger = Gal[p].TimeOfLastMinorMerger;
  MergedGal[MergedNr].OutflowRate = Gal[p].OutflowRate;
  MergedGal[MergedNr].TotalSatelliteBaryons = Gal[p].TotalSatelliteBaryons;
	// infall properties
  MergedGal[MergedNr].infallMvir = Gal[p].infallMvir;  
  MergedGal[MergedNr].infallVvir = Gal[p].infallVvir;
  MergedGal[MergedNr].infallVmax = Gal[p].infallVmax;

  if (NULL == (MergedGal[MergedNr].GridHistory = malloc(sizeof(*(MergedGal[MergedNr].GridHistory)) * MAXSNAPS))) 
  {   
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridHistory in model_mergers.c.", sizeof(*(MergedGal[MergedNr].GridHistory))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  if (NULL == (MergedGal[MergedNr].GridStellarMass = malloc(sizeof(*(MergedGal[MergedNr].GridStellarMass)) * MAXSNAPS))) 
  {   
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridStellarMass in model_mergers.c.", sizeof(*(MergedGal[MergedNr].GridStellarMass))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  if (NULL == (MergedGal[MergedNr].GridSFR = malloc(sizeof(*(MergedGal[MergedNr].GridSFR)) * MAXSNAPS))) 
  {   
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridSFR in model_mergers.c.", sizeof(*(MergedGal[MergedNr].GridSFR))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  if (NULL == (MergedGal[MergedNr].GridZ = malloc(sizeof(*(MergedGal[MergedNr].GridZ)) * MAXSNAPS))) 
  {   
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridZ in model_mergers.c.", sizeof(*(MergedGal[MergedNr].GridZ))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  if (NULL == (MergedGal[MergedNr].GridCentralGalaxyMass = malloc(sizeof(*(MergedGal[MergedNr].GridCentralGalaxyMass)) * MAXSNAPS))) 
  {   
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridCentralGalaxyMass in model_mergers.c.", sizeof(*(MergedGal[MergedNr].GridCentralGalaxyMass))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  /*
  if (NULL == (MergedGal[MergedNr].GridPhotons_HI = malloc(sizeof(*(MergedGal[MergedNr].GridPhotons_HI)) * MAXSNAPS))) 
  {   
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridPhotons_HI in model_Mergers.c.", sizeof(*(MergedGal[MergedNr].GridPhotons_HI))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  if (NULL == (MergedGal[MergedNr].GridPhotons_HeI = malloc(sizeof(*(MergedGal[MergedNr].GridPhotons_HeI)) * MAXSNAPS))) 
  {   
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridPhotons_HeI in model_mergers.c.", sizeof(*(MergedGal[MergedNr].GridPhotons_HeI))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  if (NULL == (MergedGal[MergedNr].GridPhotons_HeII = malloc(sizeof(*(MergedGal[MergedNr].GridPhotons_HeII)) * MAXSNAPS))) 
  {   
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridPhotons_HeII in model_mergers.c.", sizeof(*(MergedGal[MergedNr].GridPhotons_HeII))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }
  */

  if (NULL == (MergedGal[MergedNr].MfiltGnedin = malloc(sizeof(*(MergedGal[MergedNr].MfiltGnedin)) * MAXSNAPS)))
  {   
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate MfiltGnedin in model_mergers.c.", sizeof(*(MergedGal[MergedNr].MfiltGnedin))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  if (NULL == (MergedGal[MergedNr].MfiltSobacchi = malloc(sizeof(*(MergedGal[MergedNr].MfiltSobacchi)) * MAXSNAPS)))
  {   
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate MfiltSobacchi in model_mergers.c.", sizeof(*(MergedGal[MergedNr].MfiltSobacchi))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }
 
  if (NULL == (MergedGal[MergedNr].EjectedFraction = malloc(sizeof(*(MergedGal[MergedNr].EjectedFraction)) * MAXSNAPS)))
  {   
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate EjectedFraction in model_mergers.c.", sizeof(*(MergedGal[MergedNr].EjectedFraction))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }
  
  if (NULL == (MergedGal[MergedNr].LenHistory = malloc(sizeof(*(MergedGal[MergedNr].LenHistory)) * MAXSNAPS)))
  {   
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate LenHistory in model_mergers.c.", sizeof(*(MergedGal[MergedNr].LenHistory))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  if (NULL == (MergedGal[MergedNr].Stars = malloc(sizeof(*(MergedGal[MergedNr].Stars)) * SN_Array_Len)))
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate Stars in model_mergers.c.", sizeof(*(MergedGal[MergedNr].Stars))*SN_Array_Len);
    exit(EXIT_FAILURE);
  }

  if (NULL == (MergedGal[MergedNr].PreviousReheatedMass = malloc(sizeof(*(MergedGal[MergedNr].PreviousReheatedMass)) * MAXSNAPS)))
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate PreviousReheatedMass in model_mergers.c.", sizeof(*(MergedGal[MergedNr].PreviousReheatedMass))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  
  if (NULL == (MergedGal[MergedNr].VmaxHistory = malloc(sizeof(*(MergedGal[MergedNr].VmaxHistory)) * MAXSNAPS)))
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate VmaxHistory in model_mergers.c.", sizeof(*(MergedGal[MergedNr].VmaxHistory))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }
  

  for (j = 0; j < MAXSNAPS; ++j)
  {
    MergedGal[MergedNr].GridHistory[j] = Gal[p].GridHistory[j];

//    MergedGal[MergedNr].deltaEjectedMass[j] = Gal[p].deltaEjectedMass[j];
//    MergedGal[MergedNr].deltaMetalsEjectedMass[j] = Gal[p].deltaMetalsEjectedMass[j];
//    MergedGal[MergedNr].deltaEnergyEjected[j] = Gal[p].deltaEnergyEjected[j];
//    MergedGal[MergedNr].deltaSfr[j] = Gal[p].deltaSfr[j];
    MergedGal[MergedNr].GridStellarMass[j] = Gal[p].GridStellarMass[j];
    MergedGal[MergedNr].GridSFR[j] = Gal[p].GridSFR[j];
    MergedGal[MergedNr].GridZ[j] = Gal[p].GridZ[j];
    MergedGal[MergedNr].GridCentralGalaxyMass[j] = Gal[p].GridCentralGalaxyMass[j];
    MergedGal[MergedNr].MfiltGnedin[j] = Gal[p].MfiltGnedin[j];
    MergedGal[MergedNr].MfiltSobacchi[j] = Gal[p].MfiltSobacchi[j];
    MergedGal[MergedNr].EjectedFraction[j] = Gal[p].EjectedFraction[j]; 
    MergedGal[MergedNr].LenHistory[j] = Gal[p].LenHistory[j]; 
    MergedGal[MergedNr].PreviousReheatedMass[j] = Gal[p].PreviousReheatedMass[j];
    MergedGal[MergedNr].VmaxHistory[j] = Gal[p].VmaxHistory[j];
  }
 
  //fprintf(stderr, "Copying over Merger stars\n");
  for (j = 0; j < SN_Array_Len; ++j)
  {
//    fprintf(stderr, "j in add to merged list = %d\n", j);
    MergedGal[MergedNr].Stars[j] = Gal[p].Stars[j];
  }

  //fprintf(stderr, "Finished copying over Merger stars\n");
  free(Gal[p].GridHistory);
  free(Gal[p].GridStellarMass);
  free(Gal[p].GridSFR);
  free(Gal[p].GridZ);
  free(Gal[p].GridCentralGalaxyMass);
  free(Gal[p].MfiltGnedin);
  free(Gal[p].MfiltSobacchi);
  free(Gal[p].EjectedFraction);
  free(Gal[p].LenHistory);
  free(Gal[p].Stars);
  free(Gal[p].PreviousReheatedMass);
  free(Gal[p].VmaxHistory);

  ++MergedNr;
}
 
