#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "core_allvars_grid.h"
#include "core_proto_grid.h"

int32_t update_grid_properties(int32_t filenr)
{

  int32_t snapshot_idx, grid_num_idx;
  int64_t grid_position, gal_idx, good_gals = 0, bad_gals = 0;
  float fesc_local, Ngamma_HI, Ngamma_HeI, Ngamma_HeII;

  for (snapshot_idx = LowSnap; snapshot_idx < HighSnap + 1; ++snapshot_idx)
  { 
    grid_num_idx = snapshot_idx - LowSnap; // The grid indexing goes from 0 to NumGrids.
    for (gal_idx = 0; gal_idx < NtotGals; ++gal_idx)
    {
      grid_position = GalGrid[gal_idx].History[snapshot_idx];
      if (grid_position == -1)
      {
        continue;
      }

      if (grid_position < 0 || grid_position > Grid->NumCellsTotal)
      {
        fprintf(stderr, "The grid index must be between 0 and the GridSize cubed %d cubed = %ld.  The grid index for Galaxy %ld is %ld.\n", GridSize, (long)Grid->NumCellsTotal, (long)gal_idx, (long)grid_position);
        return EXIT_FAILURE;
      }

      if ((GalGrid[gal_idx].StellarMass[snapshot_idx] > 0.0) & (GalGrid[gal_idx].SFR[snapshot_idx] > 0.0) & (GalGrid[gal_idx].CentralGalaxyMass[snapshot_idx] > 0.0) & (GalGrid[gal_idx].LenHistory[snapshot_idx] > HaloPartCut)) // Apply some requirements for the galaxy to be included.
      {
        ++good_gals;

        if (fescPrescription == 4)
        {
          update_quasar_tracking(gal_idx, snapshot_idx); 
        }

        Grid->GridProperties[grid_num_idx].SFR[grid_position] += GalGrid[gal_idx].SFR[snapshot_idx];
        Grid->GridProperties[grid_num_idx].StellarMass[grid_position] += GalGrid[gal_idx].StellarMass[snapshot_idx];  

        if (PhotonPrescription == 1)
        {
          calculate_photons(GalGrid[gal_idx].SFR[snapshot_idx], GalGrid[gal_idx].Z[snapshot_idx], &Ngamma_HI, &Ngamma_HeI, &Ngamma_HeII); // Base number of ionizing photons
          fesc_local = calculate_fesc(gal_idx, snapshot_idx, filenr);
        }

        Grid->GridProperties[grid_num_idx].Nion_HI[grid_position] += pow(10, Ngamma_HI - 50.0)*fesc_local; // We keep these in units of 10^50 photons/s.
        if (pow(10, Ngamma_HI - 50.0) * fesc_local < 0.0)
        {
          fprintf(stderr, "For galaxy %ld, the number of HI ionizing photons is %.4f\n", (long)gal_idx, pow(10, Ngamma_HI - 50.0) * fesc_local); 
          return EXIT_FAILURE;
        }
        if (Grid->GridProperties[grid_num_idx].Nion_HI[grid_position] < 0.0 || Grid->GridProperties[grid_num_idx].Nion_HI[grid_position] > 1e100)
        {
          fprintf(stderr, "For galaxy %ld, cell %ld now has an error number of photons. This number is %.4f e50 photons/s\n", (long)gal_idx, (long)grid_position, Grid->GridProperties[grid_num_idx].Nion_HI[grid_position]);
          return EXIT_FAILURE;
        }
        Grid->GridProperties[grid_num_idx].Nion_HeI[grid_position] += pow(10, Ngamma_HeI - 50.0)*fesc_local;
        Grid->GridProperties[grid_num_idx].Nion_HeII[grid_position] += pow(10, Ngamma_HeII - 50.0)*fesc_local;
      
        ++Grid->GridProperties[snapshot_idx].GalCount[grid_position]; 
        
      }
      else
      {
        ++bad_gals;
      }
    } // Galaxy loop.
  } // Snapshot loop.
  
  return EXIT_SUCCESS;
}

// Takes the number of halos loaded in memory (totNHalos) and maps the properties onto the grid. //
/*
void update_grid_nion_halo(int GridNr) // Calculates number of ionizing photons using the halos. 
{
  int i;
  double evolve_time;

  printf("Calculating number of photons using halo-based prescription.\n");
  
 
  if (GridNr == NGrid - 1)
	  evolve_time = (Age[ListOutputGrid[GridNr]-1] - Age[ListOutputGrid[GridNr]]) * UnitTime_in_Megayears;
  else
	  evolve_time = (Age[ListOutputGrid[GridNr+1]] - Age[ListOutputGrid[GridNr]]) * UnitTime_in_Megayears;

  for (i = 0; i < CUBE(GridSize); ++i)
  {
    // Using Illiev 2012. 
//    UnitConversion = SOLAR_MASS/0.53/PROTONMASS/SEC_PER_MEGAYEAR;
    Grid[i].Nion_HI = Grid[i].HaloMass*1e10/Hubble_h*SOLAR_MASS*SourceEfficiency*BaryonFrac/0.53/PROTONMASS/evolve_time/SEC_PER_MEGAYEAR;
  }
  
}
*/
void count_grid_properties(struct GRID_STRUCT *count_grid) // Count number of galaxies/halos in the grid.
{

  int32_t snapshot_idx, grid_num_idx;

  for (snapshot_idx = LowSnap; snapshot_idx < HighSnap + 1; ++snapshot_idx)
  {
    int64_t GlobalGalCount = 0, SourcesCount = 0, cell_idx;
    float totPhotons_HI = 0, totPhotons_HeI = 0, totPhotons_HeII = 0;

    grid_num_idx = snapshot_idx - LowSnap; // The grid indexing goes from 0 to NumGrids.
    for (cell_idx = 0; cell_idx < count_grid->NumCellsTotal; ++cell_idx)
    {
      GlobalGalCount += count_grid->GridProperties[grid_num_idx].GalCount[cell_idx];
  
      totPhotons_HI += count_grid->GridProperties[grid_num_idx].Nion_HI[cell_idx];
      totPhotons_HeI += count_grid->GridProperties[grid_num_idx].Nion_HeI[cell_idx];
      totPhotons_HeII += count_grid->GridProperties[grid_num_idx].Nion_HeII[cell_idx];

      if (count_grid->GridProperties[grid_num_idx].Nion_HI[cell_idx] > 0.0)
      {
        ++SourcesCount; 
      }
    }

    printf("At redshift %.3f (Snapshot %d) there was %ld galaxies and [%.4e, %.4e, %.4e]e50 {HI, HeI, HeII}  ionizing Photons emitted per second ([%.4e %.4e %.4e]e50 s^-1 Mpc^-3), spread across %ld cells (%.4f of the total cells).\n", ZZ[snapshot_idx], snapshot_idx, GlobalGalCount, totPhotons_HI, totPhotons_HeI, totPhotons_HeII, totPhotons_HI / pow(BoxSize/Hubble_h,3), totPhotons_HeI / pow(BoxSize/Hubble_h, 3), totPhotons_HeII / pow(BoxSize/Hubble_h,3), (long)SourcesCount, (double)SourcesCount / (double)count_grid->NumCellsTotal); 

  } // Snapshot loop.

}

/*
void normalize_photon(int GridNr)
{

  int i;
  double totPhotons_HI_normgrid = 0, totPhotons_HI = 0;
  char buf[MAXLEN], tag[MAXLEN];
  double *NormGrid;
  FILE *load_fd; 

  if ((NGrid - 1 - GridNr) < 10)
    snprintf(tag, MAXLEN, "%s", "00");
  else if ((NGrid - 1 - GridNr) < 100 && (NGrid - 1 - GridNr) >= 10)
    snprintf(tag, MAXLEN, "%s", "0");
  else
    snprintf(tag, MAXLEN, "%s", "");


  NormGrid = malloc(CUBE(GridSize) * sizeof(double));

  snprintf(buf, MAXLEN, "/lustre/projects/p004_swin/jseiler/SAGE_output/1024/grid/January_input/Galaxies_noreion_z5.000_fesc0.15_nionHI_%s%d", tag, (NGrid-1) - GridNr);
  if(!(load_fd = fopen(buf, "r")))
  {
    printf("can't open file `%s'\n", buf);
    exit(0); 
  }

  for (i = 0; i < CUBE(GridSize); ++i)
  {
    fread(&NormGrid[i], 1, sizeof(double), load_fd);
    totPhotons_HI_normgrid += NormGrid[i];
    totPhotons_HI += Grid[i].Nion_HI; 
  }

  double ratio = totPhotons_HI_normgrid/totPhotons_HI;

  for (i = 0; i < CUBE(GridSize); ++i)
  {  
    Grid[i].Nion_HI = Grid[i].Nion_HI * ratio;
  }

  printf("Total photons from normalization grid is %.4e compared to the %.4e photons from fescPrescription %d giving a ratio of %.4e.\n", totPhotons_HI_normgrid, totPhotons_HI, fescPrescription, ratio); 

}

void normalize_slope_photons(int GridNr)
{
  
   int i;
   double totPhotons = 0, targetPhotons = pow(10,-0.8*ZZ[ListOutputGrid[GridNr]] + 64.3); 


   for (i = 0; i < CUBE(GridSize); ++i)
   {
      totPhotons += Grid[i].Nion_HI;
   }

   double ratio = targetPhotons/totPhotons;

   for (i = 0; i < CUBE(GridSize); ++i)
   {
     Grid[i].Nion_HI = Grid[i].Nion_HI * ratio;
   }

   printf("The total number of photons in the grid is %.4e.  We only want %.4e photons for redshift %.4e giving a ratio of %.4e.\n", totPhotons, targetPhotons, ZZ[ListOutputGrid[GridNr]], ratio);
}
*/

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

  if (SFR == 0)
  {
    *Ngamma_HI = 0;
    *Ngamma_HeI = 0;
    *Ngamma_HeII = 0; 
  }   
  else if (Z < 0.0025) // 11
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

// Here we calculate the escape fraction for a specific galaxy. 
// The user defines a value (fescPrescription) that determines how to calculate fesc.
// 0: Constant Escape fraction.  User defines the exact value.
// 1: Scaling with Halo Mass using the functional form defined by Kimm et al (2016).
// 2: Power Law as a function of halo mass. The user defines the smallest and largest halo mass in the simulation in addition to the escape fractions at these masses.
// 3: Linear relationship as a function of the ejected fraction of a galaxy.  User defines the escape fraction for ejected fractions of 0 and 1.	
//
// The values of alpha/beta are determined within the 'core_init.c' module.
//
// INPUT: The galaxy index (p). 
// 	: The snapshot index (i)
// 	: The filenumber of the galaxy (filenr); useful for debugging. 
//
// OUTPUT: The escape fraction for the galaxy. 

float calculate_fesc(int p, int i, int filenr)
{

  float fesc_local, halomass, ejectedfraction;

  halomass = GalGrid[p].CentralGalaxyMass[i];
  ejectedfraction = GalGrid[p].EjectedFraction[i];
  
  if (fescPrescription == 0) 
    fesc_local = fesc;
  else if (fescPrescription == 1)
    fesc_local = pow(10,1.0 - 0.2*log10(halomass * 1.0e10 / Hubble_h));
  else if (fescPrescription == 2)
    fesc_local = alpha * pow((halomass * 1.0e10 / Hubble_h), beta); 
  else if (fescPrescription == 3)	
    fesc_local = alpha * ejectedfraction + beta; 
	
  if (fesc_local > 1.0)
  {
  //  fprintf(stderr, "Had fesc_local = %.4f for galaxy %d in file %d with halo mass %.4e (log Msun), Stellar Mass %.4e (log Msun), SFR %.4e (log Msun yr^-1) and Ejected Fraction %.4e\n", fesc_local, p, filenr, log10(halomass * 1.0e10 / Hubble_h), log10(halomass * 1.0e10 / Hubble_h), log10(SFR), ejectedfraction);
    fesc_local = 1.0;
  
  }
   
  if (fesc_local < 0.0)
  {
    fprintf(stderr, "Had fesc_local = %.4f for galaxy %d in file %d with halo mass %.4e (log Msun), Stellar Mass %.4e (log Msun), SFR %.4e (log Msun yr^-1) and Ejected Fraction %.4e\n", fesc_local, p, filenr, log10(GalGrid[p].CentralGalaxyMass[i] * 1.0e10 / Hubble_h), log10(GalGrid[p].StellarMass[i] * 1.0e10 / Hubble_h), log10(GalGrid[p].SFR[i]), GalGrid[p].EjectedFraction[i]);
    fesc_local = 0.0; 
  }

  return fesc_local;

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

#ifdef MPI
struct GRID_STRUCT *MPI_sum_grids(void)
{

  int32_t status, grid_num_idx;
  struct GRID_STRUCT *master_grid;
  
  master_grid = malloc(sizeof(struct GRID_STRUCT));

  if (ThisTask == 0)
  {
    printf("Trying to initialize the master grid\n");
    status = init_grid(master_grid);
    if (status == EXIT_FAILURE)
    { 
      return NULL;
    }
  }

  for (grid_num_idx = 0; grid_num_idx < master_grid->NumGrids; ++grid_num_idx)
  {
    if (ThisTask == 0)
    {
      printf("Reducing grid %d\n", grid_num_idx);
    }
    MPI_Reduce(Grid->GridProperties[grid_num_idx].SFR, master_grid->GridProperties[grid_num_idx].SFR, master_grid->NumCellsTotal, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
    MPI_Reduce(Grid->GridProperties[grid_num_idx].Nion_HI, master_grid->GridProperties[grid_num_idx].Nion_HI, master_grid->NumCellsTotal, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
    MPI_Reduce(Grid->GridProperties[grid_num_idx].Nion_HeI, master_grid->GridProperties[grid_num_idx].Nion_HeI, master_grid->NumCellsTotal, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 
    MPI_Reduce(Grid->GridProperties[grid_num_idx].Nion_HeII, master_grid->GridProperties[grid_num_idx].Nion_HeII, master_grid->NumCellsTotal, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD); 

    MPI_Reduce(Grid->GridProperties[grid_num_idx].GalCount, master_grid->GridProperties[grid_num_idx].GalCount, master_grid->NumCellsTotal, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD); 
    
  }

  return master_grid;
}
#endif

int32_t update_quasar_tracking(int64_t gal_idx, int32_t snapshot_idx)
{

  float dt, substep_weight, time_into_snapshot, fraction_into_snapshot; 

  if (GalGrid[gal_idx].QuasarActivity[snapshot_idx] == 1) // A quasar has gone off during this snapshot, time to update properties.
  {
    ++QuasarActivityToggle[gal_idx]; // Note, we plus one because we want to be able to handle the case of a quasar going off when the galaxy still is being boosted. 
    QuasarSnapshot[gal_idx] = snapshot_idx; 
    TargetQuasarTime[gal_idx] = GalGrid[gal_idx].DynamicalTime[snapshot_idx];
    QuasarActivitySubstep[gal_idx] = GalGrid[gal_idx].QuasarSubstep[snapshot_idx];
    QuasarBoostActiveTime[gal_idx] = 0.0;
  }

  if (QuasarActivityToggle[gal_idx] > 0) // This galaxy is having its escape fraction boosted, check to see if we need to turn it off.
  {

    dt = (Age[snapshot_idx - 1] - Age[snapshot_idx]) * UnitTime_in_Megayears;
    if (QuasarSnapshot[gal_idx] == snapshot_idx && QuasarActivityToggle[gal_idx] == 1) // If this quasar is due to a quasar going off during this snapshot and the galaxy is NOT under the influence from a previous quasar event then we need to weight the fraction of time the photons are boosted by the substep the quasar went off in. 
    {
      substep_weight = (STEPS - QuasarActivitySubstep[gal_idx]) / STEPS;
    }
    else
    {
      substep_weight = 1.0;
    } 

    QuasarBoostActiveTime[gal_idx] += dt * substep_weight;
    QuasarFractionalPhoton[gal_idx] = substep_weight; // If the quasar turned on part-way through the snapshot, we boost the photons for the remaining time during the snapshot.

    if (QuasarBoostActiveTime[gal_idx] >= TargetQuasarTime[gal_idx]) // The boosted quasar time needs to be turned off.
    {
      time_into_snapshot = TargetQuasarTime[gal_idx] - (QuasarBoostActiveTime[gal_idx] - dt); // How much extra time into the snapshot does the quasar need to go to reach its target?
      fraction_into_snapshot = time_into_snapshot / dt; // Then what fraction of the snapshot time will this be?
      QuasarFractionalPhoton[gal_idx] = fraction_into_snapshot; 

      // Reset toggles and trackers. //
      --QuasarActivityToggle[gal_idx];
      QuasarSnapshot[gal_idx] = -1;
      TargetQuasarTime[gal_idx] = 0.0;
      QuasarBoostActiveTime[gal_idx] = 0.0;
      QuasarActivitySubstep[gal_idx] = -1;
    }

  } 

  return EXIT_SUCCESS;

}
