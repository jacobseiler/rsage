/*
This file contains function that deal with the self_consistent coupling of the RSAGE model.
Specifically, it handles the updating of the `self_consistent_grid` (commonly referred to as the Nion grid)
that tracks the number of HI ionizing photons for the specified snapshot.

This ionizing photon grid is then fed into `cifog` to determine the ionization regions.

Author: Jacob Seiler
Version: 0.0.1
*/

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

#include "../core_allvars.h"
#include "../core_proto.h"
#include "selfcon_grid.h"

#ifdef MPI
#include <mpi.h>
struct SELFCON_GRID_STRUCT *MPI_sum_grids(void);
#endif
// Local Variables //

#define STARBURSTSTEP 0.1 // This is the step size for the Starburst99 data (in Myr). 

// Local Proto-Types //

int32_t malloc_selfcon_grid(struct SELFCON_GRID_STRUCT *grid);
int32_t write_selfcon_grid(struct SELFCON_GRID_STRUCT *grid_towrite);

// External Functions //

/*
Prints out some useful information regarding `fescPrescription` and initializes the Nion grid.

Parameters
----------

None.  All variables used/adjusted are global.

Returns
----------

None. All variables adjusted are global.

Pointer Updates
----------

None.

Units  
----------

None.
*/

int32_t init_selfcon_grid(void)
{

  int32_t status;

#ifdef MPI
  if (ThisTask == 0)
#endif
  printf("Initializing the grids for the self_consistent run.\n");
  
  switch(fescPrescription)
  {

    case 0:
      if (ThisTask == 0)
      {
        printf("\n\nUsing a constant escape fraction of %.4f\n", beta); 
      }
      break;

    case 1:
      if (ThisTask == 0)
      {
        printf("\n\nUsing an fesc prescription that scales with the fraction of ejected mass in the galaxy.\nThis takes the form A*fej + B with A = %.4e and B = %.4e\n", alpha, beta);
      }
      break;

    case 2:
      if (ThisTask == 0)
      {
        printf("\n\nUsing an fesc prescription that depends upon quasar activity.\n");
        printf("\nFor a galaxy that had a quasar event within %.2f dynamical times go the escape fraction will be %.2f.  Otherwise it will have a constant value of %.2f\n", N_dyntime, quasar_boosted, quasar_baseline);
      }
      break;

    case 3:
      if (ThisTask == 0)
      {
        printf("\n\nUsing Anne's functional form for an escape fraction that decreases for increasing halo mass.\n");
        printf("MH_low = %.4e\tMH_high = %.4e\tfesc_low = %.4f\tfesc_high =  %.4f.\n", MH_low, MH_high, fesc_low, fesc_high);
      }
      XASSERT(fesc_low > fesc_high, "Input file contain fesc_low = %.2f and fesc_high = %.2f. For this prescription (fescPrescription == 3), we require fesc_low > fesc_high\n", fesc_low, fesc_high);

      break;

    case 4:
      if (ThisTask == 0)
      {
        printf("\n\nUsing Anne's functional form for an escape fraction that increases for increasing halo mass.\n");
        printf("MH_low = %.4e\tMH_high = %.4e\tfesc_low = %.4f\tfesc_high =  %.4f.\n", MH_low, MH_high, fesc_low, fesc_high);
      }
      XASSERT(fesc_low < fesc_high, "Input file contain fesc_low = %.2f and fesc_high = %.2f. For this prescription (fescPrescription == 4), we require fesc_low < fesc_high\n", fesc_low, fesc_high);
      break;

    case 5:
      if (ThisTask == 0)
      {
        printf("\n\nUsing an fesc prescription that depends upon the SFR of the galaxy.\n");
      }
      break;

    default:
      printf("\n\nOnly escape fraction prescriptions 0, 1, 2, 3, 4 or 5 are permitted.\n");
      return EXIT_FAILURE;
  }

  SelfConGrid = malloc(sizeof(*(SelfConGrid)));
  if (SelfConGrid == NULL)
  {
    fprintf(stderr, "Could not allocate memory for the high level self_consistent grid structure\n");
    return EXIT_FAILURE;
  }

  status = malloc_selfcon_grid(SelfConGrid);
  if (status != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;

}

/*
Frees the Nion grid using calls to `myfree()` to keep track of memory usage.

Parameters
----------

*grid_to_free: struct SELFCON_GRID_STRUCT pointer. See `core_allvars.h` for full struct contents.
  Pointer to the Nion grid that we are freeing.

Returns
----------

EXIT_SUCCESS.

Pointer Updates
----------

None.

Units  
----------

None.
*/

int32_t free_selfcon_grid(struct SELFCON_GRID_STRUCT *grid_to_free)
{ 

  // First free the inner arrays.
  myfree(grid_to_free->Nion_HI, sizeof(*(grid_to_free->Nion_HI)) * grid_to_free->NumCellsTotal);
  myfree(grid_to_free->GalCount, sizeof(*(grid_to_free->GalCount)) * grid_to_free->NumCellsTotal);
  
  free(grid_to_free);

  return EXIT_SUCCESS;

}

/*
Updates an Nion grid cell for the current galaxy.  

This funtion calculates the number of ionizing photons emitted, and the number that escape (using the galaxy specific of fesc).

Also saves the photon properties (e.g., value of fesc/Ngamma) of the galaxy to an ASCII file.  This file is either opened if it's
yet to be, otherwise the data is appended to the already opened file.

Parameters
----------

*g: struct GALAXY pointer.  See `core_allvars.h` for the struct architecture.
  Pointer to the galaxy that we are calculating ionizing photons for.

grid_idx: Integer.
  The 1-Dimensional grid index that corresponds to the galaxy position. 
  See `determine_1D_idx()` in `model_misc.c` for full details.

snapshot: Integer. 
  The snapshot number we are calculating for. 

Returns
----------

EXIT_SUCCESS or EXIT_FAILURE.
  If `determine_nion()` returns EXIT_FAILURE, also returns EXIT_FAILURE.
  If `determine_fesc()` returns EXIT_FAILURE, also returns EXIT_FAILURE.
  If the photon properties file cannot be opened (e.g., invalid pth), EXIT_FAILURE is returned.
  If the photon properties file is not opened prior to `fprintf` attempt, EXIT_FAILURE is returned. 
 
  Otherwise EXIT_SUCCESS is returned.

Pointer Updates
----------

None.

Units  
----------

The number of ionizing photons (Ngamma) is kept in units of 1.0e50 photons/s.
*/

int32_t update_selfcon_grid(struct GALAXY *g, int32_t grid_idx, int32_t snapshot)
{

  int32_t status;
  float Ngamma_HI, Ngamma_HeI, Ngamma_HeII, fesc_local;

  status = determine_nion(g, snapshot, &Ngamma_HI, &Ngamma_HeI, &Ngamma_HeII);
  /*
  if (Ngamma_HI < 1e-6 && g->GridStellarMass[snapshot] > 0.0 && PhotonPrescription == 1)
  {
    return EXIT_FAILURE; 
  }
  */
  if (status != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }

  status = determine_fesc(g, snapshot, &fesc_local);
  if (status != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }

  if (Ngamma_HI > 0.0)
  {
    SelfConGrid->Nion_HI[grid_idx] += Ngamma_HI * fesc_local; // Ngamma_HI is already in units of 1.0e50 photons/s. 
    SelfConGrid->Nion_HI_Total += Ngamma_HI * fesc_local; 
  }

  ++(SelfConGrid->GalCount[grid_idx]);

  return EXIT_SUCCESS;
}

/*
Wrapper function that writes out the Nion grid.
**NOTE**: The freeing of `SelfConGrid` is performed by `sage_cleanup()`.
This function is called at the end of the script, before final cleanup.

If RSAGE has been built using MPI, then each rank will have its own Nion grid.  In this case, this function calls `MPI_sum_grid()`
to sum all the grids onto rank 0 which is then subsequently written and freed.
 
Parameters
----------

None. All the grids are kept global.

Returns
----------

EXIT_SUCCESS or EXIT_FAILURE.
  All non-rank 0 processors return EXIT_SUCCESS.  
  `MPI_sum_grids()` returns a NULL pointer for all non-rank 0 processors.  However if rank 0 also returns NULL, EXIT_FAILURE is returned by this function.
  If `write_selfcon_grid()` returns EXIT_FAILURE, also returns EXIT_FAILURE.
 
  Otherwise EXIT_SUCCESS is returned.

Pointer Updates
----------

None.

Units  
----------

None.
*/

int32_t save_selfcon_grid(void)
{

  int32_t status;

#ifdef MPI

  struct SELFCON_GRID_STRUCT *master_grid;

  master_grid = MPI_sum_grids(); 
  if (ThisTask == 0 && master_grid == NULL)
  {
    return EXIT_FAILURE;
  }

  if (ThisTask == 0)
  {
    status = write_selfcon_grid(master_grid);
  }
  else
    status = EXIT_SUCCESS;
#else
  status = write_selfcon_grid(SelfConGrid);
#endif

  if (status != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }

#ifdef MPI
  if (ThisTask == 0)
    free_selfcon_grid(master_grid);    
#endif

  return EXIT_SUCCESS;
}

/*
Sets the self-consistent grid to zeros.

Parameters
----------

*my_grid: struct SELFCON_GRID_STRUCT pointer. See `core_allvars.h` for full struct architecture.
  Nion grid that is being allocated.
  **IMPORTANT** This pointer must be allocated before being passed to this function.
  If `*my_grid` is passed as NULL, EXIT_FAILURE is returned.

Returns
----------

EXIT_SUCCESS or EXIT_FAILURE.
  If `*my_grid` is passed as an unallocated NULL, EXIT_FAILURE is returned. 
  If memory cannot be allocated for `Nion_HI` or `GalCount`, EXIT_FAILURE is returned.
   
  Otherwise EXIT_SUCCESS is returned.

Pointer Updates
----------

The inner arrays of `*my_grid` are allocated memory and initialized to 0.
See `core_allvars.h` for full architecture of `SELFCON_GRID_STRUCT`.

Units  
----------

None.
*/

int32_t zero_selfcon_grid(struct SELFCON_GRID_STRUCT *my_grid)
{
  
  int32_t cell_idx;

  if (my_grid == NULL)
  {
    fprintf(stderr, "`zero_selfcon_grid` was called with a SELFCON_GRID_STRUCT pointer that has not been initialized\n");
    return EXIT_FAILURE;
  }

  my_grid->Nion_HI_Total = 0.0;

  for (cell_idx = 0; cell_idx < my_grid->NumCellsTotal; ++cell_idx)
  {
    my_grid->Nion_HI[cell_idx] = 0.0;
    my_grid->GalCount[cell_idx] = 0;
  }

  return EXIT_SUCCESS;
}

// Local Functions //

/*
Allocates memory and initializes values for the Nion grid.
 
Parameters
----------

*my_grid: struct SELFCON_GRID_STRUCT pointer. See `core_allvars.h` for full struct architecture.
  Nion grid that is being allocated.
  **IMPORTANT** This pointer must be allocated before being passed to this function.
  If `*my_grid` is passed as NULL, EXIT_FAILURE is returned.

Returns
----------

EXIT_SUCCESS or EXIT_FAILURE.
  If `*my_grid` is passed as an unallocated NULL, EXIT_FAILURE is returned. 
  If memory cannot be allocated for `Nion_HI` or `GalCount`, EXIT_FAILURE is returned.
   
  Otherwise EXIT_SUCCESS is returned.

Pointer Updates
----------

The inner arrays of `*my_grid` are allocated memory and initialized to 0.
See `core_allvars.h` for full architecture of `SELFCON_GRID_STRUCT`.

Units  
----------

None.
*/

int32_t malloc_selfcon_grid(struct SELFCON_GRID_STRUCT *my_grid)
{

#define ALLOCATE_GRID_MEMORY(name, length)  \
{                                           \
  name = mycalloc(length, sizeof(*(name))); \
  if (name == NULL)                         \
  {                                         \
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate"#name".\n", sizeof(*(name)* length)); \
    return EXIT_FAILURE;                    \
  }                                         \
}

  if (my_grid == NULL)
  {
    fprintf(stderr, "`malloc_selfcon_grid` was called with a SELFCON_GRID_STRUCT pointer that has not been initialized\n");
    return EXIT_FAILURE;
  }

  my_grid->GridSize = GridSize;
  my_grid->NumCellsTotal = CUBE(GridSize);

  ALLOCATE_GRID_MEMORY(my_grid->Nion_HI, my_grid->NumCellsTotal);
  ALLOCATE_GRID_MEMORY(my_grid->GalCount, my_grid->NumCellsTotal);
  
  return EXIT_SUCCESS;

#undef ALLOCATE_GRID_MEMORY

}

/*
Calculates the ionizing photon rate for a given galaxy at a specified snapshot.

Depending on the value of `PhotonPrescription` specified, the prescription to calculate this will change.

0: Assumes a continuous SFR over the Snapshot and assigns an ionizing photon rate proportional to log(SFR)
1: Uses the results of STARBURST99 to assign an ionizing photon rate that depends upon the age of previous SF episodes. 

**Important** See Units. 

Parameters
----------

g: struct GALAXY pointer.  See `core_allvars.h` for the struct architecture.
  Pointer to the galaxy that we are calculating ionizing photons for.

snapshot: Integer. 
  The snapshot number we are calculating for. 

*Ngamma_HI, *Ngamma_HeI, *Ngamma_HeII: Float Pointers.
  Pointers that will store the number of hydrogen, helium and helium II ionizing photons. 

Returns
----------

EXIT_SUCCESS or EXIT_FAILURE.
  For `PhotonPrescription == 0`, only EXIT_SUCCESS can be returned (no fail conditions).
  For `PhotonPrescription == 1`, if the number of ionizing photons is negative, EXIT_FAILURE is returned.  Otherwise EXIT_SUCCESS is returned. 
  For any other `PhotonPrescription`, EXIT_FAILURE is returned.  

Pointer Updates
----------

*Ngamma_HI, *Ngamma_HeI, *Ngamma_HeII.

Units  
----------

The ionizing photon rate is returned in units of 1.0e50 Photons/s.
Star formation rate is converted from internal code units to Msun/yr.
Metallicity is in absolute metallicity (not solar).
*/

int32_t determine_nion(struct GALAXY *g, int32_t snapshot, float *Ngamma_HI, float *Ngamma_HeI, float *Ngamma_HeII)
{

  switch (PhotonPrescription)
  {
    case 0: ;
      
      const double SFR_CONVERSION = UnitMass_in_g/UnitTime_in_s*SEC_PER_YEAR/SOLAR_MASS/STEPS; // Conversion from the SFR over one snapshot to Msun/yr. 

      double SFR = g->GridSFR[snapshot] * SFR_CONVERSION;
      double Z = g->GridZ[snapshot];

      if (SFR == 0)
      {
        *Ngamma_HI = 0;
        *Ngamma_HeI = 0;
        *Ngamma_HeII = 0;
        return EXIT_SUCCESS;
      }
      else if (Z < 0.0025) // 11
      {
        *Ngamma_HI = log10(SFR) + 53.154;
        *Ngamma_HeI = log10(SFR) + 52.727;
        *Ngamma_HeII = log10(SFR) + 48.941;
      }
      else if (Z >= 0.0025 && Z < 0.006) // 12
      {
        *Ngamma_HI = log10(SFR) + 53.090;
        *Ngamma_HeI = log10(SFR) + 52.583;
        *Ngamma_HeII = log10(SFR) + 49.411;
      }
      else if (Z>= 0.006 && Z < 0.014) // 13
      {
        *Ngamma_HI = log10(SFR) + 53.048;
        *Ngamma_HeI = log10(SFR) + 52.481;
        *Ngamma_HeII = log10(SFR) + 49.254;
      }
      else if (Z >= 0.014 && Z < 0.030) // 14
      {
        *Ngamma_HI = log10(SFR) + 52.966;
        *Ngamma_HeI = log10(SFR) + 52.319;
        *Ngamma_HeII = log10(SFR) + 48.596;
      }
      else // 15
      {
        *Ngamma_HI = log10(SFR) + 52.941;
        *Ngamma_HeI = log10(SFR) + 52.052;
        *Ngamma_HeII = log10(SFR) + 47.939;
      }

      *Ngamma_HI = exp10(*Ngamma_HI - 50.0);
      return EXIT_SUCCESS;

    case 1: ;

      double t;
      int32_t i, lookup_idx;
      const double lookuptable_mass = 1.0e-4*Hubble_h;

      *Ngamma_HI = 0;
      *Ngamma_HeI = 0;
      *Ngamma_HeII = 0;

      for (i = 0; i < StellarTracking_Len - 1; ++i)
      {
        if (g->Stellar_Stars[i] < 1e-10)
          continue;

        t = (i + 1) * TimeResolutionStellar; // (i + 1) because 0th entry will be at TimeResolutionSN.
        lookup_idx = (t / 0.1); // Find the index in the lookup table. 
        
        *Ngamma_HI += exp10(log10(g->Stellar_Stars[i] / lookuptable_mass) + stars_Ngamma[lookup_idx] - 50.0);        
      }

      // The units of Ngamma are 1.0e50 photons/s.  Hence a negative value is not allowed.
      if (*Ngamma_HI < 0.0)
      {
        fprintf(stderr, "Got an NgammaHI value of %.4e.  This MUST be a positive value.\nPrinting out information for every element used to calculate Ngamma.\n", *Ngamma_HI);

        // Print out information for every element of the array so we can try identify the problem.       
        for (i = 0; i < StellarTracking_Len; ++i)
        {
          t = (i + 1) * TimeResolutionStellar; // (i + 1) because 0th entry will be at TimeResolutionSN.
          lookup_idx = (t / 0.1); // Find the index in the lookup table. 

          double total = 0.0;
          total += exp10(log10(g->Stellar_Stars[i] / lookuptable_mass) + stars_Ngamma[lookup_idx] - 50.0);  
          printf("t %.4f\tlookup_idx %d\tg->Stellar_Stars[i] %.4e\tstars_Ngamma[lookup_idx] %.4e\tRunningTotal %.4e\n", t, lookup_idx, g->Stellar_Stars[i], stars_Ngamma[lookup_idx], total); 
        }
        return EXIT_FAILURE;
      }

      return EXIT_SUCCESS;

    default:
      fprintf(stderr, "The specified PhotonPrescription value is not valid.");
      return EXIT_FAILURE;
  }
}

/*
Calculates the value of fesc for this galaxy for a specific snapshot.

Depending on the value of `PhotonPrescription` specified, the prescription to calculate this will change.

0: Assigns a single value of fesc to all galaxies, regardless of properties.
   Value is given by `fesc`. 
1: DEPRECATED.
2: Power law as a function of halo mass, fesc = alpha*MH^beta. 
   The values of `alpha` and `beta` are given by specifying two fixed points: (MH_low, fesc_low) and (MH_high, fesc_high).
   These fixed points are specified in units of Msun; see `calculate_fesc_constants()` for exact equation. 
3: Linear relationship as a function of the fraction of ejected baryons in the galaxy, fesc = alpha*fej + beta.
   The values of `alpha` and `beta` are specified directly by the variables `alpha` and `beta` in the .ini file.
4: The value of fesc is boosted by recent quasar activity. Each galaxy has a baselines fesc of `quasar_baseline`.
   Following a quasar event that ejects all gas from a galaxy, the galaxy is given an fesc value `quasar_boosted` for `N_dyntime` dynamical times.
   After this, it falls back to `quasar_baseline`.
5,6: fesc either increases/decreases as a function of Halo Mass using Anne's specified functional forms.
   Ths fixed points are specified by giving (MH_low, fesc_low) and (MH_high, fesc_high).
7: Same as 3 except the relationship is a power law of the form fesc = alpha*fej^beta.
8: For galaxies with stellar mass between `fest_Mstar_low` and `fesc_Mstar_high` the value of fesc is set to `fesc_Mstar`.
   All other galaxies are set to fesc = `fesc_not_Mstar`.

Parameters
----------

g: struct GALAXY pointer.  See `core_allvars.h` for the struct architecture.
  Pointer to the galaxy that we are calculating value of fesc for.

snapshot: Integer. 
  The snapshot number we are calculating for. 

*fesc_local: Float Pointers.
  Pointer that will store the value of fesc for this galaxy at this snapshot. 

Returns
----------

EXIT_SUCCESS or EXIT_FAILURE.
  If the value of `fescPrescription` is any value other than [0, 4] , EXIT_FAILURE is returned.
  If `fesc_local` is assigned a value less than 0 or greater than 1.0 EXIT_FAILURE is returned.
  Otherwise, EXIT_SUCCESS is returned.

Pointer Updates
----------

*fesc_local.

Units  
----------

Halo masses are calculated in Msun.
Stellar masses are calculated in Msun.
Ejected fraction is the fraction of baryons (unitless).
*/

int32_t determine_fesc(struct GALAXY *g, int32_t snapshot, float *fesc_local)
{

  //float halomass = g->GridFoFMass[snapshot] * 1.0e10 / Hubble_h; // Halo Mass (of the background FoF) in Msun.
  float halomass = g->GridHaloMass[snapshot] * 1.0e10 / Hubble_h; // Halo Mass (of the host halo, NOT BACKGROUND FOF) in Msun.
  float quasarfrac = g->QuasarFractionalPhotons;  
  float ejectedfraction = 0.0; // Initialize to prevent warnings.

  const double SFR_CONVERSION = UnitMass_in_g/UnitTime_in_s*SEC_PER_YEAR/SOLAR_MASS/STEPS; // Conversion from the SFR over one snapshot to Msun/yr. 
  float sfr = g->GridSFR[snapshot] * SFR_CONVERSION; // Star formation rate in Msun yr^-1.

  switch (fescPrescription)
  {
    case 1:
      ejectedfraction = g->EjectedFraction[snapshot];
      break;

    default:
      break;
  }

  switch (fescPrescription)
  {
    case 0:
      *fesc_local = beta;
      break;

    case 1:
      *fesc_local = alpha * ejectedfraction + beta;
      if (*fesc_local < 0.0)
      {
        *fesc_local = 0.0;
      }

      if (*fesc_local > 1.0)
      {
        *fesc_local = 1.0;
      }    
      break;
    
    case 2:
      *fesc_local = quasar_baseline * (1 - quasarfrac)  + quasar_boosted * quasarfrac;
      break;
    
    case 3:
      *fesc_local = pow(fesc_low * (fesc_low/fesc_high),(-log10(halomass/MH_low)/log10(MH_high/MH_low)));
      if (*fesc_local > fesc_low)
      {
        *fesc_local = fesc_low;
      }
      break;

    case 4:
      *fesc_local = 1. - pow((1.-fesc_low) * ((1.-fesc_low)/(1.-fesc_high)),(-log10(halomass/MH_low)/log10(MH_high/MH_low)));
      if (*fesc_local < fesc_low)
      {
        *fesc_local = fesc_low;
      }
      break;

    case 5:
      if (sfr > 0.0)
      {
        *fesc_local = 1.0 / (1.0 + exp(-(log10(sfr)-1))); // -1 to have SFR = 1 corresponds to fesc = 0.5.
      }
      else
      {
        *fesc_local = 0.0;
      }
      break;

    default:
      fprintf(stderr, "The selected fescPrescription is not handled by the switch case in `determine_fesc` in `selfcon_grid.c`.\nPlease add it there.\n");
      return EXIT_FAILURE;

  }

  if (*fesc_local > 1.0 || *fesc_local < 0.0)
  {
    fprintf(stderr, "Had fesc_local = %.4f with halo mass %.4e (log Msun), Stellar Mass %.4e (log Msun), SFR %.4e (log Msun yr^-1), Ejected Fraction %.4e and QuasarFractionalPhotons %.4f\n", *fesc_local, log10(halomass * 1.0e10 / Hubble_h), log10(halomass * 1.0e10 / Hubble_h), log10(sfr), ejectedfraction, quasarfrac);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

#ifdef MPI

/*
If RSAGE is built with MPI, then each processor will have its own Nion grid.  
Since we wish to write out a single grid, we need to sum all the grids across processors.  
This function allocates a single master grid on rank 0 and then sums everything into that grid. 
 
Parameters
----------

None. All grids are kept global.

Returns
----------

All non-rank 0 processors return NULL.
 
If memory could not be allocated for the master_grid (see `malloc_selfcon_grid()`), NULL is returned for rank 0. 
Otherwise:
   *master_grid: struct SELFCON_GRID_STRUCT pointer.
    The summed master grid on rank 0. 

Pointer Updates
----------

None.

Units  
----------

None.
*/

struct SELFCON_GRID_STRUCT *MPI_sum_grids(void)
{

  int32_t status;
  struct SELFCON_GRID_STRUCT *master_grid;
  master_grid = malloc(sizeof(*(master_grid))); // Needs to be malloced for all Tasks for Reduce.

  if (ThisTask == 0)
  {
  
    status = malloc_selfcon_grid(master_grid);
    if (status != EXIT_SUCCESS)
    {
      fprintf(stderr, "Could not allocate memory when trying to sum the grids across tasks.\n");
      return NULL;
    } 
    printf("Reducing selfcon grid.\n");
  }

  MPI_Reduce(SelfConGrid->Nion_HI, master_grid->Nion_HI, SelfConGrid->NumCellsTotal, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&SelfConGrid->Nion_HI_Total, &master_grid->Nion_HI_Total, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (ThisTask == 0)
    return master_grid;
  else
    return NULL; 

}

#endif

/*
Function to write the Nion grid to disk as binary. The name of the file depends on the value of `fescPrescription` specified.
In general, the filename will capture the variables that were specified in the .ini file in addition to a tag to specify what `fescPrescription` was used. 

If RSAGE is built with MPI, this function will only be called by rank 0.
 
Parameters
----------

*grid_towrite: struct SELFCON_GRID_STRUCT pointer. See `core_allvars.h` for full struct architecture..
  Nion grid that is being written.

Returns
----------

EXIT_SUCCESS or EXIT_FAILURE.
  If the value of `fescPrescription` is any value other than [0, 4], EXIT_FAILURE is returned.
  If the output file could not be opened, EXIT_FAILURE is returned.
  If the file could not be written to, EXIT_FAILURE is returned. 

  Otherwise EXIT_SUCCESS is returned.

Pointer Updates
----------

None.

Units  
----------

None.
*/

int32_t write_selfcon_grid(struct SELFCON_GRID_STRUCT *grid_towrite)
{ 

  FILE* file_HI;
  char tag[MAX_STRING_LEN], fname_HI[MAX_STRING_LEN];
  int32_t nwritten;

  switch(fescPrescription)
  {
    case 0:
      snprintf(tag, MAX_STRING_LEN - 1, "fesc%.2f_HaloPartCut%d", beta, HaloPartCut);
      break;

    case 1:
      snprintf(tag, MAX_STRING_LEN - 1, "ejected_%.3f_%.3f_HaloPartCut%d", alpha, beta, HaloPartCut); 
      break;

    case 2:
      snprintf(tag, MAX_STRING_LEN - 1, "quasar_%.2f_%.2f_%.2f_HaloPartCut%d", quasar_baseline, quasar_boosted, N_dyntime, HaloPartCut);
      break;

    case 3:
    case 4:
      snprintf(tag, MAX_STRING_LEN - 1, "AnneMH_%.3e_%.2f_%.3e_%.2f_HaloPartCut%d", MH_low, fesc_low, MH_high, fesc_high, HaloPartCut);      
      break;

    case 5:
      snprintf(tag, MAX_STRING_LEN - 1, "SFR_HaloPartCut%d", HaloPartCut);      
      break;

    default:
      fprintf(stderr, "The selected fescPrescription is not handled by the switch case in `save_selfcon_grid` in `selfcon_grid.c`.\nPlease add it there.\n");
      return EXIT_FAILURE;

  }

  snprintf(fname_HI, MAX_STRING_LEN, "%s/%s_%s_nionHI_%03d", GridOutputDir, FileNameGalaxies, tag, ReionSnap); 

  file_HI = fopen(fname_HI, "wb");
  if (file_HI == NULL)
  {
    fprintf(stderr, "Could not open file %s.\n", fname_HI);
    return EXIT_FAILURE;
  }

  nwritten = myfwrite(grid_towrite->Nion_HI, sizeof(*(grid_towrite->Nion_HI)) * grid_towrite->NumCellsTotal, 1, file_HI);
  if (nwritten != 1)
  {
    fprintf(stderr, "Could not write 1 element of size %zu to file %s, wrote %d instead.\n", sizeof(*(grid_towrite->Nion_HI)) * grid_towrite->NumCellsTotal, fname_HI, nwritten); 
    return EXIT_FAILURE;
  }
  fclose(file_HI);

  printf("Successfully wrote Nion grid to %s\n", fname_HI);
  printf("The total number of ionizing photons emitted at Snapshot %d was %.4ee50 photons/s\n", ReionSnap, grid_towrite->Nion_HI_Total); 

  return EXIT_SUCCESS;

}
