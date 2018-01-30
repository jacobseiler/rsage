#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "core_allvars_grid.h"
#include "core_proto_grid.h"



int32_t init(void)
{
  int32_t i, status;
  

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, 42);	 // start-up seed 

  set_units(); // Change units to code units where necessary.
  srand((unsigned) time(NULL));

  read_snap_list(); // Read the snapshots required.

  for(i = 0; i < Snaplistlen; i++)
  {
    ZZ[i] = 1 / AA[i] - 1; // Redshift array.
    Age[i] = time_to_present(ZZ[i]); // Age array.
  }

  
  if (fescPrescription == 0)
  {
    fprintf(stderr, "\n\nUsing a constant escape fraction of %.4f\n", fesc); 
  }
  else if (fescPrescription == 1)
  {
    fprintf(stderr, "\n\nDeprecated! Use a new one!\n");
    return EXIT_FAILURE;
  } 
  else if (fescPrescription == 2)
  {
    fprintf(stderr,"\n\nUsing an fesc prescription that scales with halo mass.\n\n");
    determine_fesc_constants();
  }
  else if (fescPrescription == 3)
  {
    fprintf(stderr, "\n\nUsing an fesc prescription that scales with the fraction of ejected mass in the galaxy.\n");
    fprintf(stderr, "\nThis takes the form A*fej + B with A = %.4e and B = %.4e\n", alpha, beta); 
  }
  else if (fescPrescription == 4)
  {
    fprintf(stderr, "\n\nUsing an fesc prescription that depends upon quasar activity.\n");
    fprintf(stderr, "\nFor a galaxy that had a quasar event within %d dynamical times go the escape fraction will be %.2f.  Otherwise it will have a constant value of %.2f\n", N_dyntime, quasar_boosted, quasar_baseline);
  }
  else if (fescPrescription == 5)
  {
    fprintf(stderr, "\n\nUsing Anne's functional form for an escape fraction that decreases for increasing halo mass.\n");
    XASSERT(fesc_low > fesc_high, "Input file contain fesc_low = %.2f and fesc_high = %.2f. For this prescription (fescPrescription == 5), we require fesc_low > fesc_high\n", fesc_low, fesc_high);
  }
  else if (fescPrescription == 6)
  {
    fprintf(stderr, "\n\nUsing Anne's functional form for an escape fraction that increases for increasing halo mass.\n");
    XASSERT(fesc_low < fesc_high, "Input file contain fesc_low = %.2f and fesc_high = %.2f. For this prescription (fescPrescription == 6), we require fesc_low < fesc_high\n", fesc_low, fesc_high);
  }
  else
  {
    fprintf(stderr, "\n\nOnly escape fraction prescriptions 0 to 6 (exlucding 1) are permitted.\n");
    return EXIT_FAILURE;
  }

  Grid = malloc(sizeof(struct GRID_STRUCT));
  if (Grid == NULL)
  {
    fprintf(stderr, "Could not allocate memory for the high level grid structure\n");
    return EXIT_FAILURE;
  }

  status = init_grid(Grid);
  if (status == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

}


int32_t init_grid(struct GRID_STRUCT *grid)
{

  int32_t i;
  uint64_t cell_idx;

  if (grid == NULL)
  {
    fprintf(stderr, "init_grid was called with a GRID_STRUCT pointer that has not been initialized\n");
    return EXIT_FAILURE;
  } 

  
  grid->GridSize = GridSize;
  grid->NumCellsTotal = CUBE(GridSize);
  
  grid->NumGrids = NGrid; 

  grid->GridProperties = malloc(sizeof(struct GRID_PROPERTIES_STRUCT) * grid->NumGrids);
  if (grid->GridProperties == NULL)
  {
    fprintf(stderr, "Could not allocate memory for the grid properties struct\n");
    return EXIT_FAILURE;
  }

  for (i = 0; i < grid->NumGrids; ++i)
  {
    printf("Allocating grid %d\n", i);
    grid->GridProperties[i].SFR = malloc(sizeof(*(grid->GridProperties[i].SFR)) * grid->NumCellsTotal);
    if (grid->GridProperties[i].SFR == NULL)
    {
      fprintf(stderr, "Could not allocate memory for the grid SFR for grid number %d\n", i);
      return EXIT_FAILURE;
    }

    grid->GridProperties[i].StellarMass = malloc(sizeof(*(grid->GridProperties[i].StellarMass)) * grid->NumCellsTotal);
    if (grid->GridProperties[i].StellarMass == NULL)
    {
      fprintf(stderr, "Could not allocate memory for the grid StellarMass for grid number %d\n", i);
      return EXIT_FAILURE;
    }

    grid->GridProperties[i].Nion_HI = malloc(sizeof(*(grid->GridProperties[i].Nion_HI)) * grid->NumCellsTotal);
    if (grid->GridProperties[i].Nion_HI == NULL)
    {
      fprintf(stderr, "Could not allocate memory for the grid Nion_HI for grid number %d\n", i);
      return EXIT_FAILURE;
    }

    grid->GridProperties[i].Nion_HeI = malloc(sizeof(*(grid->GridProperties[i].Nion_HeI)) * grid->NumCellsTotal);
    if (grid->GridProperties[i].Nion_HeI == NULL)
    {
      fprintf(stderr, "Could not allocate memory for the grid Nion_HeI for grid number %d\n", i);
      return EXIT_FAILURE;
    }

    grid->GridProperties[i].Nion_HeII = malloc(sizeof(*(grid->GridProperties[i].Nion_HeII)) * grid->NumCellsTotal);
    if (grid->GridProperties[i].Nion_HeII == NULL)
    {
      fprintf(stderr, "Could not allocate memory for the grid Nion_HeII for grid number %d\n", i);
      return EXIT_FAILURE;
    }

    grid->GridProperties[i].GalCount = malloc(sizeof(*(grid->GridProperties[i].GalCount)) * grid->NumCellsTotal);
    if (grid->GridProperties[i].GalCount == NULL)
    {
      fprintf(stderr, "Could not allocate memory for the grid GalCount for grid number %d\n", i);
      return EXIT_FAILURE;
    }

    // All the memory has been allocated for the inner grid.
    // Now let's initialize the values.

    for (cell_idx = 0; cell_idx < Grid->NumCellsTotal; ++cell_idx)
    {
      grid->GridProperties[i].SFR[cell_idx] = 0.0;  
      grid->GridProperties[i].StellarMass[cell_idx] = 0.0;  
      grid->GridProperties[i].Nion_HI[cell_idx] = 0.0;  
      grid->GridProperties[i].Nion_HeI[cell_idx] = 0.0;  
      grid->GridProperties[i].Nion_HeII[cell_idx] = 0.0;  

      grid->GridProperties[i].GalCount[cell_idx] = 0; 
    } 

  } 
 
  printf("All grids have been initialized\n");
  return EXIT_SUCCESS; 
}

int32_t free_grid(void)
{

  int32_t i;

  for (i = 0; i < Grid->NumGrids; ++i)
  {
    free(Grid->GridProperties[i].GalCount);

    free(Grid->GridProperties[i].Nion_HeII);
    free(Grid->GridProperties[i].Nion_HeI);
    free(Grid->GridProperties[i].Nion_HI);
    free(Grid->GridProperties[i].StellarMass);
    free(Grid->GridProperties[i].SFR);
  }

  free(Grid->GridProperties);
  free(Grid);

  return EXIT_SUCCESS;

}


void set_units(void)
{

  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitTime_in_Megayears = UnitTime_in_s / SEC_PER_MEGAYEAR;
  G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);
  UnitDensity_in_cgs = UnitMass_in_g / pow(UnitLength_in_cm, 3);
  UnitPressure_in_cgs = UnitMass_in_g / UnitLength_in_cm / pow(UnitTime_in_s, 2);
  UnitCoolingRate_in_cgs = UnitPressure_in_cgs / UnitTime_in_s;
  UnitEnergy_in_cgs = UnitMass_in_g * pow(UnitLength_in_cm, 2) / pow(UnitTime_in_s, 2);

  // convert some physical input parameters to internal units 
  Hubble = HUBBLE * UnitTime_in_s;

  // compute a few quantitites 
  RhoCrit = 3 * Hubble * Hubble / (8 * M_PI * G);

}

void read_snap_list(void)
{
  FILE *fd;
  char fname[MAXLEN];

  snprintf(fname, MAXLEN, "%s", FileWithSnapList);

  if(!(fd = fopen(fname, "r")))
  {
    printf("can't read output list in file '%s'\n", fname);
    exit(0);
  }

  Snaplistlen = 0;
  do
  {
    if(fscanf(fd, " %lg ", &AA[Snaplistlen]) == 1)
      Snaplistlen++;
    else
      break;
  }
  while(Snaplistlen < MAXSNAPS);

  fclose(fd);

#ifdef MPI
  if(ThisTask == 0)
#endif
    printf("found %d defined times in snaplist\n", Snaplistlen);
}

double time_to_present(double z)
{
#define WORKSIZE 1000
  gsl_function F;
  gsl_integration_workspace *workspace;
  double time, result, abserr;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &integrand_time_to_present;

  gsl_integration_qag(&F, 1.0 / (z + 1), 1.0, 1.0 / Hubble,
    1.0e-8, WORKSIZE, GSL_INTEG_GAUSS21, workspace, &result, &abserr);

  time = 1 / Hubble * result;

  gsl_integration_workspace_free(workspace);

  // return time to present as a function of redshift 
  return time;
}



double integrand_time_to_present(double a, void *param)
{
  return 1 / sqrt(Omega / a + (1 - Omega - OmegaLambda) + OmegaLambda * a * a);
}

void determine_fesc_constants(void)
{

  double A, B, log_A;

  if (fescPrescription == 2)
  {
    log_A = (log10(fesc_high) - (log10(fesc_low)*log10(MH_high)/log10(MH_low))) * pow(1 - (log10(MH_high) / log10(MH_low)), -1);
    B = (log10(fesc_low) - log_A) / log10(MH_low);
    A = pow(10, log_A);

    alpha = A;
    beta = B;

    printf("Fixing the points (%.4e, %.2f) and (%.4e, %.2f)\n", MH_low, fesc_low, MH_high, fesc_high);
    printf("This gives a power law with constants A = %.4e, B = %.4e\n", alpha, beta);
  }

  
}
/*
void estimate_grid_memory(void)
{

  if (Verbose == 1)
  {
    printf("The size of the Grid Struct (that contains all the pointers) is %lu bytes\nEach grid cell has a struct, in addition to data that the pointers point to.\n", sizeof(struct GRID));
    printf("The total memory taken up by the Grid will be the sizeof a grid cell (%lu bytes) times the number of grid cells (%d) times the number of output snapshots (%d) times 2 because we need memory for both the pointers for the grid and the actual numbers themselves.\n", sizeof(struct GRID), CUBE(GridSize), NGrid);
  }

  printf("Approximate total memory for the grid ===== %lu mb\n", sizeof(struct GRID)*CUBE(GridSize)*NGrid*2/1024/1024);

}

void estimate_gal_memory(int NtotGals)
{
  if (Verbose == 1)
  {
    printf("The size of the Galaxy struct is %lu bytes and the size of the GalaxyGrid struct (that contains all the pointers for the SnapShot history) is %lu\n", sizeof(struct GALAXY_INPUT), sizeof(struct GALAXY_GRID));
    printf("The total memory taken up by this file of galaxies will be the sizeof the galaxy input (%lu bytes) times the number of galaxies (%d) plus the sizeof the galaxygrid (%lu bytes) times the number of snapshots in the simulation (%d) times the number of galaxies (%d)\n", sizeof(struct GALAXY_INPUT), NtotGals, sizeof(struct GALAXY_GRID), LastSnapShotNr, NtotGals); 
  }

  printf("Approximate memory occupied for this file of galaxies ===== %lu mb\n", (sizeof(struct GALAXY_INPUT)*NtotGals + sizeof(struct GALAXY_GRID)*NtotGals*LastSnapShotNr)/1024/1024);

}
*/
