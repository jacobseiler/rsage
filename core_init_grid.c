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



void init(void)
{
  int i;
  double evolve_time;

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

  GalaxyCount = malloc(sizeof(int)*NGrid);
  MajorMergerCount = malloc(sizeof(int)*NGrid); 
  MinorMergerCount = malloc(sizeof(int)*NGrid); 
  DiskCount = malloc(sizeof(int)*NGrid); 
  ICSCount = malloc(sizeof(int)*NGrid); 

  for(i = 0; i < NGrid; ++i)
  { 
    if (Verbose == 1)
    {
      if (i == NGrid - 1) {
        evolve_time = (Age[ListOutputGrid[i]-1] - Age[ListOutputGrid[i]])*UnitTime_in_Megayears / Hubble_h; 
//        printf("Evolve time for Snapshots %d and %d is %.3e Myr.\n", ListOutputGrid[i], ListOutputGrid[i]-1, evolve_time);
      }
      else {  
        evolve_time = (Age[ListOutputGrid[i+1]] - Age[ListOutputGrid[i]])*UnitTime_in_Megayears / Hubble_h;
//        printf("Evolve time for Snapshots %d and %d is %.3e Myr.\n", ListOutputGrid[i], ListOutputGrid[i+1], evolve_time);
      }
    }
    GalaxyCount[i] = 0;
    MajorMergerCount[i] = 0;
    MinorMergerCount[i] = 0;
    DiskCount[i] = 0;
    ICSCount[i] = 0;

  }

  if (fescPrescription == 2)
  {
    beta = (delta - kappa)/(MH_max - MH_min); // Index of the power law.	
    alpha = pow(10, kappa - beta * MH_min);  // Constant out the front of the power law.
    printf("kappa = %.4f \t MH_min = %.4f \t pow(10,MH_min = %.4e \t -beta = %.4f\n", kappa, MH_min, pow(10, MH_min), beta);
    printf("\n\nUsing an fesc prescription that scales with halo mass.\n");
    printf("The smallest halo mass is %.4e and is fixed to have fesc = %.4e, and the largest halo mass is %.4e and is fixed to have fesc = %.4e\n", pow(10,MH_min),  pow(10, kappa), pow(10, MH_max), pow(10, delta)); 
    printf("This results in a power law of the form fesc = A*MH^B with A = %.4e and B = %.4e\n\n", alpha, beta);
  }
  else if (fescPrescription == 3)
  {
    beta = pow(10, delta) - pow(10, kappa); // Slope of the linear relation. 
    alpha = pow(10, kappa);  // Intercept.
    printf("\n\nUsing an fesc prescription that scales with the fraction of ejected mass in the galaxy.\n"); 
    printf("We are setting an Ejected Fraction of 0 to have escape fraction %.4e and an Ejected Fraction of 1 to have escape fraction %.4e\n", pow(10, kappa), pow(10, delta)); 
    printf("This results in a linear relationship of the form fesc = B*Ejected + A with A = %.4e and B = %.4e\n\n", alpha, beta);
  }
}


void init_grid(int GridNr)
{
  int i;

  if (GridNr == 0) // Only need to malloc the grid once.
    Grid = mymalloc(sizeof(struct GRID)*CUBE(GridSize));

  //estimate_grid_memory();
  
  for(i = 0; i < CUBE(GridSize); ++i)
  {
    Grid[i].Sfr = 0.0;
    Grid[i].StellarMass = 0.0;

    Grid[i].Density = 0.0;
    Grid[i].Nion_HI = 0.0;
    Grid[i].Nion_HeI = 0.0;
    Grid[i].Nion_HeII = 0.0;

    Grid[i].Diffuse = 0;
    Grid[i].Count = 0;

    Grid[i].HaloMass = 0.0;
    Grid[i].HaloCount = 0;
 
  }
 
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
  char fname[1000];

  sprintf(fname, "%s", FileWithSnapList);

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
