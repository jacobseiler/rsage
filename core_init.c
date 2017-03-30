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

#include "core_allvars.h"
#include "core_proto.h"

void init(void)
{
  int i;

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, 42);	 // start-up seed 

  set_units();
  srand((unsigned) time(NULL));

  read_snap_list();

  for(i = 0; i < Snaplistlen; i++)
  {
    ZZ[i] = 1 / AA[i] - 1;
    Age[i] = time_to_present(ZZ[i]);

  }

  a0 = 1.0 / (1.0 + Reionization_z0);
  ar = 1.0 / (1.0 + Reionization_zr);

  read_cooling_functions();

  count_onehalo = 0;  

  zeromass_count = 0;
  suppression_count = 0;
  previous_tree = 0;

  smallest_mass = 100000;
  lowmass_halo = 0;

  count = 0;

  IMF_norm = 0.380417; // IMF Norm calculated from Chabrier IMF.
  //IMF_norm = 0.1706; // IM Normalization calculated from Salpeter IMF.

  Eta_SNII = IMF_norm/0.1706 * 7.432e-3; //Number fraction of stars that will end their lives as type II supernovae.
  m_SNII = IMF_norm/0.1706 * 0.144; // Mass fraction of stars that will end their lives as type II supernovae.
  // This number was obtained by using the Mutch et al. (2016) value and then rescaling using the IMF normalization for a Chabrier IMF.
  // They used a Salpeter IMF and gotten (Phi_norm, Eta_SNII) = (0.1706, 7.432e-3).
  // As I calculated Phi_norm = 0.380417 this then corresponds to an Eta_SNII that is ~2.299 times larger than theirs.

  alpha_energy = 0.18;
  beta_energy = 3.2;
  V_energy = 70.0;
  alpha_mass = 4.0;
  beta_mass = 3.2;
  V_mass = 70.0;
  epsilon_mass_max = 10.0;

}


void init_grid()
{

  int i,j;
  FILE *load_fd;
  char buf[1000];

  PhotoGrid = mymalloc(sizeof(struct PHOTO_GRID)*CUBE(GridSize)); 
  
  printf("Reading the photoionization and redshift of reionization grids.\n");
  for(i = 0; i < CUBE(GridSize); ++i)
  {
    if (NULL == (PhotoGrid[i].PhotoRate = malloc(sizeof(double)*MAXSNAPS)))  
    {   
      fprintf(stderr, "Out of memoery allocating %ld bytes, could not allocate PhotoGrid[i].PhotoRate.", sizeof(double)*MAXSNAPS);
      exit(EXIT_FAILURE);
    }
    
  }

  printf("Maxsnaps = %d\n", MAXSNAPS); 
  for(j = 0; j < MAXSNAPS; ++j)
  {

    if (j < 10)
      sprintf(buf, "%s/%s_0%d", PhotoionDir, PhotoionName, j);
    else
      sprintf(buf, "%s/%s_%d", PhotoionDir, PhotoionName, j);

    if(!(load_fd = fopen(buf, "r")))
    {
      printf("can't open file `%s'\n", buf);
      ABORT(0);
    }
   
    for(i = 0; i < CUBE(GridSize); ++i)
    {
      fread(&PhotoGrid[i].PhotoRate[j], sizeof(double), 1, load_fd);
    } 
 
    fclose(load_fd);


  }

  sprintf(buf, "%s/%s", PhotoionDir, ReionRedshiftName);

  if(!(load_fd = fopen(buf, "r")))
  {
    printf("can't open file `%s'\n", buf);
    ABORT(0);
  }

  for(i = 0; i < CUBE(GridSize); ++i)
  {   
    fread(&PhotoGrid[i].ReionRedshift, sizeof(double), 1, load_fd);
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

  EnergySNcode = EnergySN / UnitEnergy_in_cgs * Hubble_h;
  EtaSNcode = EtaSN * (UnitMass_in_g / SOLAR_MASS) / Hubble_h;

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
    ABORT(0);
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



