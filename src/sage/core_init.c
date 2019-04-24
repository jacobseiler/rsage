#define _GNU_SOURCE
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
#include "UVmag/UVmag.h"

// Local Proto-Types //

int32_t init_delayedSN(void);
int32_t init_nionlookup(void);
int32_t init_metalcooling(void);
void read_snap_list(void);
void set_units(void);
int32_t init_fesc(void);
int32_t determine_fescMH_constants(void);

// External Functions //

void sage_init(void)
{
  int32_t i, status;

  count_gal = 0;  
  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, 42);	 // start-up seed 

  set_units();
  srand((unsigned) time(NULL));

  read_snap_list();
  for(i = 0; i < MAXSNAPS; i++)
  {
    ZZ[i] = 1 / AA[i] - 1;
    Age[i] = time_to_present(ZZ[i]);
  }

  a0 = 1.0 / (1.0 + Reionization_z0);
  ar = 1.0 / (1.0 + Reionization_zr);

  read_cooling_functions();

  status = init_fesc();
  if (status != EXIT_SUCCESS)
  {
    ABORT(0);
  }

  count_onehalo = 0;  

  zeromass_count = 0;
  suppression_count = 0;
  previous_tree = 0;

  smallest_mass = 100000;
  lowmass_halo = 0;

  outside_box = 0;
  inside_box = 0;

  count_Mvir = 0;
  count_Len = 0;

  if (IMF == 0)
  {
    // Salpeter IMF //
    IMF_norm =0.1706;
    IMF_slope = -2.35;
    Eta_SNII = 7.432e-3; //Number fraction of stars that will end their lives as type II supernovae.
    m_SNII = 0.144; // Mass fraction of stars that will end their lives as type II supernovae.

  } else if (IMF == 1)
  {
    // Chabrier IMF //
    IMF_norm = 0.23638; 
    IMF_slope = -2.3;
    Eta_SNII = 1.1893e-2; //Number fraction of stars that will end their lives as type II supernovae.
    m_SNII = 0.23638; // Mass fraction of stars that will end their lives as type II supernovae.

  }
 
  // Tiamat Parameters //

  V_energy = 70.0;
  alpha_energy = 0.5;
  beta_energy = 10.0;

  V_mass = 70.0; 
  alpha_mass = 6.0; 
  beta_mass = 10.0; 

  epsilon_mass_max = 30.0;

  if (IRA == 0)
  {
    status = init_delayedSN();
    if (status != EXIT_SUCCESS)
    {
      ABORT(EXIT_FAILURE);
    }  
  }

  if (PhotonPrescription == 1)
  {
    status = init_nionlookup();
    if (status != EXIT_SUCCESS)
    {
      ABORT(EXIT_FAILURE);
    }
  }

  if (calcUVmag == 1)
  {
    status = init_UVlookup();
    if (status != EXIT_SUCCESS)
    {
      ABORT(EXIT_FAILURE);
    }
  }

  if ((PhotonPrescription == 1) || (calcUVmag == 1))
  {
    StellarTracking_Len = STELLAR_TRACKING_TIME / TimeResolutionStellar;
  }
 
  mergedgal_mallocs = 0;
  gal_mallocs = 0 ;

  mergedgal_frees = 0;
  gal_frees = 0;

  Ngamma_HI_Total = 0.0;
 
}

// Local Functions //

int32_t init_delayedSN(void)
{ 

  // We keep the past 50 Myrs of star formation stored for each galaxy.
  // Based on the TimeResolution specified need to find out how many elements this corresponds to. 

  if(TimeResolutionSN > 50)
  {
    fprintf(stderr, "The selected time resolution for SN feedback (TimeResolutionSN) is set too high (%d Myr).  Using TimeResolutionSN > 50Myr is the same as using the instantaneous recycling approximation; set 'IRA' to 1 instead!\n", TimeResolutionSN); 
    ABORT(EXIT_FAILURE);  
  } else if(TimeResolutionSN > 35)
  {
    fprintf(stderr, "Your selected time resolution for SN feedback (TimeResolutionSN) is quite high (%d Myr).  Beyond 50Myr the instantaneous recycling approximation is valid hence with your value it would likely be correct to set 'IRA' to 1.\n", TimeResolutionSN);
  } else
  {
    Time_SFH = 0;
    SN_Array_Len = 0;
    while(Time_SFH < 50)
    {
      Time_SFH += TimeResolutionSN;
      ++SN_Array_Len;
    }

  }

  // For the delayed SN scheme the number of SN for each time step depends on the IMF.
  // As we will need to calculate this for every substep we create look-up tables for efficiency. 

  // First we need to know given a timestep, what is the minimum mass star that will go supernova? 

  double a = 0.7473; // Fits from Portinari et al. (1998). 
  double b = -2.6979;
  double c = -4.7659;
  double d = 0.5934;

  coreburning_tbins_low = 2; // Any time below 2Myr will result in stars > 120Myr to explode. This is beyond the limits of the IMF.
  coreburning_tbins_high = 45; // Any time above 45Myr will result in stars > 8Myr to explode. Can just return 8.0 in this case. 
  coreburning_tbins_delta = 0.00001;

  N_tbins = (coreburning_tbins_high - coreburning_tbins_low) / (coreburning_tbins_delta);
  int32_t bin_idx;

  coreburning_times = malloc(sizeof(*(coreburning_times)) * N_tbins);
  if (coreburning_times == NULL)
  {
    fprintf(stderr, "Could not allocate memory for the grid to contains the coreburning times.\n");
    return EXIT_FAILURE;
  }

  for (bin_idx = 0; bin_idx < N_tbins; ++bin_idx)
  {
    double t = coreburning_tbins_low + ((double)bin_idx * coreburning_tbins_delta); 
    coreburning_times[bin_idx] = exp10(a/log10(t) + b * exp(c/log10(t)) + d); 
  }

  // Now need to know given a mass range (defined by the timestep) what is the number and mass fraction of SN that go supernova. 

  m_IMFbins_low = 8.0; // The IMF range is from 8.0 to 120.0 Msun.
  m_IMFbins_high = 120.0;
  m_IMFbins_delta = 0.00001;

  N_massbins = (m_IMFbins_high - m_IMFbins_low) / (m_IMFbins_delta);
  
  IMF_massgrid_eta = malloc(sizeof(*(IMF_massgrid_eta)) * N_massbins); 
  if (IMF_massgrid_eta == NULL)
  {
    fprintf(stderr, "Could not allocate memory for the grid to contain the IMF eta masses for delayed supernova.\n");
    return EXIT_FAILURE;
  }

  IMF_massgrid_m = malloc(sizeof(*(IMF_massgrid_m)) * N_massbins); 
  if (IMF_massgrid_m == NULL)
  {
    fprintf(stderr, "Could not allocate memory for the grid to contain the IMF m masses for delayed supernova.\n");
    return EXIT_FAILURE;
  }

  for (bin_idx = 0; bin_idx < N_massbins; ++bin_idx)
  {
    IMF_massgrid_eta[bin_idx] = pow(m_IMFbins_low + ((double)bin_idx * m_IMFbins_delta), IMF_slope + 1.0);
    IMF_massgrid_m[bin_idx] = pow(m_IMFbins_low + ((double)bin_idx * m_IMFbins_delta), IMF_slope + 2.0);
  }

  return EXIT_SUCCESS;

}

int32_t init_nionlookup(void)
{

  // For PhotonPrescription == 1 we wish to explicitly track the stellar ages of a galaxy.
  // Then we determine the number of ionizing photons using the age of the stellar population.

  // Using STARBURST99 it was determined that the number of ionizing photons emitted from an instantaneous starburst depends on the 
  // mass of stars formed and the time since the starburst.  Furthermore, the number of ionizing photons scales linearly with the mass of stars formed.
  // That is, a starburst that forms 7.0e10Msun worth of stars will emit 10x as many photons as a starburst that forms 6.0e10Msun worth of stars.

  // So if we read in a table that contains the number of ionizing photons emitted from a starburst for a 6.0e10Msun episode, then we can scale our values to this lookup table
  // using log10 Ngamma(Msun, t) = (log10 M* - 6.0) + log10 Ngamma(6.0, t). 

#define MAXBINS 10000

  char buf[MAX_STRING_LEN], fname[MAX_STRING_LEN];
  FILE *niontable;
  int32_t i = 0, num_lines = 0;
  float t, HI, HI_L, HeI, HeI_L, HeII, HeII_L, L;

  snprintf(fname, MAX_STRING_LEN - 1, ROOT_DIR "/extra/nion_table.txt");
  niontable = fopen(fname, "r");
  if (niontable == NULL)
  {
    fprintf(stderr, "Could not open file %s\n", fname);
    return EXIT_FAILURE;
  }

  stars_tbins = calloc(MAXBINS, sizeof(*(stars_tbins)));
  if (stars_tbins == NULL)
  {
    fprintf(stderr, "Could not allocate memory for the time bins for the tracking of stellar populations.\n");
    return EXIT_FAILURE;
  }

  stars_Ngamma = calloc(MAXBINS, sizeof(*(stars_Ngamma)));
  if (stars_Ngamma == NULL)
  {
    fprintf(stderr, "Could not allocate memory for the Ngamma HI bins for the tracking of stellar populations.\n");
    return EXIT_FAILURE;
  }

  // The first 8 lines are all header information. So move past them.
  while (i < 8)
  {
    fgets(buf, MAX_STRING_LEN, niontable);
    ++i;
  }

  while (fscanf(niontable, "%f %f %f %f %f %f %f %f", &t, &HI, &HI_L, &HeI, &HeI_L, &HeII, &HeII_L, &L) == 8) 
  {

    stars_tbins[num_lines] = t;
    stars_Ngamma[num_lines] = HI;

    ++num_lines;
    if (num_lines == MAXBINS - 1)
    {
      fprintf(stderr, "Exceeding the maximum bins for the tracking of stellar populations.\n");
      return EXIT_FAILURE;
    }  
  }
  fclose(niontable);

  // Check that the Nion lookup table had enough datapoints to cover the time we're tracking the ages for. 
  if (stars_tbins[num_lines - 1] / 1.0e6 < STELLAR_TRACKING_TIME)
  {
    fprintf(stderr, "The final time specified in the Nion lookup table is %.4f Myr. However we specified to track stellar ages over %d Myr.\n", stars_tbins[num_lines - 1] / 1.0e6, STELLAR_TRACKING_TIME);
    fprintf(stderr, "Either update the Nion lookup table or reduce the value of STELLAR_TRACKING_TIME in `core_allvars.h`.\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

#undef MAXBINS
}

void set_units(void)
{
  
  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;  
  UnitTime_in_Megayears = UnitTime_in_s / SEC_PER_MEGAYEAR;
  UnitDensity_in_cgs = UnitMass_in_g / pow(UnitLength_in_cm, 3);
  UnitPressure_in_cgs = UnitMass_in_g / UnitLength_in_cm / pow(UnitTime_in_s, 2);
  UnitCoolingRate_in_cgs = UnitPressure_in_cgs / UnitTime_in_s;
  UnitEnergy_in_cgs = UnitMass_in_g * pow(UnitLength_in_cm, 2) / pow(UnitTime_in_s, 2);

  EnergySNcode = EnergySN / UnitEnergy_in_cgs * Hubble_h;

  sage_G = GRAVITY / pow(UnitLength_in_cm, 3) * UnitMass_in_g * pow(UnitTime_in_s, 2);
  
  // convert some physical input parameters to internal units
  sage_Hubble = HUBBLE * UnitTime_in_s;

  // compute a few quantitites 
  RhoCrit = 3 * sage_Hubble * sage_Hubble / (8 * M_PI * sage_G);
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
}

double time_to_present(double z)
{
#define WORKSIZE 1000
  gsl_function F;
  gsl_integration_workspace *workspace;
  double time, result, abserr;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &integrand_time_to_present;

  gsl_integration_qag(&F, 1.0 / (z + 1), 1.0, 1.0 / sage_Hubble,
    1.0e-8, WORKSIZE, GSL_INTEG_GAUSS21, workspace, &result, &abserr);

  time = 1 / sage_Hubble * result;

  gsl_integration_workspace_free(workspace);

  // return time to present as a function of redshift 
  return time;
}

double integrand_time_to_present(double a, void *param)
{
  (void)(param); // Avoids triggering unused-parameter warning.
  return 1 / sqrt(Omega / a + (1 - Omega - OmegaLambda) + OmegaLambda * a * a);
}


int32_t init_fesc(void)
{
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
        printf("\n\nUsing Anne's functional form for an escape fraction that increases for descreasing halo mass.\n");
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
        printf("This takes the function form of a logistic function, fesc = delta / (1.0 + exp(-alpha*(log10(SFR)-beta))), with range [0.0, delta].  `alpha` controls the steepness of the curve and beta defines the value of log10(SFR) gives fesc = delta/2.\n");
        printf("alpha = %.4f\tbeta = %.4f\tdelta = %.4f\n", alpha, beta, delta);
      }
      break;

    default:
      printf("\n\nOnly escape fraction prescriptions 0, 1, 2, 3, 4 or 5 are permitted.\n");
      return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

/*
For some escape fraction prescriptions, the functional form is a power law with the constants calculated using two fixed points.
This function determines these constants based on the fixed points specified.
Refer to the .ini file or `determine_fesc()` function for full details on each `fescPrescription`.

Parameters
----------
None.  All variables used/adjusted are global.

Returns
----------
None. All variables adjusted are global
.
Pointer Updates
----------
None.

Units  
----------
The Halo Masses used to specify the fixed points (MH_low and MH_high) are in Msun.
*/

int32_t determine_fescMH_constants(void)
{ 
  
  double A, B, log_A;
 
  switch(fescPrescription)
  {
    case 3:
      log_A = (log10(fesc_high) - (log10(fesc_low)*log10(MH_high)/log10(MH_low))) * pow(1 - (log10(MH_high) / log10(MH_low)), -1);
      B = (log10(fesc_low) - log_A) / log10(MH_low);
      A = pow(10, log_A);
      
      fescMH_alpha = A;
      fescMH_beta = B;

      printf("Fixing the points (%.4e, %.2f) and (%.4e, %.2f)\n", MH_low, fesc_low, MH_high, fesc_high);
      printf("This gives a power law with constants A = %.4e, B = %.4e\n", fescMH_alpha, fescMH_beta);

      break;
      
  }

  return EXIT_SUCCESS;

}
