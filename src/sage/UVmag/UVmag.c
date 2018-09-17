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

#include "../core_allvars.h"

#define LUV_LOOKUPTABLE_MASS 1e6

int32_t init_UVlookup(void)
{

  // For calcUVmag == 1 we must explicitly track the stellar ages of a galaxy.
  // Then we determine the UV magnitude using the age of the stellar population.

  // Using STARBURST99 it was determined that the 1600A Luminosty depends on the mass
  // of stars formed and the time since the starburst. Furthermore, the UV luminosity linearly with the mass of stars formed.
  // That is, a starburst that forms 7.0e10Msun worth of stars will exhibit 10 times the UV luminosity as a starburst that forms 6.0e10Msun worth of stars.

  // So if we read in a table that contains the number of ionizing photons emitted from a starburst for a 6.0e10Msun episode, then we can scale our values to this lookup table
  // using log10 LUV(Msun, t) = (log10 M* - 6.0) + log10 LUV(6.0, t). 

#define MAXBINS 10000

  char fname[MAX_STRING_LEN];
  FILE *LUVtable;
  int32_t num_lines = 0;
  float t, lambda, LUV, norm_spec;

  snprintf(fname, MAX_STRING_LEN - 1, ROOT_DIR "/extra/LUV_table.txt");
  LUVtable = fopen(fname, "r");
  if (LUVtable == NULL)
  {
    fprintf(stderr, "Could not open file %s\n", fname);
    return EXIT_FAILURE;
  }

  stars_LUV = calloc(MAXBINS, sizeof(*(stars_LUV)));
  if (stars_LUV == NULL)
  {
    fprintf(stderr, "Could not allocate memory for the UV luminosity bins for the tracking of stellar populations.\n");
    return EXIT_FAILURE;
  }

  while (fscanf(LUVtable, "%f %f %f %f", &t, &lambda, &LUV, &norm_spec) == 4) 
  {

    stars_LUV[num_lines] = LUV;

    ++num_lines;
    if (num_lines == MAXBINS - 1)
    {
      fprintf(stderr, "Exceeding the maximum bins for the tracking of stellar populations.\n");
      return EXIT_FAILURE;
    }  
  }
  fclose(LUVtable);

  // Check that the Nion lookup table had enough datapoints to cover the time we're tracking the ages for. 
  if (t / 1.0e6 < STELLAR_TRACKING_TIME)
  {
    fprintf(stderr, "The final time specified in the UV luminosity lookup table is %.4f Myr. However we specified to track stellar ages over %d Myr.\n", t / 1.0e6, STELLAR_TRACKING_TIME);
    fprintf(stderr, "Either update the UV luminosity lookup table or reduce the value of STELLAR_TRACKING_TIME in `core_allvars.h`.\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

#undef MAXBINS
}

/*
Calculates the UV Luminosity (i.e., at 1600A) for a given galaxy at a specified snapshot.

**Important** See Units. 

Parameters
----------

g: struct GALAXY pointer.  See `core_allvars.h` for the struct architecture.
  Pointer to the galaxy that we are calculating UV luminosity. 

*LUV: Float Pointer.
  Pointer that will store the UV luminosity. 

Returns
----------

EXIT_SUCCESS or EXIT_FAILURE.
  If the UV luminosity is negative, EXIT_FAILURE is returned.  Otherwise EXIT_SUCCESS is returned.

Pointer Updates
----------

*LUV.

Units  
----------

The UV luminosty is returned in units of 1.0e50 erg s^-1 A^-1 
*/

int32_t calc_LUV(struct GALAXY *g, float *LUV)
{

  double t;
  int32_t i, lookup_idx;

  // We ran STARBURST99 for a single mass and scale our results.
  // First convert the mass to internal code units. 
  const double code_LUV_lookuptable_mass = LUV_LOOKUPTABLE_MASS*1.0e-10*Hubble_h;

  *LUV = 0.0;

  // Then go through the previous 100Myr worth of SF and sum up the UV luminosity.
  for (i = 0; i < StellarTracking_Len - 1; ++i)
  {
    if (g->Stellar_Stars[i] < 1e-10)
      continue;

    t = (i + 1) * TimeResolutionStellar; // (i + 1) because 0th entry will be at TimeResolutionSN.
    lookup_idx = (t / 0.1); // Find the index in the lookup table. 

    *LUV += exp10(log10(g->Stellar_Stars[i] / code_LUV_lookuptable_mass) + stars_LUV[lookup_idx] - 50.0);        
  }

  // The units of LUV are 1.0e50 erg s^-1 A^-1.  Hence a negative value is not allowed.
  if (*LUV < 0.0)
  {
    fprintf(stderr, "Got an LUV value of %.4e.  This MUST be a positive value.\nPrinting out information for every element used to calculate LUV.\n", *LUV);

    // Print out information for every element of the array so we can try identify the problem.       
    for (i = 0; i < StellarTracking_Len; ++i)
    {
      t = (i + 1) * TimeResolutionStellar; // (i + 1) because 0th entry will be at TimeResolutionSN.
      lookup_idx = (t / 0.1); // Find the index in the lookup table.

      double total = 0.0;
      total += exp10(log10(g->Stellar_Stars[i] / code_LUV_lookuptable_mass) + stars_LUV[lookup_idx] - 50.0);  
      printf("t %.4f\tlookup_idx %d\tg->Stellar_Stars[i] %.4e\tstars_LUV[lookup_idx] %.4e\tRunningTotal %.4e\n", t, lookup_idx, g->Stellar_Stars[i], stars_LUV[lookup_idx], total); 
    }
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
