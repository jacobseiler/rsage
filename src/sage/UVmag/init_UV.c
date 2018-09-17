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

  stars_LUV = calloc(MAXBINS, sizeof(*(stars_Ngamma)));
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
