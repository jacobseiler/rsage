#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>

#include "main.h"

#include "../grid-model/src/confObj.h"

int32_t check_ini_files(confObj_t *simParam, char *PhotoionDir, char *PhotoionName, int32_t fescPrescription, double alpha, double beta, double delta,
                        double quasar_baseline, double quasar_boosted, double N_dyntime, double MH_low, double fesc_low,
                        double MH_high, double fesc_high, int32_t HaloPartCut)
{

  // First check that the name of the photoionization file for SAGE and cifog matches.
  char SAGE_photHI_path[MAX_STRING_LEN]; 


  snprintf(SAGE_photHI_path, MAX_STRING_LEN - 1, "%s/%s", PhotoionDir, PhotoionName);

  // strncmp returns 0 if the two strings are equal.
  if (strncmp(SAGE_photHI_path, (*simParam)->out_photHI_file, MAX_STRING_LEN - 1) != 0)
  {
    printf("The SAGE and cifog ini files had different paths for the photoionization rate.\n");
    printf("Updating the cifog simParam variable to match the SAGE path.\n");

    (*simParam)->out_photHI_file = SAGE_photHI_path;
  }

  // Next adjust the nion_file variable using the parameters from the SAGE ini file.
  

  return EXIT_FAILURE;
  return EXIT_SUCCESS;
}
