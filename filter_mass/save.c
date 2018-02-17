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

#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

#define MAXLEN 1024

int32_t bubble_sort(int64_t *HaloID, float *ReionMod, int32_t NHalos_Ionized)
{

  int32_t swapped, i;

  printf("BRBEURBQ\n");
  do
  {
    swapped = 0;

    for (i = 0; i < NHalos_Ionized-1; ++i)
    {
      if ((HaloID)[i] > (HaloID)[i+1])
      {
        int64_t temp = (HaloID)[i];
        (HaloID)[i] = (HaloID)[i+1];
        (HaloID)[i+1] = temp;

        float tmp = (ReionMod)[i];
        (ReionMod)[i] = (ReionMod)[i+1];
        (ReionMod)[i+1] = tmp;

        swapped = 1;
      }          
    }
  } while (swapped == 1);

  return EXIT_SUCCESS;;
}

int32_t save_arrays(int64_t *HaloID, float *ReionMod, SAGE_params params, int32_t NHalos_Ionized, int32_t filenr, int32_t first_run)
{

  FILE *outfile;
  char outfile_name[MAXLEN];
  int32_t status;

  printf("Saving the arrays now.\n");

  snprintf(outfile_name, MAXLEN, "%s/reionization_modifiers/treefile_%03d", params->PhotoionDir, filenr);

  // If this is the first time the code is executed, need to create a new file. Otherwise append to the end.
  if (first_run == 0)
  {
    outfile = fopen(outfile_name, "ab");
  }
  else
  {
    outfile = fopen(outfile_name, "wb");
  }

  if (outfile == NULL)
  {
    fprintf(stderr, "Could not open file %s\n", outfile_name);
    return EXIT_FAILURE;
  }

  printf("Saving to %s\n", outfile_name);
  
  status = bubble_sort(HaloID, ReionMod, NHalos_Ionized);
 
  return EXIT_SUCCESS;

}
