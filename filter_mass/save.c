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

int32_t save_arrays(int64_t *HaloID, float *ReionMod, SAGE_params params, int32_t NHalos_Ionized, int32_t filenr, int32_t ThisSnap, int32_t first_run)
{

  FILE *outfile;
  char outfile_name[MAXLEN];
  int32_t tmp_snap, tmp_NHalos, i;

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


  // If this is the first time the code is executed, create a pad for the snapshots below the current one. 
  if (first_run == 1)
  {
    tmp_NHalos = 0;
    for (i = 0; i < ThisSnap; ++i)
    {
      tmp_snap = i;

      fwrite(&tmp_snap, sizeof(int32_t), 1, outfile);
      fwrite(&tmp_NHalos, sizeof(int32_t), 1, outfile);
    }
  } 

  fwrite(&ThisSnap, sizeof(int32_t), 1, outfile); // Write out header info.
  fwrite(&NHalos_Ionized, sizeof(int32_t), 1, outfile);
  
  fwrite(HaloID, sizeof(int64_t), NHalos_Ionized, outfile); // Write out arrays.
  fwrite(ReionMod, sizeof(float), NHalos_Ionized, outfile);

  fclose(outfile);

  printf("All saved!\n");
 
  return EXIT_SUCCESS;

}

int32_t read_arrays(SAGE_params params, int32_t filenr, int32_t ThisSnap)
{

  FILE *infile;
  char infile_name[MAXLEN];
  int32_t snapread, NHalos_read, treenr, halonr, i, j;

  int64_t *HaloID;
  float *ReionMod;

  snprintf(infile_name, MAXLEN, "%s/reionization_modifiers/treefile_%03d", params->PhotoionDir, filenr);
  infile = fopen(infile_name, "rb");

  if (infile == NULL)
  {
    fprintf(stderr, "Could not open file %s\n", infile_name);
    return EXIT_FAILURE;
  }

  for (i = 0; i < ThisSnap + 1; ++i)
  {
    fread(&snapread, sizeof(int32_t), 1, infile);
    fread(&NHalos_read, sizeof(int32_t), 1, infile);

    printf("Snapshot %d has %d in region\n", snapread, NHalos_read);
    if (NHalos_read == 0)
    {
      continue;
    }

    HaloID = malloc(sizeof(int64_t) * NHalos_read);
    if (HaloID == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for HaloID read\n");
      return EXIT_FAILURE; 
    }  

    ReionMod = malloc(sizeof(float) * NHalos_read);   
    if (ReionMod == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for ReionMod read\n");
      return EXIT_FAILURE; 
    }  

    fread(HaloID, sizeof(int64_t), NHalos_read, infile);
    fread(ReionMod, sizeof(float), NHalos_read, infile);
      
    printf("After read\n"); 
    for (j = 0; j < NHalos_read; ++j)
    {
     
      treenr = (int32_t)(HaloID[j] >> 32);
      halonr = (int32_t)HaloID[j];

      printf("Unique HaloID %ld represents treenr %d halonr %d with reionization modifier %.4f\n", HaloID[j], treenr, halonr, ReionMod[j]);
 
    } 
    free(HaloID);
    free(ReionMod);
  }

  fclose(infile);
  return EXIT_SUCCESS;

}
