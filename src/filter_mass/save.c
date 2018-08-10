#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>

#include "filter_mass.h"

#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

#define MAXLEN 1024

int32_t save_arrays(int64_t *HaloID, float *ReionMod, struct SAGE_PARAMETERS *params, int32_t NHalos_Ionized, int32_t filenr, int32_t ThisSnap, int32_t first_update_flag)
{

  FILE *outfile;
  char outfile_name[MAXLEN];
  int32_t tmp_snap, tmp_NHalos, halo_idx, nwritten;
  
  snprintf(outfile_name, MAXLEN - 1, "%s/reionization_modifiers/%s_treefile_%03d", params->PhotoionDir, params->FileNameGalaxies, filenr);

  // If this is the first time the code is executed, need to create a new file. Otherwise append to the end.
  if (first_update_flag == 0)
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
  if (first_update_flag == 1)
  {
    tmp_NHalos = 0;
    for (halo_idx = 0; halo_idx < ThisSnap; ++halo_idx)
    {
      tmp_snap = halo_idx;

      nwritten = fwrite(&tmp_snap, sizeof(int32_t), 1, outfile);
      if (nwritten != 1)
      {
        fprintf(stderr, "Only wrote %d elements for snapshot %d when trying to create the snapshot pad.  Needed to write 1 elements.\n", nwritten, halo_idx);
      }
 
      nwritten = fwrite(&tmp_NHalos, sizeof(int32_t), 1, outfile);
      if (nwritten != 1)
      {
        fprintf(stderr, "Only wrote %d elements for snapshot %d when trying to create the NHalo pad.  Needed to write 1 elements.\n", nwritten, halo_idx);
      } 
    }
  } 

  printf("Writing the snapshot number %d\n", ThisSnap);
  nwritten = fwrite(&ThisSnap, sizeof(int32_t), 1, outfile); // Write out header info.
  if (nwritten != 1)
  {
    fprintf(stderr, "Only wrote %d elements when writing out the current snapshot number.  Needed to write 1 elements.\n", nwritten);
  } 

  printf("Writing the number of ionized halos %d\n", NHalos_Ionized);
  nwritten = fwrite(&NHalos_Ionized, sizeof(int32_t), 1, outfile); 
  if (nwritten != 1)
  {
    fprintf(stderr, "Only wrote %d elements when writing out the number of ionized halos.  Needed to write 1 elements.\n", nwritten);
  } 
 
  nwritten = fwrite(HaloID, sizeof(int64_t), NHalos_Ionized, outfile); // Write out arrays.
  if (nwritten != NHalos_Ionized)
  {
    fprintf(stderr, "Only wrote %d elements when writing out the IDs of the ionized halos.  Needed to write %d elements.\n", nwritten, NHalos_Ionized);
  } 

  nwritten = fwrite(ReionMod, sizeof(float), NHalos_Ionized, outfile);
  if (nwritten != NHalos_Ionized)
  {
    fprintf(stderr, "Only wrote %d elements when writing out the reionization modifier of the ionized halos.  Needed to write %d elements.\n", nwritten, NHalos_Ionized);
  } 

  fclose(outfile);
 
  return EXIT_SUCCESS;

}

int32_t read_arrays(struct SAGE_PARAMETERS *params, int32_t filenr, int32_t ThisSnap)
{

  FILE *infile;
  char infile_name[MAXLEN];
  int32_t snapread, NHalos_read, treenr, halonr, i, j;

  int64_t *HaloID;
  float *ReionMod;

  snprintf(infile_name, MAXLEN, "%s/reionization_modifiers/%s_treefile_%03d", params->PhotoionDir, params->FileNameGalaxies, filenr);
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
