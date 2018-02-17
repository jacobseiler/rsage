#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>

#include "tree_io.h"

FILE *forest_file = NULL;
// keep a static file handle to remove the need to do constant seeking

int32_t load_tree_table(int32_t filenr, SAGE_params params, int32_t *Ntrees, int32_t *totNHalos, int32_t **TreeNHalos)
{
  char forestname[MAXLEN]; 

  snprintf(forestname, MAXLEN, "%s/%s_%03d.dat", params->TreeDir, params->TreeName, filenr);

  fprintf(stderr, "Reading file %s\n", forestname);
  if(!(forest_file= fopen(forestname, "rb")))
  {
    printf("can't open file `%s'\n", forestname);
    return EXIT_FAILURE; 
  }

  fread(Ntrees, 1, sizeof(int32_t), forest_file);
  fread(totNHalos, 1, sizeof(int32_t), forest_file);

  *TreeNHalos = malloc(sizeof(int) * (*Ntrees));
  if (*TreeNHalos == NULL)
  {
    fprintf(stderr, "Could not allocate memory for TreeNHalos.\n");
    return EXIT_FAILURE;
  }
  
  fread(*TreeNHalos, (*Ntrees), sizeof(int), forest_file); 

  fprintf(stderr, "Read the table\n");
    
  return EXIT_SUCCESS;

}

int32_t load_halos(int32_t treenr, int32_t NHalos_ThisTree, halo_t *Halos)
{

  *Halos = malloc(sizeof(struct HALO_STRUCT) * NHalos_ThisTree);
  if (*Halos == NULL)
  {
    fprintf(stderr, "Could not allocate memory for Halos in tree %d\n", treenr);
    return EXIT_FAILURE;
  }

  fread(*Halos, NHalos_ThisTree, sizeof(struct HALO_STRUCT), forest_file); 
  return EXIT_SUCCESS;

} 


int32_t allocate_array_memory(int32_t totNHalos, int64_t **HaloID, float **ReionMod)
{

  *HaloID = malloc(sizeof(*(*HaloID)) * totNHalos);
  if (*HaloID == NULL)
  {
    fprintf(stderr, "Could not allocate memory for HaloID.\n");
    return EXIT_FAILURE;
  }

  *ReionMod = malloc(sizeof(*(*ReionMod)) * totNHalos);
  if (*ReionMod == NULL)
  {
    fprintf(stderr, "Could not allocate memory for ReionMod.\n");
    return EXIT_FAILURE;
  } 

  return EXIT_SUCCESS;

}

int32_t free_memory(int32_t **TreeNHalos, int64_t **HaloID, float **ReionMod)
{

  free(*TreeNHalos);
  free(*HaloID);
  free(*ReionMod);

  if (forest_file)
  {
    fclose(forest_file);
    forest_file = NULL;
  }
    

  return EXIT_SUCCESS;

}

int32_t trim_arrays(int64_t *HaloID, float *ReionMod, int32_t NHalos_Ionized)
{

  int64_t *HaloID_tmp;
  float *ReionMod_tmp;
  int32_t i;

  printf("Trimming the arrays to the number of those in the ionized regions.\n");

  HaloID_tmp = malloc(sizeof(*(HaloID_tmp)) * NHalos_Ionized);
  if (HaloID_tmp == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for the tempory HaloID array (for the trimming\n");
    return EXIT_FAILURE;
  }

  ReionMod_tmp = malloc(sizeof(*(ReionMod_tmp)) * NHalos_Ionized);
  if (ReionMod_tmp == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for the tempory ReionMod array (for the trimming\n");
    return EXIT_FAILURE;
  }

  for (i = 0; i < NHalos_Ionized; ++i)
  {       
    HaloID_tmp[i] = (HaloID)[i];  
    ReionMod_tmp[i] = (ReionMod)[i];
  } 

  printf("Freeing the originals\n"); 
  free(HaloID);
  free(ReionMod);

  HaloID = HaloID_tmp;
  ReionMod = ReionMod_tmp;

  printf("Successfully reassigned, freeing the temp.\n");
  free(HaloID_tmp);
  free(ReionMod_tmp);

  return EXIT_SUCCESS;

} 

