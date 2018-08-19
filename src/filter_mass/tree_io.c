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

int32_t filter_mass_load_tree_table(int32_t filenr, struct SAGE_PARAMETERS *params, int32_t *Ntrees, int32_t *totNHalos, int32_t **TreeNHalos)
{
  char forestname[MAXLEN]; 

  snprintf(forestname, MAXLEN, "%s/%s_%03d.dat", params->TreeDir, params->TreeName, filenr);

  if(!(forest_file= fopen(forestname, "rb")))
  {
    printf("can't open file `%s'\n", forestname);
    return EXIT_FAILURE; 
  }

  fread(Ntrees, 1, sizeof(int32_t), forest_file);
  fread(totNHalos, 1, sizeof(int32_t), forest_file);

  *TreeNHalos = malloc(sizeof(int) * (*Ntrees));
  if ((*TreeNHalos) == NULL)
  {
    fprintf(stderr, "Could not allocate memory for TreeNHalos.\n");
    return EXIT_FAILURE;
  }
  
  fread(*TreeNHalos, (*Ntrees), sizeof(int), forest_file); 
  
  return EXIT_SUCCESS;

}

int32_t load_halos(int32_t treenr, int32_t NHalos_ThisTree, struct HALO_STRUCT *Halos)
{
  int32_t Halos_Read;
 
  Halos_Read = fread(Halos, sizeof(struct HALO_STRUCT), NHalos_ThisTree, forest_file);  
  if (NHalos_ThisTree != Halos_Read)
  {
    fprintf(stderr, "I tried to read %d halos but I only read %d\n", NHalos_ThisTree, Halos_Read);
    return EXIT_FAILURE;
  }

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
    printf("Closing forest file\n");
  }

  return EXIT_SUCCESS;

}

