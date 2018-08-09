#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"
#include "self_consistent/selfcon_grid.h"

#ifdef MPI
#include <mpi.h>
#endif

#define MEMORY_THRESHOLD 10.0 // If the allocated memory is under this value, we won't track it. 
// This prevents the user from being spammed by output messages for small allocs.

static double TotMem = 0, HighMarkMem = 0;

void track_memory(double size);

void *mymalloc(size_t size)
{

  void *ptr;

  ptr = malloc(size);
  if (ptr == NULL)
  {
    fprintf(stderr, "Could not allocate %.4f MB of memory.\n", size / (1024.0 * 1024.0));
    return NULL; 
  }

  track_memory(size);

  return ptr;
}

void *mycalloc(size_t n, size_t size)
{

  void *ptr;

  ptr = calloc(n, size);
  if (ptr == NULL)
  {
    fprintf(stderr, "Could not allocate %.4f MB of memory.\n", (n * size) / (1024.0 * 1024.0));
    return NULL; 
  }

  track_memory(n*size);

  return ptr;
}

void track_memory(double size)
{
 
  TotMem += size; 

  if (TotMem > HighMarkMem && (TotMem - HighMarkMem > (MEMORY_THRESHOLD * 1024.0 * 1024.0)))
  {
#ifdef MPI
    printf("New highmark memory (Task %d): %.4f MB (was %.4f MB)\n", ThisTask, TotMem / (1024.0 * 1024.0), HighMarkMem / (1024.0 * 1024.0));
#else
    printf("New highmark memory: %.4f MB (was %.4f MB)\n", TotMem / (1024.0 * 1024.0), HighMarkMem / (1024.0 * 1024.0));
#endif
    HighMarkMem = TotMem;
  }
}

void *myrealloc(void *p, size_t new_n, size_t old_n)
{
  
  void *newp = realloc(p, new_n);
  if(newp == NULL) 
  {
    printf("Failed to re-allocate memory for %.4f MB (old memory size was %.4f MB).\n",  new_n / (1024.0 * 1024.0), old_n / (1024.0 * 1024.0));
    ABORT(EXIT_FAILURE);
  }

  TotMem += (new_n - old_n);
 
  return newp; 
}

void myfree(void *p, size_t size)
{

  TotMem -= size; 
  free(p);

}

void print_final_memory(void)
{

#ifdef MPI

  double master_HighMarkMem = 0.0;

  MPI_Reduce(&HighMarkMem, &master_HighMarkMem, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (ThisTask == 0)
  { 
    printf("The highest peak memory usage was %.4f MB (%.4f GB)\n", master_HighMarkMem / (1024.0 * 1024.0), master_HighMarkMem / (1024.0 * 1024.0 * 1024.0));
    printf("Note: This is only updated every %.4f MB.\nTo change this, adjust the variable 'MEMORY_THRESHOLD' in `core_mymallo.c`.\n", MEMORY_THRESHOLD);
  }
#else
  printf("The highest peak memory usage was %.4f MB (%.4f GB)\n", HighMarkMem / (1024.0 * 1024.0), HighMarkMem / (1024.0 * 1024.0 * 1024.0));
  printf("Note: This is only updated every %.4f MB.\nTo change this, adjust the variable 'MEMORY_THRESHOLD' in `core_mymallo.c`.\n", MEMORY_THRESHOLD);

#endif

}

int32_t final_cleanup(char **argv)
{

  int32_t status;

  XASSERT((gal_mallocs == gal_frees) && (mergedgal_mallocs == mergedgal_frees), "We had %d Galaxy Mallocs and %d Galaxy Frees\n We had %d MergedGalaxy Mallocs and %d MergedGalaxy Frees.\n", gal_mallocs, gal_frees, mergedgal_mallocs, mergedgal_frees);  
 
  gsl_rng_free(random_generator); 

  if (ReionizationOn == 2 )
  {
    status = free_grid();
  } 

  // Copy the parameter file to the output directory. 
  char copy_command[MAX_STRING_LEN];
  snprintf(copy_command, MAX_STRING_LEN - 1, "cp %s %s", argv[1], OutputDir); 
  system(copy_command);
  
  if ((self_consistent == 1) && (ReionizationOn == 3 || ReionizationOn == 4))
  {
    status = save_selfcon_grid();
    if (status != EXIT_SUCCESS)
    {
      return EXIT_FAILURE; 
    }

    free_selfcon_grid(SelfConGrid);
  }

  if (IRA == 0)
  {
    free(coreburning_times);
    free(IMF_massgrid_eta);
    free(IMF_massgrid_m);
  }

  if (PhotonPrescription == 1)
  {
    free(stars_tbins);
    free(stars_Ngamma);
  }

  print_final_memory();
  printf("There was a total of %.4ee50 Photon/s emitted (intrinsic, not necessarily escaped) over the entire simulation.\n", Ngamma_HI_Total);


  return EXIT_SUCCESS;
}
