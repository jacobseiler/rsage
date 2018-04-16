#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"

#ifdef MPI
#include <mpi.h>
#endif

#define MEMORY_THRESHOLD 10.0 // If the allocated memory is under this value, we won't track it. 
// This prevents the user from being spammed by output messages for small allocs.

static double TotMem = 0, HighMarkMem = 0;

void *mycalloc(size_t n, size_t size)
{

  void *ptr;

  ptr = calloc(n, size);
  if (ptr == NULL)
  {
    fprintf(stderr, "Could not allocate %.4f MB of memory.\n", n / (1024.0 * 1024.0));
    return NULL; 
  }

  TotMem += n * size; 

  if (TotMem > HighMarkMem && (TotMem - HighMarkMem > (MEMORY_THRESHOLD * 1024.0 * 1024.0)))
  {
#ifdef MPI
    printf("New highmark memory (Task %d): %.4f MB (was %.4f MB)\n", ThisTask, TotMem / (1024.0 * 1024.0), HighMarkMem / (1024.0 * 1024.0));
#else
    printf("New highmark memory: %.4f MB (was %.4f MB)\n", TotMem / (1024.0 * 1024.0), HighMarkMem / (1024.0 * 1024.0));
#endif
    HighMarkMem = TotMem;
  } 
  return ptr; 
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
    printf("The highest peak memory usage was %.4f MB (%.4f GB)\n", master_HighMarkMem / (1024.0 * 1024.0), master_HighMarkMem / (1024.0 * 1024.0 * 1024.0));
#else
  printf("The highest peak memory usage was %.4f MB (%.4f GB)\n", HighMarkMem / (1024.0 * 1024.0), HighMarkMem / (1024.0 * 1024.0 * 1024.0));
  printf("Note: Memory allocations below %.4f MB are not tracked.\nTo change this, adjust the variable 'MEMORY_THRESHOLD' in `core_mymallo.c`.\n", MEMORY_THRESHOLD);

#endif

}
