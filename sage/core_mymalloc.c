#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"

#define MAXBLOCKS 256

static size_t TotMem = 0, HighMarkMem = 0;

void *mymalloc(size_t n)
{

  void *ptr;

  ptr = malloc(n);
  if (ptr == NULL)
  {
    fprintf(stderr, "Could not allocate %g MB of memory.\n", n / (1024.0 * 1024.0));
    ABORT(EXIT_FAILURE);
  }

  TotMem += n; 

  if (TotMem > HighMarkMem)
  {
    printf("New highmark memory: %g MB (was %g MB)\n", TotMem / (1024.0 * 1024.0), HighMarkMem / (1024.0 * 1024.0));
    HighMarkMem = TotMem;
  } 
  return ptr; 
}

void *myrealloc(void *p, size_t new_n, size_t old_n)
{
  
  void *newp = realloc(p, new_n);
  if(newp == NULL) 
  {
    printf("Failed to re-allocate memory for %g MB (old memory size was %g MB).\n",  new_n / (1024.0 * 1024.0), old_n / (1024.0 * 1024.0));
    ABORT(EXIT_FAILURE);
  }

  TotMem += (new_n - old_n);
 
  return newp; 
}

void myfree(void *p, size_t n)
{

  TotMem -= n;
  free(p);

}
