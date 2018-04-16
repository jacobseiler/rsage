#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>

#include "../core_allvars.h"
#include "../core_proto.h"
#include "tree_binary.h"

// Local Variables //

FILE *load_fd; 

// Local Proto-Types //

// External Functions //

void load_tree_table_binary(char *fname) 
{
  int32_t i, totNHalos;

  // We've already checked that the tree exists.

  load_fd = fopen(fname, "r");
  if (load_fd == NULL)
  {
    printf("Could not locate %s\n", fname);
    ABORT(EXIT_FAILURE);
  }

  myfread(&Ntrees, 1, sizeof(Ntrees), load_fd);
  myfread(&totNHalos, 1, sizeof(totNHalos), load_fd);

  TreeNHalos = mycalloc(Ntrees, sizeof(*(TreeNHalos)));  
  if (TreeNHalos == NULL)
  {
    fprintf(stderr, "Could not allocate memory for `TreeNHalos`.\n");
    ABORT(EXIT_FAILURE);
  }
  myfread(TreeNHalos, Ntrees, sizeof(*(TreeNHalos)), load_fd);

  TreeFirstHalo = mycalloc(Ntrees, sizeof(*(TreeFirstHalo)));
  if (TreeFirstHalo == NULL)
  {
    fprintf(stderr, "Could not allocate memory for `TreeFirstHalo`.\n");
    ABORT(EXIT_FAILURE);
  }

  if(Ntrees)
    TreeFirstHalo[0] = 0;
  for(i = 1; i < Ntrees; i++)
    TreeFirstHalo[i] = TreeFirstHalo[i - 1] + TreeNHalos[i - 1];

}

void load_tree_binary(int32_t treenr)
{
  // must have an FD
  assert(load_fd );

  Halo = mycalloc(TreeNHalos[treenr], sizeof(*(Halo)));
  if (Halo == NULL)
  {
    fprintf(stderr, "Could not allocate memory for `Halo`.\nTried to allocate memory for %d Halos,\n", TreeNHalos[treenr]);
    ABORT(EXIT_FAILURE);
  }

  myfread(Halo, TreeNHalos[treenr], sizeof(*(Halo)), load_fd);

}

void close_binary_file(void)
{
  if(load_fd) 
  {
    fclose(load_fd);
    load_fd = NULL;
  }
}
// Local Functions //

