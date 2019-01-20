#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "core_allvars.h"
#include "core_proto.h"
#include "self_consistent/selfcon_grid.h"

#define TREEPROGRESSBAR 100000 // Number of trees processed before printing out status.

char bufz0[1000], bufmergedz0[1000];
int exitfail = 1;

struct sigaction saveaction_XCPU;
volatile sig_atomic_t gotXCPU = 0;

/* Local Proto-Types */

int32_t check_tree_file(int32_t filenr, int32_t *treestyle);

/**/

void bye()
{
#ifdef MPI
  MPI_Finalize();
  free(ThisNode);
#endif
}

void myexit(int signum)
{
#ifdef MPI
  fprintf(stderr, "Task: %d\tnode: %s\tis exiting\n\n\n", ThisTask, ThisNode);
  MPI_Abort(MPI_COMM_WORLD, signum);
#else
  fprintf(stderr, "We're exiting\n\n\n");
	exit(signum);
#endif

}

int32_t check_tree_file(int32_t filenr, int32_t *treestyle)
{
  struct stat filestatus;

  snprintf(bufz0, MAXLEN, "%s/%s_%03d%s", SimulationDir, TreeName, filenr, TreeExtension);
  if (stat(bufz0, &filestatus) != 0)
  {
    printf("Could not locate %s\n", bufz0);
  }
  else
  {    
    *treestyle = 0;
    return EXIT_SUCCESS;
  }  

  snprintf(bufz0, MAXLEN, "%s/%s.%d%s", SimulationDir, TreeName, filenr, TreeExtension);
  if (stat(bufz0, &filestatus) != 0)
  {
    printf("Could not locate %s either. Skipping this tree!\n", bufz0);
    return EXIT_FAILURE; 
  }
  else
  {
    *treestyle = 1;
    return EXIT_SUCCESS;
  }  

}

int32_t sage(void)
{

  int32_t filenr, status, treestyle, tree, halonr;
  struct stat filestatus;

#ifdef MPI
  for(filenr = FirstFile+ThisTask; filenr <= LastFile; filenr += NTask)
#else
  for(filenr = FirstFile; filenr <= LastFile; filenr++)
#endif
  {
    if (ReionizationOn == 3 || ReionizationOn == 4)
    {
      status = init_reion_lists(filenr);
      if (status == EXIT_FAILURE)
      {
        return EXIT_FAILURE;
      }
    }
      
    // Some tree files are padded to three digits whereas some are not (e.g., tree.4 vs tree.004)
    // Try both styles to see if the tree can be located.  If not, skip this tree. 
    status = check_tree_file(filenr, &treestyle);
    if (status != EXIT_SUCCESS) 
    {
      continue;
    }

    if (treestyle == 0)
      sprintf(bufz0, "%s/%s_z%1.3f_%03d", GalaxyOutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[0]], filenr);
    else
      sprintf(bufz0, "%s/%s_z%1.3f.%d", GalaxyOutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[0]], filenr);
  
    if(stat(bufz0, &filestatus) == 0 && self_consistent == 0)
    {
      printf("-- output for tree %s already exists ... skipping\n", bufz0);
      continue;  // output seems to already exist, dont overwrite, move along
    }

    if (treestyle == 0)
      sprintf(bufmergedz0, "%s/%s_MergedGalaxies_%03d", GalaxyOutputDir, FileNameGalaxies, filenr);
    else
      sprintf(bufmergedz0, "%s/%s_MergedGalaxies.%d", GalaxyOutputDir, FileNameGalaxies, filenr);
    
    FileNum = filenr;
    load_tree_table(filenr, treestyle);
    
    for(tree = 0; tree < Ntrees; tree++)
    {      
			assert(!gotXCPU);

      if((tree+1) % TREEPROGRESSBAR == 0)
      {
#ifdef MPI
        printf("\ttask: %d\tnode: %s\tfile: %i\ttree: %i of %i\n", ThisTask, ThisNode, filenr, tree, Ntrees);
#else
				printf("\tfile: %d\ttree: %d of %d\n", filenr, tree, Ntrees);
#endif
				fflush(stdout);
      }

      TreeID = tree;
      load_tree(tree);
      gsl_rng_set(random_generator, filenr * 100000 + tree);
      NumGals = 0;
      GalaxyCounter = 0;
      MergedNr = 0;

      for(halonr = 0; halonr < TreeNHalos[tree]; halonr++)
        if(HaloAux[halonr].DoneFlag == 0)	
        construct_galaxies(halonr, tree, filenr);

      save_galaxies(filenr, tree);
      save_merged_galaxies(filenr, tree);    
      free_galaxies_and_tree(tree);      
    }

    finalize_galaxy_file();  
    finalize_merged_galaxy_file();
    
    free_tree_table();

    printf("Done file %d\n", filenr);

    if (self_consistent == 1 && (ReionizationOn == 3 || ReionizationOn == 4))
    {
      status = free_reion_lists(filenr);
      if (status != EXIT_SUCCESS)
      {
        return EXIT_FAILURE;
      }
    }
  } // Filenr loop

  if ((self_consistent == 1) && (ReionizationOn == 3 || ReionizationOn == 4))
  {
    status = save_selfcon_grid();
    if (status != EXIT_SUCCESS)
    {
      return EXIT_FAILURE; 
    }
  }

  return EXIT_SUCCESS;
}
