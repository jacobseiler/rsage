#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "sage/core_allvars.h"
#include "sage/core_proto.h"

#ifdef RSAGE
#ifdef MPI
#include <fftw3-mpi.h>
#else
#include <fftw3.h>
#endif

#include "sage/self_consistent/selfcon_grid.h"

#include "../grid-model/src/confObj.h"
#include "../grid-model/src/grid.h"
#include "../grid-model/src/photion_background.h"
#include "../grid-model/src/sources.h"
#include "../grid-model/src/recombination.h"

#include "../grid-model/src/init.h"
#include "../grid-model/src/cifog.h"
#endif

#define MAXLEN 1024
#define	CUBE(x) (x*x*x)

// Local Structs //

// Local Variables //

// Proto-types //

/*
void myexit(int signum);
*/
// Functions //

void my_bye()
{
#ifdef MPI
  MPI_Finalize();
  free(ThisNode);
#endif
}

void main_myexit(int signum)
{
#ifdef MPI
  fprintf(stderr, "Task: %d\tnode: %s\tis exiting\n\n\n", ThisTask, ThisNode);
  MPI_Abort(MPI_COMM_WORLD, signum);
#else
  fprintf(stderr, "We're exiting\n\n\n");
	exit(signum);
#endif

}

int main(int argc, char **argv)
{

#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  ThisNode = malloc(MPI_MAX_PROCESSOR_NAME * sizeof(char));

  MPI_Get_processor_name(ThisNode, &nodeNameLen);
  if (nodeNameLen >= MPI_MAX_PROCESSOR_NAME)
  {
    printf("Node name string not long enough!...\n");
    ABORT(EXIT_FAILURE);
  }
#else
  ThisTask = 1;
  NTask = 1;
#endif

  int32_t status;

  atexit(my_bye);

#ifdef RSAGE
  if(argc != 3)
  {
    printf("\n  usage: rsage <sage_parameterfile> <cifog_parameterfile>\n\n");
    ABORT(EXIT_FAILURE);
  }
#else
  if(argc != 2)
  {
    printf("\n  usage: rsage <sage_parameterfile>\n\n");
    ABORT(EXIT_FAILURE);
  }
#endif

  status = read_parameter_file(argv[1]);
  if (status == EXIT_FAILURE)
  {
    ABORT(EXIT_FAILURE);
  }

  printf("Initing Sage\n");
  sage_init();

#ifdef RSAGE
  // First initialize the grid for the self-consistent part of SAGE.
  if (self_consistent == 1 && (ReionizationOn == 3 || ReionizationOn == 4))
  {
    status = init_selfcon_grid();
    if (status != EXIT_SUCCESS)
    {
      ABORT(EXIT_FAILURE);
    }
  
    status = zero_selfcon_grid(SelfConGrid);
    if (status != EXIT_SUCCESS)
    {
      ABORT(EXIT_FAILURE);
    }
  }

  // Then initialize all the variables for cifog.
  //int32_t RestartMode, num_cycles = 1;
  int32_t num_cycles = 1;
  double *redshift_list = NULL;

  confObj_t simParam;  
  grid_t *grid = NULL;
  //sourcelist_t *sourcelist = NULL;
  integral_table_t *integralTable = NULL;
  photIonlist_t *photIonBgList = NULL;
  
  status = init_cifog(argv[2], &simParam, &redshift_list, &grid, &integralTable, &photIonBgList, &num_cycles, ThisTask);
  if (status !=  EXIT_SUCCESS)
  {
    exit(EXIT_FAILURE);
  }

#endif

#ifdef RSAGE
  int32_t SnapNum; 

  printf("LowSnap %d\n", LowSnap);
  printf("HighSnap %d\n", HighSnap);
  for (SnapNum = LowSnap; SnapNum < HighSnap + 1; ++SnapNum)
  {
    printf("Running SAGE\n");
    sage();

    /*
    cifog_zero_grids(grid, simParam);

    // When we run cifog we want to read the output of the previous snapshot and save it at the end.
    // For the first snapshot we only save, otherwise we read and save.
    if (SnapNum == LowSnap)
    {
      RestartMode = 1;
    }
    else
    {
      RestartMode = 3;
    }
    cifog(simParam, redshift_list, grid, sourcelist, integralTable, photIonBgList, num_cycles, ThisTask, RestartMode);
    */
  } 

#else
  sage();
#endif
  //parse_params(argc, argv);

  // First read the snapshot list and initalize the lookup tables for delayedSN
  // and ionizing photon prescription (if requested).
  /*
  sage_init();


  for (snapshot that we're looping over)
  {
    if (self_consistent == 1 && (ReionizationOn == 3 || ReionizationOn == 4))
      status = init_selfcon_grid();
    if(ReionizationOn == 2)
    {
      status = init_grid();
      if (status != EXIT_SUCCESS)
      {
        ABORT(EXIT_FAILURE);
      }
    }
    sage;
    cifog;
    gen_reion_redshift;
    gen_halolists;
   
  }
  */

  status = sage_cleanup(argv);
  if (status != EXIT_SUCCESS)
  {
    ABORT(EXIT_FAILURE);
  }

#ifdef RSAGE
  status = cleanup_cifog(simParam, integralTable, photIonBgList, grid, redshift_list, ThisTask);
  if (status !=  EXIT_SUCCESS)
  {
    exit(EXIT_FAILURE);
  }

#endif

  return EXIT_SUCCESS;
} 
