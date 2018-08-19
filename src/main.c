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

// cifog Declarations //
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

// Reionization Redshift Declarations //
#include "reion_redshift.h"

// Filter Mass Declarations //
#include "filter_mass/filter_mass.h"

#endif

#define MAXLEN 1024

// Local Structs //

// Local Variables //

// Proto-types //

// Functions //

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
  ThisTask = 0;
  NTask = 1;
#endif

  int32_t status;

#ifdef RSAGE
  if(argc != 3)
  {
    fprintf(stderr, "usage: rsage <sage_parameterfile> <cifog_parameterfile>\n");
    fprintf(stderr, "If you wish to only run SAGE without self-consistent reionization, recompile "
                    "with `BUILD_RSAGE` undefined.\n"); 
    ABORT(EXIT_FAILURE);
  }
  else
  {
    if (ThisTask == 0)
    {
      printf("\n\n======================================\n");
      printf("Running Full RSAGE Pipeline\n");
      printf("SAGE Parameter file: %s\n", argv[1]); 
      printf("cifog Parameter file: %s\n", argv[2]); 
      printf("Git Version: %s\n", VERSION);
      printf("======================================\n\n");
    }

  }
#else
  if(argc != 2)
  {
    printf("\n  usage: rsage <sage_parameterfile>\n\n");
    ABORT(EXIT_FAILURE);
  }
  else
  {
    if (ThisTask == 0)
    {
      printf("\n\n======================================\n");
      printf("Running SAGE\n");
      printf("SAGE Parameter file: %s\n", argv[1]);
      printf("Git Version: %s\n", VERSION);
      printf("======================================\n\n");
    }
  }
#endif

  status = read_parameter_file(argv[1]);
  if (status == EXIT_FAILURE)
  {
    ABORT(EXIT_FAILURE);
  }  
  sage_init();

#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD); 
#endif

#ifdef RSAGE
  // First initialize the grid for the self-consistent part of SAGE.
  if (self_consistent == 1 && (ReionizationOn == 3 || ReionizationOn == 4))
  {
    status = init_selfcon_grid();
    if (status != EXIT_SUCCESS)
    {
      ABORT(EXIT_FAILURE);
    }
  }

  // Then initialize all the variables for cifog.
  int32_t RestartMode, num_cycles;
  double *redshift_list = NULL;

  confObj_t simParam;  
  grid_t *grid = NULL;
  sourcelist_t *sourcelist = NULL;
  integral_table_t *integralTable = NULL;
  photIonlist_t *photIonBgList = NULL;

  if (self_consistent == 1 && (ReionizationOn == 3 || ReionizationOn == 4))
  {
    status = init_cifog(argv[2], &simParam, &redshift_list, &grid, &integralTable, &photIonBgList, &num_cycles, ThisTask);
    if (status !=  EXIT_SUCCESS)
    {
      goto err;
    }
  }
#endif

#ifdef RSAGE

  if (self_consistent == 1 && (ReionizationOn == 3 || ReionizationOn == 4))
  {
    int32_t loop_SnapNum, first_update_flag, count; 

    for (loop_SnapNum = LowSnap, count = 0; loop_SnapNum < LastSnapShotNr + 1; ++loop_SnapNum, ++count)
    {
      printf("\n\n======================================\n");
      printf("(Rank %d) Applying reionization to Snapshot %d\n", ThisTask, loop_SnapNum); 
      printf("======================================\n\n");
      ReionSnap = loop_SnapNum;

      printf("\n\n================================\n");
      printf("(Rank %d) Running SAGE\n", ThisTask);
      printf("================================\n");

      status = zero_selfcon_grid(SelfConGrid);
      if (status != EXIT_SUCCESS)
      {
        goto err; 
      }

      status = sage();
      if (status !=  EXIT_SUCCESS)
      {
        goto err;
      }

      status = cifog_zero_grids(grid, simParam);
      if (status !=  EXIT_SUCCESS)
      {
        goto err;
      }

#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD); 
#endif

      // When we run cifog we want to read the output of the previous snapshot and save it at the end.
      // For the first snapshot we only save, otherwise we read and save.
      if (loop_SnapNum == LowSnap)
      {
        RestartMode = 1;
        first_update_flag = 1;
      }
      else
      {
        RestartMode = 3;
        first_update_flag = 0;
      }

      printf("\n\n================================\n");
      printf("(Rank %d) Running cifog\n", ThisTask);
      printf("================================\n");
      num_cycles = count + 1;
      status = cifog(simParam, redshift_list, grid, sourcelist, integralTable, photIonBgList, num_cycles, ThisTask, RestartMode);
      if (status !=  EXIT_SUCCESS)
      {
        goto err;
      }

      // Because of how cifog handles numbering, need to pass the Snapshot Number + 1.
      // This file only needs to be created once, so do it on Root node. 
      if (ThisTask == 0)
      {        
        printf("\n\n================================\n");
        printf("(Rank %d) Updating the Reionization Redshift file.\n", ThisTask);
        printf("================================\n");
        status = update_reion_redshift(loop_SnapNum+1, ZZ[loop_SnapNum], GridSize, first_update_flag,
                                       PhotoionDir, FileNameGalaxies, ReionRedshiftName);

        if (status !=  EXIT_SUCCESS)
        {
          goto err;
        }
      }
      
#ifdef MPI
      MPI_Barrier(MPI_COMM_WORLD); 
#endif

      printf("\n\n================================\n");
      printf("(Rank %d) Updating the filter mass values.\n", ThisTask);
      printf("================================\n");
      status = filter_masses(FileNameGalaxies, SimulationDir, TreeName, PhotoionDir, PhotoionName, ReionRedshiftName,
                             FirstFile, LastFile, GridSize, BoxSize, Hubble_h, loop_SnapNum, ZZ[loop_SnapNum], 
                             first_update_flag);

      if (status !=  EXIT_SUCCESS)
      {
        goto err;
      }

    } // Self-Consistent Snapshot Loop.
  } // ReionizationOption condition.
  else
    sage();
#else
  sage();
#endif

  status = sage_cleanup(argv);
  if (status != EXIT_SUCCESS)
  {
    goto err;
  }

#ifdef RSAGE

  if (self_consistent == 1 && (ReionizationOn == 3 || ReionizationOn == 4))
  {
    status = cleanup_cifog(simParam, integralTable, photIonBgList, grid, redshift_list, ThisTask);
    if (status !=  EXIT_SUCCESS)
    {
      goto err;
    }
  }

#endif

   err:    
#ifdef MPI        
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    MPI_Finalize();
#endif

  return status;
} 
