#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>
#include <time.h>

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

// Misc Declarations //
#include "check_ini_files.h"


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
  double SAGE_time = 0.0, cifog_time = 0.0, filter_time = 0.0, misc_time = 0.0;
  clock_t before;

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

  // Some ini file parameters are allowed to be specifically ``None``.
  // For these, we want to update their values depending upon the recipe flags
  // and constants chosen.
  status = check_ini_files(FileNameGalaxies, OutputDir, GalaxyOutputDir, GridOutputDir,
                           PhotoionDir, PhotoionName, ReionRedshiftName, argv, ThisTask);
  if (status !=  EXIT_SUCCESS)
  {
    goto err;
  }

#ifdef RSAGE

  if (self_consistent == 1 && (ReionizationOn == 3 || ReionizationOn == 4))
  {
    int32_t loop_SnapNum, first_update_flag, count; 

    for (loop_SnapNum = LowSnap, count = 0; loop_SnapNum < LastSnapShotNr; ++loop_SnapNum, ++count)
    {
      if (ThisTask == 0)
      {
        printf("\n\n======================================\n");
        printf("(Rank %d) Applying reionization to Snapshot %d\n", ThisTask, loop_SnapNum); 
        printf("======================================\n\n");
      }
      ReionSnap = loop_SnapNum;

      if (ThisTask == 0)
      {
        printf("\n\n================================\n");
        printf("(Rank %d) Running SAGE\n", ThisTask);
        printf("================================\n");
      }

      before = clock();
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
      SAGE_time = SAGE_time + (clock() - before); 
 
      before = clock(); 
      status = cifog_zero_grids(grid, simParam);
      if (status !=  EXIT_SUCCESS)
      {
        goto err;
      }

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

      if (ThisTask == 0)
      {
        printf("\n\n================================\n");
        printf("(Rank %d) Running cifog\n", ThisTask);
        printf("================================\n");
      }

      num_cycles = count + 1;
      status = cifog(simParam, redshift_list, grid, sourcelist, integralTable, photIonBgList, num_cycles, ThisTask, RestartMode);
      if (status !=  EXIT_SUCCESS)
      {
        goto err;
      }

#ifdef MPI
      MPI_Barrier(MPI_COMM_WORLD); 
#endif
      cifog_time = cifog_time + (clock() - before);

      // Because of how cifog handles numbering, need to pass the Snapshot Number + 1.
      // This file only needs to be created once, so do it on Root node. 
      if (ThisTask == 0)
      {        
        printf("\n\n================================\n");
        printf("(Rank %d) Updating the Reionization Redshift file.\n", ThisTask);
        printf("================================\n");

        before = clock();
        status = update_reion_redshift(loop_SnapNum+1, ZZ[loop_SnapNum], GridSize, first_update_flag,
                                       PhotoionDir, FileNameGalaxies, ReionRedshiftName);
        if (status !=  EXIT_SUCCESS)
        {
          goto err;
        }
        misc_time = misc_time + (clock() - before);

      }
      
#ifdef MPI
      MPI_Barrier(MPI_COMM_WORLD); 
#endif

      if (ThisTask == 0)
      {
        printf("\n\n================================\n");
        printf("(Rank %d) Updating the filter mass values.\n", ThisTask);
        printf("================================\n");
      }

      before = clock();
      status = filter_masses(FileNameGalaxies, SimulationDir, TreeName, PhotoionDir, PhotoionName, ReionRedshiftName,
                             FirstFile, LastFile, GridSize, BoxSize, Hubble_h, loop_SnapNum, ZZ[loop_SnapNum], 
                             first_update_flag, ThisTask, NTask);
      filter_time = filter_time + (clock() - before);

      if (status !=  EXIT_SUCCESS)
      {
        goto err;
      }

    } // Self-Consistent Snapshot Loop.
  } // ReionizationOption condition.
  else
  {
    before = clock();
    sage();
    SAGE_time = clock() - before;
  }
#else
  before = clock();
  sage();
  SAGE_time = clock() - before;
#endif

  printf("About to cleanup\n");
  status = sage_cleanup();
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

  if (ThisTask == 0)
  {
    printf("\n\n================================\n");
    printf("Runtime Statistics\n");
    printf("================================\n");
  }

#ifdef MPI
  double master_SAGE_time, master_cifog_time, master_filter_time, master_misc_time;

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Reduce(&SAGE_time, &master_SAGE_time, NTask, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
  MPI_Reduce(&cifog_time, &master_cifog_time, NTask, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
  MPI_Reduce(&filter_time, &master_filter_time, NTask, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
  MPI_Reduce(&misc_time, &master_misc_time, NTask, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 

  if (ThisTask == 0)
  {
    master_SAGE_time = master_SAGE_time / NTask;
    master_cifog_time = master_cifog_time / NTask;
    master_filter_time = master_filter_time / NTask;
    master_misc_time = master_misc_time / NTask;
  }

#else
  double master_SAGE_time = SAGE_time;
  double master_cifog_time = cifog_time;
  double master_filter_time = filter_time;
  double master_misc_time = misc_time;
#endif

  if (ThisTask == 0)
  {
    printf("SAGE took an average time of %.4f seconds to execute\n", master_SAGE_time / CLOCKS_PER_SEC);
    printf("cifog took an average time of %.4f seconds to execute\n", master_cifog_time / CLOCKS_PER_SEC);
    printf("Creation of Halo filter masses took an average time of %.4f seconds to execute\n", master_filter_time / CLOCKS_PER_SEC);
    printf("Creation of reionization redshift grid took an average time of %.4f seconds to execute\n", master_misc_time / CLOCKS_PER_SEC * NTask); // Reionization redshift only done on root node.
  }

#ifdef MPI
  MPI_Finalize();
#endif

  return EXIT_SUCCESS;

 err:    
#ifdef MPI        
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    MPI_Finalize();
#endif
    return status;
} 
