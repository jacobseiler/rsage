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

#include "read_parameter_file.h"
#include "tree_io.h"
#include "reionization.h"

#define MAXLEN 1024
#define	CUBE(x) (x*x*x)
#define ABSOLUTEMAXSNAPS 999

// Local Structs //

// Local Variables //

SAGE_params params;
int32_t SnapNum, first_run;
#ifdef MPI
int ThisTask, NTask, nodeNameLen;
char *ThisNode;
#endif

// Proto-types //

int32_t parse_params(int32_t argc, char **argv);
int32_t read_snap_list(SAGE_params params);

// Functions //

void bye()
{
#ifdef MPI
  MPI_Finalize();
  free(ThisNode);
#endif
}

int32_t parse_params(int32_t argc, char **argv)
{

  int32_t status;

  if (argc != 4)
  {
    fprintf(stderr, "\n\n"); 
    fprintf(stderr, "This code reads in photoionization grids and pre-processes the dark matter halo trees to create a list of filtering masses.\n");
    fprintf(stderr, "./create_filter mass <SAGE Input Parameter File> <Snapshot Number> <First Run Flag>\n");
    fprintf(stderr, "The snapshot number should be the snapshot we're creating the filtering masses for.\n");
    fprintf(stderr, "First run flag is 0 or 1 denoting whether this is the first snapshot we're creating the list for. Since we want one list for each tree file that contains all snapshots, the first time the code is run we need to create the file and the other times we simply need to append.\n");
    fprintf(stderr, "\n\n"); 
    return EXIT_FAILURE; 
  }

  status = read_parameter_file(argv[1], &params); 
  if (status == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  } 

  SnapNum = atoi(argv[2]);

  first_run = atoi(argv[3]);
  if (first_run < 0 || first_run > 1)
  {
    fprintf(stderr, "The first run flag can only be 0 or 1.\n");
    return EXIT_FAILURE;
  }

  printf("\n\n");
  printf("==================================================\n");
  printf("Executing with\nTree Directory: %s\nTree Name: %s\nNumber of Tree Files: %d\nPhotoionization Directory: %s\nPhotoionization Name: %s\nReionization Redshift Name: %s\nGridSize: %d\nSnapshot Number: %d\nFirst Run Flag: %d\n", params->TreeDir, params->TreeName, params->LastFile - params->FirstFile + 1, params->PhotoionDir, params->PhotoionName, params->ReionRedshiftName, params->GridSize, SnapNum, first_run); 
  printf("==================================================\n\n");
  printf("\n\n");
 
  return EXIT_SUCCESS;
}

int32_t read_snap_list(SAGE_params params) 
{

  FILE *fd;
  char fname[1000];
  double AA[ABSOLUTEMAXSNAPS];

  int32_t i, Snaplistlen = 0;

  snprintf(fname, MAXLEN, "%s", params->SnapListFile);

  printf("%s\n", fname);
  if(!(fd = fopen(fname, "r")))
  {
    printf("can't read output list in file '%s'\n", fname);
    return EXIT_FAILURE;  
  }

  Snaplistlen = 0;
  do
  {
    if(fscanf(fd, " %lg ", &AA[Snaplistlen]) == 1)
      Snaplistlen++;
    else
      break;
  }
  while(Snaplistlen < (params->LastSnapshotNr + 1));

  fclose(fd);

#ifdef MPI
  if(ThisTask == 0)
#endif
    printf("found %d defined times in snaplist\n", Snaplistlen);

  for (i = 0; i < Snaplistlen; ++i)
  {  
    params->ZZ[i] = 1 / AA[i] - 1;  
  } 

  return EXIT_SUCCESS;

}

int main(int argc, char **argv)
{

  int32_t status, filenr, Ntrees, totNHalos, *TreeNHalos, treenr, NHalos_Ionized = 0;
  int64_t *HaloID;
  float *reion_mod, sum_reion_mod; 

  halo_t Halos;
  grid_t Grid;

#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  ThisNode = malloc(MPI_MAX_PROCESSOR_NAME * sizeof(char));

  MPI_Get_processor_name(ThisNode, &nodeNameLen);
  if (nodeNameLen >= MPI_MAX_PROCESSOR_NAME)
  {
    printf("Node name string not long enough!...\n");
    exit(EXIT_FAILURE);
  }
#endif

  atexit(bye);
 
  status = parse_params(argc, argv); // Set the input parameters.
  if (status == EXIT_FAILURE)
  {
    exit(EXIT_FAILURE);
  }

  status = read_snap_list(params);
  if (status == EXIT_FAILURE)
  {
    exit(EXIT_FAILURE);
  }

  status = read_grid(SnapNum, params, &Grid); 
  if (status == EXIT_FAILURE)
  {
    exit(EXIT_FAILURE);
  }
  
  
  params->LastFile = 0;
#ifdef MPI  
  for(filenr = params->FirstFile + ThisTask; filenr < params->LastFile + 1 ; filenr += NTask)
#else
  for(filenr = params->FirstFile; filenr < params->LastFile + 1; filenr++)
#endif
  { 

    status = load_tree_table(filenr, params, &Ntrees, &totNHalos, &TreeNHalos);
    if (status == EXIT_FAILURE)
    {
      exit(EXIT_FAILURE);
    }       

    status = allocate_array_memory(totNHalos, &HaloID, &reion_mod);   
    if (status == EXIT_FAILURE)
    {
      exit(EXIT_FAILURE);
    }

    for (treenr = 0; treenr < Ntrees; ++treenr)
    {

      status = load_halos(treenr, TreeNHalos[treenr], &Halos);
      if (status == EXIT_FAILURE)
      {
        exit(EXIT_FAILURE);
      }

      status = populate_halo_arrays(filenr, treenr, TreeNHalos[treenr], SnapNum, Halos, Grid, params, &HaloID, &reion_mod, &NHalos_Ionized, &sum_reion_mod);
      if (status == EXIT_FAILURE)
      {
        exit(EXIT_FAILURE);
      }

      free(Halos);
    } 

    free_memory(&TreeNHalos, &HaloID, &reion_mod);
    
    printf("For file %d there were %d total halos within ionized regions (out of %d halos, a ratio of %.4f). The average ionization modifier for these is %.4f\n", filenr, NHalos_Ionized, totNHalos, (float)NHalos_Ionized / (float)totNHalos, sum_reion_mod / NHalos_Ionized);
    
  }

  status = free_grid(&Grid);
  if (status == EXIT_FAILURE)
  {
    exit(EXIT_FAILURE);
  }

  free(params); 
  return EXIT_SUCCESS;
} 
