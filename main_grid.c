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

#include "core_allvars_grid.h"
#include "core_proto_grid.h"

char buf[MAXLEN];

struct sigaction saveaction_XCPU;
volatile sig_atomic_t gotXCPU = 0;
int exitfail = 1;

void termination_handler(int signum)
{
  gotXCPU = 1;
  sigaction(SIGXCPU, &saveaction_XCPU, NULL);
  if(saveaction_XCPU.sa_handler != NULL)
    (*saveaction_XCPU.sa_handler) (signum);
}

void myexit(int signum)
{
#ifdef MPI
  printf("Task: %d\tnode: %s\tis exiting\n\n\n", ThisTask, ThisNode);
#else
  printf("We're exiting\n\n\n");
#endif    
  exit(signum);
}

void bye()
{
#ifdef MPI
  MPI_Finalize();
  free(ThisNode);
#endif

  if(exitfail)
  {
#ifdef MPI
    if(ThisTask == 0 && gotXCPU == 1)
      printf("Received XCPU, exiting. But we'll be back.\n");
#endif
  }
}


int main(int argc, char **argv)
{

  struct sigaction current_XCPU;

  int32_t status;

#ifdef MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  ThisNode = malloc(MPI_MAX_PROCESSOR_NAME * sizeof(char));

  MPI_Get_processor_name(ThisNode, &nodeNameLen);
  if (nodeNameLen >= MPI_MAX_PROCESSOR_NAME)
  {
    printf("Node name string not long enough!...\n");
    ABORT(0);
  }
#endif

  atexit(bye);

  sigaction(SIGXCPU, NULL, &saveaction_XCPU);
  current_XCPU = saveaction_XCPU;
  current_XCPU.sa_handler = termination_handler;
  sigaction(SIGXCPU, &current_XCPU, NULL);


  int filenr, i;

  
  read_parameter_file(argv[1]);

  status = init(); // Initialize all the parameters (set units, create scale factor/age arrays etc).  
  if (status == EXIT_FAILURE)
  {
    exit(EXIT_FAILURE);
  } 
 
#ifdef MPI
  for(filenr = FirstFile+ThisTask; filenr <= LastFile; filenr += NTask)
#else
  for(filenr = FirstFile; filenr <= LastFile; filenr++)
#endif
  { 
#ifdef MPI     
    printf("I am Task %d and I'm doing file %d.\n", ThisTask, filenr);
#endif     

    for (i = 0; i < 2; ++i) // i = 0 does the normal galaxies, i = 1 does the merged galaxies.
    {
      if(i == 0)      
        snprintf(buf, MAXLEN, "%s/%s_z%1.3f_%d", GalaxiesInputDir, FileNameGalaxies, ZZ[LastSnapShotNr], filenr);
      else       
        snprintf(buf, MAXLEN, "%s/%s_MergedGalaxies_%d", GalaxiesInputDir, FileNameGalaxies, filenr);

      if ( access(buf, F_OK ) == -1) // Sanity check.
      {
        printf("-- input for file %s does not exist, exiting now.\n", buf);
        exit(EXIT_FAILURE); 
      }
      if (Verbose == 1)
      {
        printf("Loading galaxies for file %d, name '%s'\n", filenr, buf); 
      }
      status = load_gals(buf);    
      if (status == EXIT_FAILURE)
      {
        exit(EXIT_FAILURE);
      }

      status = update_grid_properties(filenr); // Go through each galaxy and read it's grid history and grid the properties.
      if (status == EXIT_FAILURE)
      {
        exit(EXIT_FAILURE);
      }
  

      if (fescPrescription != 1)
      {
        //status = save_fesc_properties(filenr, i);
        if (status == EXIT_FAILURE)
        {
          exit(EXIT_FAILURE);
        }
      } 

   
      free_gals();	

    } // Galaxy/Merger loop.

    printf("Done File %d.\n\n", filenr);

  } // File Loop.
  
#ifdef MPI
  
  struct GRID_STRUCT *master_grid;

  master_grid = MPI_sum_grids();
  if (master_grid == NULL)
  {
    exit(EXIT_FAILURE);
  }

  if (ThisTask == 0)
  {
    count_grid_properties(master_grid); // Counts how many halos/galaxies/Photons are in the grid at each redshift.
    status = save_grid(master_grid); // Saves grid.
  }
#else 
  count_grid_properties(Grid); // Counts how many halos/galaxies/Photons are in the grid at each redshift.
  status = save_grid(Grid); // Saves grid.
#endif
  if (status == EXIT_FAILURE)
  {
    exit(EXIT_FAILURE);
  }
 
  free_grid();

  exitfail = 0;
  gsl_rng_free(random_generator); 
  return 0;

} 
