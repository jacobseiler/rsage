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


char bufz0[1000], bufmergedz0[1000];
int exitfail = 1;

struct sigaction saveaction_XCPU;
volatile sig_atomic_t gotXCPU = 0;



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
    unlink(bufz0);

#ifdef MPI
    if(ThisTask == 0 && gotXCPU == 1)
      printf("Received XCPU, exiting. But we'll be back.\n");
#endif
	  }
}



int main(int argc, char **argv)
{
  int filenr, tree, halonr;
  struct sigaction current_XCPU;

  struct stat filestatus;
  FILE *fd;
  char tag[MAXLEN];

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

  if(argc != 2)
  {
    printf("\n  usage: sage <parameterfile>\n\n");
    ABORT(0);
  }

  atexit(bye);

  sigaction(SIGXCPU, NULL, &saveaction_XCPU);
  current_XCPU = saveaction_XCPU;
  current_XCPU.sa_handler = termination_handler;
  sigaction(SIGXCPU, &current_XCPU, NULL);

  read_parameter_file(argv[1]);
  init();

  if(ReionizationOn == 2 || ReionizationOn == 3)
  {
    status = init_grid();
    if (status == EXIT_FAILURE)
    {
      exit(EXIT_FAILURE);
    }
  }

  fprintf(stderr, "ListoutputSnaps[0] = %d\n", ListOutputSnaps[0]);

#ifdef MPI
  for(filenr = FirstFile+ThisTask; filenr <= LastFile; filenr += NTask)
#else
  for(filenr = FirstFile; filenr <= LastFile; filenr++)
#endif
  {
    
   //if(filenr == 33 || filenr == 49 || filenr == 80)
//	continue;

    if(use_tiamat == 1)
    {
      if(filenr < 10)
        snprintf(tag, MAXLEN, "00%d", filenr);
      else if (filenr >= 10 && filenr < 100)
        snprintf(tag, MAXLEN, "0%d", filenr);
      else
        snprintf(tag, MAXLEN, "%d", filenr);
      snprintf(bufz0, MAXLEN, "%s/%s_%s.dat", SimulationDir, TreeName, tag);
    }
    else
    {
      snprintf(bufz0, MAXLEN, "%s/%s.%d", SimulationDir, TreeName, filenr);
    }

    if(!(fd = fopen(bufz0, "r")))
    {
      printf("-- missing tree %s ... skipping\n", bufz0);
      continue;  // tree file does not exist, move along
    }
    else
      fclose(fd);


    sprintf(bufz0, "%s/%s_z%1.3f_%d", OutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[0]], filenr);
    if(stat(bufz0, &filestatus) == 0)
    {
      printf("-- output for tree %s already exists ... skipping\n", bufz0);
      continue;  // output seems to already exist, dont overwrite, move along
    }

    if((fd = fopen(bufz0, "w")))
      fclose(fd);

    sprintf(bufmergedz0, "%s/%s_MergedGalaxies_%d", OutputDir, FileNameGalaxies, filenr);
    
    if((fd = fopen(bufmergedz0, "w")))
      fclose(fd);

    FileNum = filenr;
    load_tree_table(filenr);    
    
    for(tree = 0; tree < Ntrees; tree++)
    {      
			assert(!gotXCPU);

      if(tree % 10000 == 0)
      {
#ifdef MPI
        printf("\ttask: %d\tnode: %s\tfile: %i\ttree: %i of %i\n", ThisTask, ThisNode, filenr, tree, Ntrees);
#else
				printf("\tfile: %i\ttree: %i of %i\n", filenr, tree, Ntrees);
#endif
				fflush(stdout);
      }

      TreeID = tree;
      load_tree(filenr, tree);

      gsl_rng_set(random_generator, filenr * 100000 + tree);
      NumGals = 0;
      GalaxyCounter = 0;
      MergedNr = 0;

      for(halonr = 0; halonr < TreeNHalos[tree]; halonr++)
        if(HaloAux[halonr].DoneFlag == 0)	
        construct_galaxies(halonr, tree, filenr);
   
      save_galaxies(filenr, tree);
      save_merged_galaxies(filenr, tree);    
      free_galaxies_and_tree();
      //break;            
    }

    finalize_galaxy_file(filenr);  
    finalize_merged_galaxy_file(filenr);
    
    free_tree_table();
    printf("\ndone file %d\n\n", filenr);
    //break; 
  }
  XASSERT((gal_mallocs == gal_frees) && (mergedgal_mallocs == mergedgal_frees), "We had %d Galaxy Mallocs and %d Galaxy Frees\n We had %d MergedGalaxy Mallocs and %d MergedGalaxy Frees.\n", gal_mallocs, gal_frees, mergedgal_mallocs, mergedgal_frees);
  exitfail = 0;
  printf("There was %d firstSF and %d notfirstSF\n", count_firstSF, count_notfirstSF);
  fprintf(stderr, "There was %d Star formation steps less than 50Myr and %d steps above 50Myr giving a ratio of %.4e \n", small_steps, large_steps, (float) large_steps/(float) small_steps); 
  fprintf(stderr, "Write grid array was called %d times.\n", times_written);
  fprintf(stderr, "There were %d galaxies outside box compared to %d inside\n", outside_box, inside_box);
  gsl_rng_free(random_generator); 

  fprintf(stderr, "Returned Mvir %d times compared to Len %d times\n", count_Mvir, count_Len);
  return 0;
}

