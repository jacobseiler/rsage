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

// Local Variables //

char *treedir, *treename, *photoiondir, *photoionname, *reionredshiftname;
int32_t snapnum, first_run, numfiles;
#ifdef MPI
int ThisTask, NTask, nodeNameLen;
char *ThisNode;
#endif

// Proto-types //

int32_t parse_params(int32_t argc, char **argv);
int32_t init(void);

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

  status = read_parameter_file(argv[1], &treedir, &treename, &photoiondir, &photoionname, &reionredshiftname, &numfiles);  
  if (status == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  } 

  snapnum = atoi(argv[2]);

  first_run = atoi(argv[3]);
  if (first_run < 0 || first_run > 1)
  {
    fprintf(stderr, "The first run flag can only be 0 or 1.\n");
    return EXIT_FAILURE;
  }

  printf("\n\n");
  printf("==================================================\n");
  printf("Executing with\nTree Directory: %s\nTree Name: %s\nNumber of Tree Files: %d\nPhotoionization Directory: %s\nPhotoionization Name: %s\nReionization Redshift Name: %s\nSnapshot Number: %d\nFirst Run Flag: %d\n", treedir, treename, numfiles, photoiondir, photoionname, reionredshiftname, snapnum, first_run); 
  printf("==================================================\n\n");
  printf("\n\n");
 
  return EXIT_SUCCESS;
}

int32_t init()
{

  return EXIT_SUCCESS;

}

int main(int argc, char **argv)
{

  int32_t status, filenr, Ntrees, totNHalos, *TreeNHalos, treenr, NHalos_ThisSnap = 0;
  int64_t *HaloID;
  float *Mfilt; 

  Halo_t Halos;
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

  status = read_grid(snapnum, photoiondir, photoionname, &Grid); 
  
  numfiles = 1;
#ifdef MPI  
  for(filenr = 0 + ThisTask; filenr < numfiles; filenr += NTask)
#else
  for(filenr = 0; filenr < numfiles; filenr++)
#endif
  { 

    status = load_tree_table(filenr, treedir, treename, &Ntrees, &totNHalos, &TreeNHalos);
    if (status == EXIT_FAILURE)
    {
      exit(EXIT_FAILURE);
    }       

    status = allocate_array_memory(totNHalos, &HaloID, &Mfilt);   
    if (status == EXIT_FAILURE)
    {
      exit(EXIT_FAILURE);
    }

    for (treenr = 0; treenr < Ntrees; ++treenr)
    {

      status = load_halos(treenr, TreeNHalos[treenr], &Halos, &NHalos_ThisSnap);
      if (status == EXIT_FAILURE)
      {
        exit(EXIT_FAILURE);
      }

      free(Halos);
    } 

    free_memory(&TreeNHalos, &HaloID, &Mfilt);
        
  }


 
  return EXIT_SUCCESS;
} 
