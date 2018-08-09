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

#define MAXLEN 1024
#define	CUBE(x) (x*x*x)
#define ABSOLUTEMAXSNAPS 999

// Local Structs //

// Local Variables //

// Proto-types //

int32_t parse_params(int32_t argc, char **argv, struct SAGE_PARAMETERS *params);
int32_t read_snap_list(struct SAGE_PARAMETERS *params);
void myexit(int signum);

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
    return EXIT_FAILURE; 
  }
  
  status = read_parameter_file(argv[1]);
  if (status != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  } 
 
  return EXIT_SUCCESS;
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
#endif

  atexit(bye);

  parse_params(argc, argv);
  sage_init();
 
  if (self_consistent == 1 && (ReionizationOn == 3 || ReionizationOn == 4))
    status = init_selfcon_grid();
  if(ReionizationOn == 2)
  {
    status = init_grid();
    if (status == EXIT_FAILURE)
    {
      ABORT(EXIT_FAILURE);
    }
  }

  return EXIT_SUCCESS;
} 
