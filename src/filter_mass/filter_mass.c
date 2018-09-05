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

#include "tree_io.h"
#include "reionization.h"
#include "save.h"

#define MAXLEN 1024
#define	CUBE(x) (x*x*x)
#define ABSOLUTEMAXSNAPS 999

// Local Structs //

// Local Variables //

// Proto-types //

// Functions //

int32_t filter_masses(char *FileNameGalaxies, char *TreeDir, char *TreeName, 
                      char *PhotoionDir, char *PhotoionName, char *ReionRedshiftName,
                      int32_t FirstFile, int32_t LastFile, int32_t GridSize, double BoxSize,
                      double Hubble_h, int32_t SnapNum, double Redshift, int32_t first_update_flag,
                      int32_t ThisTask, int32_t NTask)
{

  int32_t status, filenr, Ntrees, totNHalos, *TreeNHalos, treenr, NHalos_Ionized, NHalos_In_Regions, NHalos_ThisSnap;
  int64_t *HaloID;
  float *ReionMod, sum_ReionMod; 

  struct HALO_STRUCT *Halos; 
  struct GRID_STRUCT *Grid;

  struct SAGE_PARAMETERS *params;

  params = malloc(sizeof(struct SAGE_PARAMETERS));
  if (params == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for parameter struct in filter_mass.\n");
    return EXIT_FAILURE; 
  }

  params->FileNameGalaxies = FileNameGalaxies;
  params->TreeDir = TreeDir;
  params->TreeName = TreeName;
  params->PhotoionDir = PhotoionDir;
  params->PhotoionName = PhotoionName;
  params->ReionRedshiftName = ReionRedshiftName;
  params->FirstFile = FirstFile;
  params->LastFile = LastFile;
  params->BoxSize = BoxSize;
  params->Hubble_h = Hubble_h;  
  params->Redshift = Redshift;
  
  Grid = malloc(sizeof(struct GRID_STRUCT));
  if (Grid == NULL)
  {
    fprintf(stderr, "Could not allocate memory for Grid struct.\n");
    return EXIT_FAILURE; 
  }

  status = filter_mass_read_grid(SnapNum, first_update_flag, GridSize, BoxSize,
                                 PhotoionDir, ReionRedshiftName, PhotoionName, Grid, ThisTask);
  if (status != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }
   
#ifdef MPI 
  if (NTask > LastFile + 1)
  {
    fprintf(stderr, "Attempted to run filter_mass with %d Tasks.  However the last file is %d.\n", NTask, LastFile);
    return EXIT_FAILURE;
  }

  for(filenr = FirstFile + ThisTask; filenr < LastFile + 1 ; filenr += NTask)
#else
  for(filenr = FirstFile; filenr < LastFile + 1; filenr++)
#endif
  {
    NHalos_Ionized = 0;
    NHalos_In_Regions = 0;
    NHalos_ThisSnap = 0;
    
    status = filter_mass_load_tree_table(filenr, params, &Ntrees, &totNHalos, &TreeNHalos); // Loads the table for this file.
    if (status != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }       

#ifdef DEBUG_TREES
    int32_t i;
    int32_t maxhalos = -1, minhalos = 1e5;
    for (i = 0; i < Ntrees; ++i)
    {
      if (TreeNHalos[i] > maxhalos)
        maxhalos = TreeNHalos[i];
    
      if (TreeNHalos[i] < minhalos)
        minhalos = TreeNHalos[i];
    }
    printf("The smallest tree had %d halos and the largest had %d halos\n", minhalos, maxhalos);
#endif

    status = allocate_array_memory(totNHalos, &HaloID, &ReionMod); // Memory for the output arrays. 
    if (status != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }

    for (treenr = 0; treenr < Ntrees; ++treenr)
    {
      Halos = malloc(sizeof(struct HALO_STRUCT) * TreeNHalos[treenr]); // Allocate the memory in main to keep consistency.
      if (Halos == NULL)
      {
        fprintf(stderr, "Could not allocate memory for Halos in tree %d\n", treenr);
        return EXIT_FAILURE; 
      }

      status = load_halos(treenr, TreeNHalos[treenr], Halos); // Loads the halos for this tree.
      if (status != EXIT_SUCCESS)
      {
        return EXIT_FAILURE;
      }



      // Now time to go through all the halos in this tree, determine those at the Snapshot specified and the associate reionization modifier (if it's within an ionized cell).
      status = populate_halo_arrays(filenr, treenr, TreeNHalos[treenr], SnapNum, first_update_flag, Halos, Grid, params, &HaloID, &ReionMod, &NHalos_ThisSnap, &NHalos_Ionized, &NHalos_In_Regions, &sum_ReionMod);
      if (status != EXIT_SUCCESS)
      {
        return EXIT_FAILURE;
      }

      free(Halos);
    } 

#ifdef DEBUG
    printf("For file %d there were %d total halos within ionized regions (out of %d halos in this snapshot, a ratio of %.4f). There were %d total halos with a reionization modifier lower than 1.0 (a ratio of %.4f to the total number of halos in this snapshot). The average ionization modifier for these is %.4f\n", filenr, NHalos_In_Regions, NHalos_ThisSnap, (float)NHalos_In_Regions / (float)NHalos_ThisSnap, NHalos_Ionized, (float)NHalos_Ionized / (float)NHalos_ThisSnap, sum_ReionMod / NHalos_Ionized);
#endif

    status = save_arrays(HaloID, ReionMod, params, NHalos_Ionized, filenr, SnapNum, first_update_flag);
    if (status != EXIT_SUCCESS)
    {
      return EXIT_FAILURE;
    }

    free_memory(&TreeNHalos, &HaloID, &ReionMod);
  }

  // Everything done, time to free!

  status = filter_mass_free_grid(Grid);
  if (status != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }

  free(params); 
  return EXIT_SUCCESS;
} 
