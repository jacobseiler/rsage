#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>

#include "../core_allvars.h"
#include "../core_proto.h"
#include "selfcon.h"

// Local Variables //

// Local Proto-Types //

// External Functions //

int32_t init_grid()
{

  // There are two modes of operation for accounting for the photoionization feedback. //
  // In the first we read the photoionization rates and redshift of reionization for ALL redshifts. //
  // Feedback is then applied to all of redshifts.  This has extensive memory requirements as assuming a 512^3 grid to float precision, then each grid will be ~0.5GB in size. //

  int32_t i; 
  FILE *photoion, *reionredshift;
  char buf[MAXLEN];

  printf("Reading in the data for the reionization feedback using cifog.\n");
  Grid = calloc(1, sizeof(*(Grid)));
  if (Grid == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for the high level grid struct.\n");
    return EXIT_FAILURE;
  }

  Grid->GridSize = GridSize;
  Grid->NumCellsTotal = CUBE(GridSize); 

  Grid->ReionRedshift = mycalloc(Grid->NumCellsTotal, sizeof(*(Grid->ReionRedshift))); 
  if (Grid->ReionRedshift == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for the reionization redshift grid.\n");
    return EXIT_FAILURE;
  } 

  snprintf(buf, MAXLEN, "%s/%s", PhotoionDir, ReionRedshiftName); 
  if(!(reionredshift= fopen(buf, "rb")))
  {
    fprintf(stderr, "Cannot open file %s\n", buf);
    return EXIT_FAILURE;
  }

  printf("Reading the reionization redshift grid from %s\n", buf);
  fread(Grid->ReionRedshift, sizeof(*(Grid->ReionRedshift)), Grid->NumCellsTotal, reionredshift);
  fclose(reionredshift);
   
  Grid->NumGrids = MAXSNAPS;
    
  Grid->PhotoGrid = malloc(sizeof(struct PHOTO_GRID) * Grid->NumGrids); // Allocate enough memory to hold the photoionization grids.
  if (Grid->PhotoGrid == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for the photoionization grid struct\n");
    return EXIT_FAILURE;
  }

  for (i = 0; i < Grid->NumGrids; ++i)
  {
    Grid->PhotoGrid[i].PhotoRate = mycalloc(sizeof(*(Grid->PhotoGrid[i].PhotoRate)), Grid->NumCellsTotal); // Then allocate memory for each photoionization grid.
 
    // For some of the early snapshots we don't have photoionization grids (because ionization hasn't started yet at z=100).  
    // Let's put a flag to know whether we have any valid data for this snapshot so we don't have to create empty grids.
    if (i > LowSnap && i <= HighSnap)
    {
      Grid->PhotoGrid[i].valid_grid = 1;
    }
    else
    { 
      Grid->PhotoGrid[i].valid_grid = 0;
      printf("Snapshot %d is not a valid snapshot for reionization -- SKIPPING! --\n", i);
      continue;
    }
    snprintf(buf, MAXLEN, "%s/%s_%03d", PhotoionDir, PhotoionName, i);

    if(!(photoion = fopen(buf, "rb")))
    {
      fprintf(stderr, "Cannot open file %s\n", buf);
      return EXIT_FAILURE;
    }

    printf("Reading photoionization grid %s\n", buf);    
    fread(Grid->PhotoGrid[i].PhotoRate, sizeof(*(Grid->PhotoGrid[i].PhotoRate)), Grid->NumCellsTotal, photoion);
    fclose(photoion);

  }

  printf("All reionization feedback stuff read in successfully\n");

  return EXIT_SUCCESS;   
}

int32_t free_grid()
{

  int32_t i;

  for (i = 0; i < Grid->NumGrids; ++i)
  {
    myfree(Grid->PhotoGrid[i].PhotoRate, sizeof(*(Grid->PhotoGrid[i].PhotoRate)) * Grid->NumCellsTotal);
  }

  free(Grid->PhotoGrid);
  myfree(Grid->ReionRedshift, sizeof(*(Grid->ReionRedshift)));
  free(Grid);

  return EXIT_SUCCESS;

} 

double do_grid_reionization(int gal, double Zcurr, double *Mfilt)
{
  double z_reion;

  double M = 3.0e9; // Fits from Sobbachi 2015.
  double a = 0.17;
  double b = -2.1;
  double c = 2.0;
  double d = 2.3;
  double my_Mfilt, Mvir, PhotHI, reionization_modifier;

  int32_t grid_position, status;

  if (Grid->PhotoGrid[Gal[gal].SnapNum].valid_grid == 0) // In this instance reionization has not started yet so we don't have any data loaded.
  {
    return 1.0;
  }

  status= determine_1D_idx(Gal[gal].Pos[0], Gal[gal].Pos[1], Gal[gal].Pos[2], &grid_position); 
  if (status == EXIT_FAILURE)
  {
    ABORT(EXIT_FAILURE);
  }

  z_reion = Grid->ReionRedshift[grid_position]; // This is the redshift the cell was ionized at. 

  if(ReionizationOn == 2 && Zcurr < z_reion) // Has the cell been reionized yet? 
  {
    PhotHI = Grid->PhotoGrid[Gal[gal].SnapNum].PhotoRate[grid_position]/1.0e-12; // Photoionization Rate (in units of 1e-12).
    
    my_Mfilt = M * pow(PhotHI,a) * pow((1.0 + Zcurr)/10.0,b) * pow(1.0 - pow((1.0 + Zcurr)/(1.0 + z_reion), c), d);
    Mvir = Gal[gal].Mvir * 1.0e10 / Hubble_h;
 
    reionization_modifier = pow(2.0, -my_Mfilt/Mvir);
  }
  else
  {
    reionization_modifier = 1.0;
    my_Mfilt = 0.0; 
  }

  *Mfilt = my_Mfilt;
  return reionization_modifier;
  
}


// Local Functions //

