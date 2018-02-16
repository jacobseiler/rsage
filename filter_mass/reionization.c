#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>

#include "reionization.h"

#define	CUBE(x) (x*x*x)
#define MAXLEN 1024

int32_t read_grid(int32_t SnapNum, SAGE_params params, grid_t *Grid)
{

  FILE *ReionRedshiftFile, *PhotoionFile;
  char RedshiftGridName[MAXLEN], PhotoionGridName[MAXLEN]; 

  *Grid = malloc(sizeof(struct GRID_STRUCT));
  if (*Grid == NULL)
  {
    fprintf(stderr, "Could not allocate memory for Grid struct.\n");
    return EXIT_FAILURE;
  }

  (*Grid)->GridSize = params->GridSize;
  (*Grid)->NumCellsTotal = CUBE(params->GridSize);
  (*Grid)->BoxSize = params->BoxSize;


  (*Grid)->ReionRedshift = malloc(sizeof(*((*Grid)->ReionRedshift)) * (*Grid)->NumCellsTotal);   
  if ((*Grid)->ReionRedshift == NULL)
  {
    fprintf(stderr, "Coult not allocate memory for the reionization redshift grid.\n");
    return EXIT_FAILURE;
  }

  (*Grid)->PhotoRate = malloc(sizeof(*((*Grid)->PhotoRate)) * (*Grid)->NumCellsTotal);   
  if ((*Grid)->PhotoRate == NULL)
  {
    fprintf(stderr, "Coult not allocate memory for the PhotoionFileization rate grid.\n");
    return EXIT_FAILURE;
  }

  snprintf(RedshiftGridName, MAXLEN, "%s/%s", params->PhotoionDir, params->ReionRedshiftName); 
  if (!(ReionRedshiftFile = fopen(RedshiftGridName, "rb")))
  { 
    fprintf(stderr, "Could not open file %s\n", RedshiftGridName);
    return EXIT_FAILURE;
  }

  printf("Reading reionization redshift file.\n");
  fread((*Grid)->ReionRedshift, sizeof(*((*Grid)->ReionRedshift)), (*Grid)->NumCellsTotal, ReionRedshiftFile);
  fclose(ReionRedshiftFile);

  snprintf(PhotoionGridName, MAXLEN, "%s/%s_%03d", params->PhotoionDir, params->PhotoionName, SnapNum); 
  if (!(PhotoionFile = fopen(PhotoionGridName, "rb")))
  { 
    fprintf(stderr, "Could not open file %s\n", PhotoionGridName);
    return EXIT_FAILURE;
  }

  printf("Reading photoionization rate.\n"); 
  fread((*Grid)->PhotoRate, sizeof(*((*Grid)->PhotoRate)), (*Grid)->NumCellsTotal, PhotoionFile);
  fclose(PhotoionFile);

  return EXIT_SUCCESS;

}

int32_t free_grid(grid_t *Grid)
{

  free((*Grid)->PhotoRate);
  free((*Grid)->ReionRedshift);

  free(*Grid);

  return EXIT_SUCCESS;

}
 
int32_t populate_halo_arrays(int32_t filenr, int32_t treenr, int32_t NHalos_ThisSnap, int32_t ThisSnap, halo_t Halos, grid_t Grid, int64_t **HaloID, float **Mfilt, int32_t *NHalos_Ionized)
{

  int32_t halonr, status;
  float Mfilt_tmp;

  for (halonr = 0; halonr < NHalos_ThisSnap; ++halonr)
  {

    if (Halos[halonr].SnapNum == ThisSnap)
    {
      status = determine_Mfilt(Halos[halonr], Grid, &Mfilt_tmp);
      if (status == EXIT_FAILURE)
      {
        return EXIT_FAILURE;
      }
    }
  } 

  return EXIT_SUCCESS;

}

int32_t determine_Mfilt(struct HALO_STRUCT Halo, grid_t Grid, float *filtering_mass)
{ 

  int32_t status, grid_idx;

  status = determine_1D_idx(Halo.Pos[0], Halo.Pos[1], Halo.Pos[2], Grid->GridSize, Grid->BoxSize, &grid_idx);
  if (status == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  } 

  return EXIT_SUCCESS;

}

int32_t determine_1D_idx(float pos_x, float pos_y, float pos_z, int32_t GridSize, double BoxSize, int32_t *grid_1D)
{

  int32_t x_grid, y_grid, z_grid;

  x_grid = pos_x * GridSize/BoxSize;
  y_grid = pos_y * GridSize/BoxSize;
  z_grid = pos_z * GridSize/BoxSize;

  *grid_1D = (z_grid*GridSize+y_grid)*GridSize+x_grid; // Convert the grid (x,y,z) to a 1D value.

  if(*grid_1D > CUBE(GridSize) || *grid_1D < 0) // Sanity check to ensure that no Grid Positions are outside the box.
  {
    fprintf(stderr, "Found a Grid Position outside the bounds of the box or negative; grid_position = %d\nPos[0] = %.4f\t Pos[1] = %.4f\tPos[2] = %.4f", *grid_1D, pos_x, pos_y, pos_z);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

