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

// Local Proto-Type //

int32_t count_ionized_cells(grid_t Grid);

int32_t read_grid(int32_t SnapNum, int32_t first_run, SAGE_params params, grid_t *Grid)
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
  count_ionized_cells(*Grid);
  fclose(ReionRedshiftFile);
  
  snprintf(PhotoionGridName, MAXLEN, "%s/%s_%03d", params->PhotoionDir, params->PhotoionName, SnapNum + 1); 
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

int32_t count_ionized_cells(grid_t Grid)
{

  int32_t cell_idx, num_ionized_cells = 0;

  for (cell_idx = 0; cell_idx < Grid->NumCellsTotal; ++cell_idx)
  {
  
    if (Grid->ReionRedshift[cell_idx] > -1.0)
    {
      ++num_ionized_cells;
    } 

  }

  printf("There were %d Cells marked as ionized from the reionization redshift grid.\n", num_ionized_cells);

  return EXIT_SUCCESS;

}

int32_t free_grid(grid_t *Grid)
{

  free((*Grid)->PhotoRate);
  free((*Grid)->ReionRedshift);

  free(*Grid);

  return EXIT_SUCCESS;

}
 
int32_t populate_halo_arrays(int32_t filenr, int32_t treenr, int32_t NHalos_ThisTree, int32_t ThisSnap, int32_t first_run, halo_t Halos, grid_t Grid, SAGE_params params, int64_t **HaloID, float **ReionMod, int32_t *NHalos_ThisSnap, int32_t *NHalos_Ionized, int32_t *NHalos_In_Regions, float *sum_ReionMod)
{

  int32_t halonr, status;
  float ReionMod_tmp;
  int64_t unique_ID;

  if (first_run == 1)
  {
    return EXIT_SUCCESS;
  }

  for (halonr = 0; halonr < NHalos_ThisTree; ++halonr)
  {
    if (Halos[halonr].SnapNum == ThisSnap) // Only care about halos at the snapshot specified.
    {
      ++(*NHalos_ThisSnap);

      status = determine_Mfilt(Halos[halonr], Grid, params, &ReionMod_tmp, NHalos_In_Regions); // Determine the reionization modifier for this halo.
      if (status == EXIT_FAILURE)
      {
        return EXIT_FAILURE;
      }
  
      if (ReionMod_tmp > 0.9999) // Only want to save those halos with a reionization modifier (appreciably) below 1.
      {
        continue;
      }

      unique_ID = (int64_t) treenr << 32 | halonr; // We create a unique ID for each halo within the file by generating a 64 bit number with the left-most 32 bits being the tree number and the right-most bits being the halo number.     
                                                   // As the tree number can only increase, this creates an ascending list without the need to sort. 
      printf("Unique ID for tree number %d and halo number %d is %ld with ReionMod %.4f\n", treenr, halonr, (long)unique_ID, ReionMod_tmp); 
      (*HaloID)[(*NHalos_Ionized)] = unique_ID;
      (*ReionMod)[(*NHalos_Ionized)] = ReionMod_tmp;
      (*sum_ReionMod) += ReionMod_tmp;

      ++(*NHalos_Ionized);


    }
  } 

  return EXIT_SUCCESS;

}

int32_t determine_Mfilt(struct HALO_STRUCT Halo, grid_t Grid, SAGE_params params, float *reionization_modifier, int32_t *NHalos_In_Regions) 
{ 

  int32_t status, grid_idx;

  double z_reion;

  double M = 3.0e9; // Fits from Sobbachi 2015.
  double a = 0.17;
  double b = -2.1;
  double c = 2.0;
  double d = 2.3;
  double Mfilt, Mvir, PhotHI;
  double Zcurr = params->ZZ[Halo.SnapNum]; 
  
  status = determine_1D_idx(Halo.Pos[0], Halo.Pos[1], Halo.Pos[2], Grid->GridSize, Grid->BoxSize, &grid_idx);
  if (status == EXIT_FAILURE)
  {
    return EXIT_FAILURE;
  } 

  z_reion = Grid->ReionRedshift[grid_idx]; // This is the redshift the cell was ionized at. 

  if(Zcurr < z_reion) // Has the cell been reionized yet? 
  {
    ++(*NHalos_In_Regions);
    PhotHI = Grid->PhotoRate[grid_idx]/1.0e-12; // Photoionization Rate (in units of 1e-12).

    Mfilt = M * pow(PhotHI,a) * pow((1.0 + Zcurr)/10.0,b) * pow(1.0 - pow((1.0 + Zcurr)/(1.0 + z_reion), c), d);
    Mvir = Halo.Mvir * 1.0e10 / params->Hubble_h;

    *reionization_modifier = pow(2.0, -Mfilt/Mvir);
  }
  else
  {
    *reionization_modifier = 1.0;
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

