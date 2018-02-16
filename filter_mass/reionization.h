#ifndef REIONIZATION_H 
#define REIONIZATION_H 

#include <stdint.h>
#include "tree_io.h"

// Structs //
struct GRID_STRUCT
{
  int32_t GridSize;
  uint64_t NumCellsTotal;
  double BoxSize;

  double *ReionRedshift; // This is the redshift the grid cell was ionized at.
  double *PhotoRate;  
};

typedef struct GRID_STRUCT *grid_t;

// Proto-Types //

int32_t read_grid(int32_t SnapNum, SAGE_params params, grid_t *Grid); 
int32_t free_grid(grid_t *Grid);
int32_t populate_halo_arrays(int32_t filenr, int32_t treenr, int32_t NHalos_ThisSnap, int32_t ThisSnap, halo_t Halos, grid_t Grid, int64_t **HaloID, float **Mfilt, int32_t *NHalos_Ionized);
int32_t determine_Mfilt(struct HALO_STRUCT Halo, grid_t Grid, float *filtering_mass);
int32_t determine_1D_idx(float pos_x, float pos_y, float pos_z, int32_t GridSize, double BoxSize, int32_t *grid_1D);
#endif
