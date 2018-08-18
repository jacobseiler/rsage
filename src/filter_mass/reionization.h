#ifndef REIONIZATION_H 
#define REIONIZATION_H 

#include <stdint.h>
#include "tree_io.h"
#include "filter_mass.h"

// Structs //
struct GRID_STRUCT
{
  int32_t GridSize;
  uint64_t NumCellsTotal;
  double BoxSize;

  double *ReionRedshift; // This is the redshift the grid cell was ionized at.
  double *PhotoRate;  
} *Grid;

// Proto-Types //

int32_t filter_mass_read_grid(int32_t SnapNum, int32_t first_update_flag, int32_t GridSize, double BoxSize,
                              char *PhotoionDir, char *ReionRedshiftName, char *PhotoionName,
                              struct GRID_STRUCT *Grid);
int32_t filter_mass_free_grid(struct GRID_STRUCT *Grid);
int32_t populate_halo_arrays(int32_t filenr, int32_t treenr, int32_t NHalos_ThisTree, int32_t ThisSnap, int32_t first_update_flag, struct HALO_STRUCT *Halos, struct GRID_STRUCT *Grid, struct SAGE_PARAMETERS *params, int64_t **HaloID, float **ReionMod, int32_t *NHalos_ThisSnap, int32_t *NHalos_Ionized, int32_t *NHalos_In_Region, float *sum_ReionMod);
int32_t determine_Mfilt(struct HALO_STRUCT Halo, struct GRID_STRUCT *Grid, struct SAGE_PARAMETERS *params, float *filtering_mass, int32_t *NHalos_In_Regions);
int32_t filter_mass_determine_1D_idx(float pos_x, float pos_y, float pos_z, int32_t GridSize, double BoxSize, int32_t *grid_1D);
#endif
