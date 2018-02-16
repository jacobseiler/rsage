#ifndef REIONIZATION_H 
#define REIONIZATION_H 

#include <stdint.h>

// Structs //
struct GRID_STRUCT
{
  int32_t GridSize;
  uint64_t NumCellsTotal;

  double *ReionRedshift; // This is the redshift the grid cell was ionized at.
  double *PhotoionizationRate;  
};

typedef struct GRID_STRUCT *grid_t;

// Proto-Types //

int32_t read_grid(int32_t snapnum, char *photoiondir, char *photoionname, grid_t *Grid); 

#endif
