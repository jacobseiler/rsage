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
#include "selfcon_grid.h"

// Local Variables //

// Local Proto-Types //

// External Functions //

int32_t init_selfcon_grid(void)
{

  int32_t status;

#ifdef MPI
  if (ThisTask == 0)
#endif
  printf("Initializing the grids for the self_consistent run.\n");

  if (fescPrescription == 0)
  {
#ifdef MPI
    if (ThisTask == 0)
#endif
    printf("\n\nUsing a constant escape fraction of %.4f\n", fesc); 
  }
  else if (fescPrescription == 1)
  {
    fprintf(stderr, "\n\nDeprecated! Use a different value for fescPrescription one!\n");
    return EXIT_FAILURE;
  } 
  else if (fescPrescription == 2)
  {
#ifdef MPI
    if (ThisTask == 0)
#endif
    printf("\n\nUsing an fesc prescription that scales with halo mass.\n\n");
    determine_fesc_constants();
  }
  else if (fescPrescription == 3)
  {
#ifdef MPI
    if (ThisTask == 0)
#endif
    printf("\n\nUsing an fesc prescription that scales with the fraction of ejected mass in the galaxy.\nThis takes the form A*fej + B with A = %.4e and B = %.4e\n", alpha, beta); 
  }
  else if (fescPrescription == 4)
  {
#ifdef MPI
    if (ThisTask == 0)
#endif
    {
      printf("\n\nUsing an fesc prescription that depends upon quasar activity.\n");
      printf("\nFor a galaxy that had a quasar event within %.2f dynamical times go the escape fraction will be %.2f.  Otherwise it will have a constant value of %.2f\n", N_dyntime, quasar_boosted, quasar_baseline);
    }
  }
  else if (fescPrescription == 5)
  {    
#ifdef MPI
    if (ThisTask == 0)
#endif
    printf("\n\nUsing Anne's functional form for an escape fraction that decreases for increasing halo mass.\n");
    XASSERT(fesc_low > fesc_high, "Input file contain fesc_low = %.2f and fesc_high = %.2f. For this prescription (fescPrescription == 5), we require fesc_low > fesc_high\n", fesc_low, fesc_high);
  }
  else if (fescPrescription == 6)
  {
#ifdef MPI
    if (ThisTask == 0)
#endif
    printf("\n\nUsing Anne's functional form for an escape fraction that increases for increasing halo mass.\n");
    XASSERT(fesc_low < fesc_high, "Input file contain fesc_low = %.2f and fesc_high = %.2f. For this prescription (fescPrescription == 6), we require fesc_low < fesc_high\n", fesc_low, fesc_high);
  }
  else
  {
    printf("\n\nOnly escape fraction prescriptions 0 to 6 (exlucding 1) are permitted.\n");
    return EXIT_FAILURE;
  }


  SelfConGrid = malloc(sizeof(*(SelfConGrid)));
  if (SelfConGrid == NULL)
  {
    fprintf(stderr, "Could not allocate memory for the high level self_consistent grid structure\n");
    return EXIT_FAILURE;
  }

  status = malloc_selfcon_grid(SelfConGrid);
  if (status != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

}

int32_t free_selfcon_grid(void)
{ 

  myfree(SelfConGrid->Nion_HI, sizeof(*(SelfConGrid->Nion_HI)) * Grid->NumCellsTotal);
  myfree(SelfConGrid->GalCount, sizeof(*(SelfConGrid->GalCount)) * Grid->NumCellsTotal);
  
  free(Grid);

  return EXIT_SUCCESS;

}


// Local Functions //

void determine_fesc_constants(void)
{ 
  
  double A, B, log_A;
  
  if (fescPrescription == 2)
  { 
    log_A = (log10(fesc_high) - (log10(fesc_low)*log10(MH_high)/log10(MH_low))) * pow(1 - (log10(MH_high) / log10(MH_low)), -1);
    B = (log10(fesc_low) - log_A) / log10(MH_low);
    A = pow(10, log_A);
    
    alpha = A;
    beta = B;
    
    printf("Fixing the points (%.4e, %.2f) and (%.4e, %.2f)\n", MH_low, fesc_low, MH_high, fesc_high);
    printf("This gives a power law with constants A = %.4e, B = %.4e\n", alpha, beta);
  }
 
}

int32_t malloc_selfcon_grid(struct SELFCON_GRID_STRUCT *grid)
{

#define ALLOCATE_GRID_MEMORY(name, length) \
{                                          \
  name = mycalloc(length, sizeof(*(name)));  \
  if (name == NULL)                        \
  {                                        \
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate"#name".\n", sizeof(*(name)* length)); \
    return EXIT_FAILURE;                   \
  }                                        \
}

  uint64_t cell_idx;

  if (grid == NULL)
  {
    fprintf(stderr, "`init_selfcon_grid` was called with a SELFCON_GRID_STRUCT pointer that has not been initialized\n");
    return EXIT_FAILURE;
  }

  grid->GridSize = GridSize;
  grid->NumCellsTotal = CUBE(GridSize);

  ALLOCATE_GRID_MEMORY(grid->Nion_HI, grid->NumCellsTotal);
  ALLOCATE_GRID_MEMORY(grid->GalCount, grid->NumCellsTotal);

 for (cell_idx = 0; cell_idx < Grid->NumCellsTotal; ++cell_idx)
 { 
      grid->Nion_HI[cell_idx] = 0.0;
      grid->GalCount[cell_idx] = 0;
 }

  return EXIT_SUCCESS;

#undef ALLOCATE_GRID_MEMORY

}
