#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>

#define reionization_threshold 0.9

int32_t update_reion_redshift(int32_t SnapNum, double redshift, int32_t GridSize, int32_t first_update_flag,
                              char *PhotionDir, char *FileNameGalaxies, char *ReionredshiftName) 
{

  int32_t grid_idx, num_read, num_write;

  double *reion_grid;
  reion_grid = malloc(sizeof(*(*reion_grid)) * CUBE(GridSize)); 
  if (reion_grid == NULL)
  {
    fprintf(stderr, "Could not allocate memory for the reionization redshift grid at snapshot %d.\n", SnapNum);
    return EXIT_FAILURE;
  }

  char *reion_reion_fname[MAX_STRING_LEN];
  snprintf(reion_fname, MAX_STRING_LEN - 1, "%s/%s", PhotionDir, ReionredshiftName);

  FILE *reion_grid_file;

  // If this is the first iteration of RSAGE we first need to initialize the reionization redshift
  // grid. Otherwise we read in the previous values. 
  if (first_update_flag == 1)
  {
    for (grid_idx = 0; grid_idx < CUBE(GridSize); ++grid_idx)
    {
      reion_grid[grid_idx] = -1.0;
    }
  }
  else
  {
    reion_grid_file = fopen(reion_fname, "rb");
    if (reion_grid_file == NULL)
    {
      fprintf(stderr, "Could not open file %s for reading.\n", reion_fname);
      return EXIT_FAILURE;
    }
    
    num_read = fread(reion_grid, sizeof(*(*reion_grid)), CUBE(GridSize), reion_grid_file);
    if (num_read != CUBE(GridSize))
    {
      fprintf(stderr, "Attempted to read %d elements but only read %d elements from file %s\n", 
              CUBE(GridSize), num_read, reion_fname);
      return EXIT_FAILURE;
    }  

    fclose(reion_grid_file);
  }

  // Now that the grid is initialized, read the current XHII grid.
  double *XHII_grid;
  XHII_grid = malloc(sizeof(*(*XHII_grid)) * CUBE(GridSize)); 
  if (XHII_grid == NULL)
  {
    fprintf(stderr, "Could not allocate memory for reading the XHII grid of snapshot %d\n", SnapNum); 
    return EXIT_FAILURE;
  }

  char *XHII_fname[MAX_STRING_LEN];
  snprintf(XHII_fname, MAX_STRING_LEN - 1, "%s/%s_%03d", PhotionDir, FileNameGalaxies, SnapNum);

  FILE *XHII_file;
  XHII_file = fopen(XHII_fname, "rb");
  if (XHII_file == NULL)
  {
    fprintf(stderr, "Could not open file %s\n", XHII_fname);
    return EXIT_FAILURE;
  }

  num_read = fread(XHII_grid, sizeof(*(*XHII_grid)), CUBE(GridSize), XHII_grid_file);  
  if (num_read != CUBE(GridSize))
  {
    fprintf(stderr, "Attempted to read %d elements but only read %d elements from file %s\n", 
            CUBE(GridSize), num_read, XHII_fname);
    return EXIT_FAILURE;
  }  

  // Now go through the XHII grid and update the reionization redshift grid.
  int32_t count = 0;
  for (grid_idx = 0; grid_idx < CUBE(GridSize); ++grid_idx)
  {
    if (reion_grid[grid_idx] < 0.0) && (XHII_grid[grid_idx] > reionization_threshold)
    {
      reion_grid[grid_idx] = redshift;
      ++count;
    }
  }

  printf("For Snapshot %d there were %d cells ionized.\n", SnapNum, count);

  // Time to save the updated reionization redshift.
  reion_grid_file = fopen(reion_fname, "wb"); 
  if (reion_grid_file == NULL)
  {
    fprintf(stderr, "Could not open file %s for writing.\n", reion_fname);
    return EXIT_FAILURE;
  }

  num_write = fwrite(reion_grid, sizeof(*(reion_grid)), CUBE(GridSize), reion_grid_file);
  if (num_write != CUBE(GridSize))
  {
    fprintf(stderr, "Attempted to write %d elements but only wrote %d elements from file %s\n", 
            CUBE(GridSize), num_write, reion_fname);
    return EXIT_FAILURE;
  }

  fclose(reion_grid_file);

  free(XHII_grid);
  free(reion_grid);

  return EXIT_SUCCESS;
}
