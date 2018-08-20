#ifndef FILTER_MASS_H 
#define FILTER_MASS_H 

#include <stdint.h>

// Structs //

struct SAGE_PARAMETERS
{

  char *FileNameGalaxies; 
  char *TreeDir;
  char *TreeName;
  char *PhotoionDir;
  char *PhotoionName;
  char *ReionRedshiftName;

  int32_t GridSize;
  double BoxSize;

  double Hubble_h;

  int32_t FirstFile;
  int32_t LastFile;

  double Redshift;
 
} *SAGE_params;

// Proto-Types //

int32_t filter_masses(char *FileNameGalaxies, char *TreeDir, char *TreeName, 
                      char *PhotoionDir, char *PhotoionName, char *ReionRedshiftName,
                      int32_t FirstFile, int32_t LastFile, int32_t GridSize, double BoxSize,
                      double Hubble_h, int32_t SnapNum, double Redshift, int32_t first_update_flag,
                      int32_t ThisTask, int32_t NTask);
#endif
