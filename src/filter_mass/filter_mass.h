#ifndef MAIN_H 
#define MAIN_H 

#include <stdint.h>

#define ABSOLUTEMAXSNAPS 999

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

#endif
