#ifndef MAIN_H 
#define MAIN_H 

#include <stdint.h>

#define ABSOLUTEMAXSNAPS 999

// Structs //
struct SAGE_PARAMETERS
{

  char *TreeDir;
  char *TreeName;
  char *PhotoionDir;
  char *PhotoionName;
  char *ReionRedshiftName;
  char *SnapListFile;

  int32_t GridSize;
  double BoxSize;

  double Hubble_h;

  int32_t LastSnapshotNr;
  int32_t FirstFile;
  int32_t LastFile;

  double ZZ[ABSOLUTEMAXSNAPS];
 
};
typedef struct SAGE_PARAMETERS *SAGE_params;

#endif
