#ifndef MAIN_H 
#define MAIN_H 

#include <stdint.h>

#define ABSOLUTEMAXSNAPS 999

#define ABORT(sigterm)                                                  \
do {                                                                \
  printf("Error in file: %s\tfunc: %s\tline: %i\n", __FILE__, __FUNCTION__, __LINE__); \
  myexit(sigterm);                                                \
} while(0)

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
