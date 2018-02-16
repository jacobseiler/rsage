#ifndef MAIN_H 
#define MAIN_H 

#include <stdint.h>

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

  int32_t FirstFile;
  int32_t LastFile;
 
};
typedef struct SAGE_PARAMETERS *SAGE_params;

#endif
