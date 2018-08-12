#ifndef REION_REDSHIFT_H 
#define REION_REDSHIFT_H 

#include <stdint.h>
#include "main.h"

// Structs //

// Proto-Types //
int32_t update_reion_redshift(int32_t SnapNum, double redshift, int32_t GridSize, int32_t first_update_flag,
                              char *PhotionDir, char *FileNameGalaxies, char *ReionredshiftName); 
#endif
