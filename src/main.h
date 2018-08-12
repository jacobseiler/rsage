#ifndef MAIN_H 
#define MAIN_H 

#define	CUBE(x) (x*x*x)
#define MAX_STRING_LEN 1024

#include <stdint.h>

// Structs //

// Proto-Types //
int32_t update_reion_redshift(int32_t SnapNum, double redshift, int32_t GridSize, int32_t first_update_flag,
                              char *PhotionDir, char *FileNameGalaxies, char *ReionredshiftName); 
#endif
