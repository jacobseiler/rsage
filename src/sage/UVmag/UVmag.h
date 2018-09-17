#ifndef UVMAG_H 
#define UVMAG_H 

#include <stdint.h>

// Structs //

// Proto-Types //
int32_t init_UVlookup(void);
int32_t calc_LUV(struct GALAXY *g, float *LUV);
#endif
