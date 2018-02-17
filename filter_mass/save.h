#ifndef SAVE_H 
#define SAVE_H 

#include <stdint.h>
#include "main.h"

// Structs //


// Proto-Types //

int32_t save_arrays(int64_t *HaloID, float *ReionMod, SAGE_params params, int32_t NHalos_Ionized, int32_t filenr, int32_t ThisSnap, int32_t first_run);
int32_t read_arrays(SAGE_params params, int32_t filenr, int32_t ThisSnap);
#endif
