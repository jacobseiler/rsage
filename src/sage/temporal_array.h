#ifndef COMMON_H 
#define COMMON_H 

#include "core_allvars.h"

// Proto-Types //

void update_temporal_array(int p, int halonr, int steps_completed);
void free_temporal_arrays(struct GALAXY *g);
int32_t malloc_temporal_arrays(struct GALAXY *g);
void write_temporal_arrays(struct GALAXY *g, FILE *fp);
#endif
