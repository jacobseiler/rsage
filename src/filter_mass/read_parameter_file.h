#ifndef PARAMETER_FILE_H
#define PARAMETER_FILE_H

#include <stdint.h>
#include "main.h"

// Proto-Types //

int32_t read_parameter_file(char *fname, struct SAGE_PARAMETERS *params); 
#endif
