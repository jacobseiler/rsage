#ifndef PARAMETER_FILE_H
#define PARAMETER_FILE_H

#include <stdint.h>

// Proto-Types //

int32_t read_parameter_file(char *fname, char **treedir, char **treename, char **photoiondir, char **photoionname, char **reionredshiftname, int32_t *numfiles);
#endif
