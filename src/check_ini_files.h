#ifndef CHECK_INI_H 
#define CHECK_INI_H 

#ifdef RSAGE

#include "../grid-model/src/confObj.h"

int32_t check_cifog_ini(confObj_t *simParam, char *OutputDir, char *RunPrefix,
                        char **argv, int32_t ThisTask);

#endif

int32_t check_sage_ini(char *RunPrefix, char *OutputDir, char *GalaxyOutputDir,
                       char *GridOutputDir, char *PhotoionDir, char *PhotoionName,
                       char *ReionRedshiftName, char **argv, int32_t ThisTask);

#endif
