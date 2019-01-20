#ifndef CHECK_INI_H 
#define CHECK_INI_H 

#ifdef RSAGE

#include "../grid-model/src/confObj.h"

int32_t check_cifog_ini(char *OutputDir, confObj_t *simParam, char **argv,
                        char *FileNameGalaxies, int32_t ThisTask);

#endif

int32_t check_sage_ini(char *FileNameGalaxies, char *OutputDir, char *GalaxyOutputDir,
                       char *GridOutputDir, char *PhotoionDir, char *PhotoionName,
                       char *ReionRedshiftName, char **argv, int32_t ThisTask);

#endif
