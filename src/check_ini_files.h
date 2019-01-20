#ifndef CHECK_INI_H 
#define CHECK_INI_H 

#ifdef RSAGE

#include "../grid-model/src/confObj.h"

#endif

int32_t check_ini_files(char *FileNameGalaxies, char *OutputDir, char *GalaxyOutputDir,
                        char *GridOutputDir, char *PhotoionDir, char *PhotoionName,
                        char *ReionRedshiftName, char **argv, int32_t ThisTask);

#endif
