#ifndef CHECK_INI_H 
#define CHECK_INI_H 

#ifdef RSAGE

#include "../grid-model/src/confObj.h"
int32_t check_ini_files(confObj_t *simParam, char *PhotoionName, char *PhotoionDir, int32_t fescPrescription, double alpha, double beta, double delta,
                        double quasar_baseline, double quasar_boosted, double N_dyntime, double MH_low, double fesc_low,
                        double MH_high, double fesc_high, int32_t HaloPartCut);

#endif

#endif
