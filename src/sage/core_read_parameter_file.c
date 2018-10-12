#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"

#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300


// Local Proto-Types //

int32_t check_ini_parameters(void);

// Extern Functions // 

int32_t read_parameter_file(char *fname)
{
  FILE *fd;
  char buf[400], buf1[400], buf2[400], buf3[400];
  int i, j, nt = 0;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int errorFlag = 0; 

  strcpy(tag[nt], "FileNameGalaxies");
  addr[nt] = FileNameGalaxies;
  id[nt++] = STRING;

  strcpy(tag[nt], "OutputDir");
  addr[nt] = OutputDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "GridOutputDir");
  addr[nt] = GridOutputDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "FirstFile");
  addr[nt] = &FirstFile;
  id[nt++] = INT;

  strcpy(tag[nt], "LastFile");
  addr[nt] = &LastFile;
  id[nt++] = INT;

  strcpy(tag[nt], "TreeName");
  addr[nt] = TreeName;
  id[nt++] = STRING;

  strcpy(tag[nt], "TreeExtension");
  addr[nt] = TreeExtension;
  id[nt++] = STRING;

  strcpy(tag[nt], "SimulationDir");
  addr[nt] = SimulationDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileWithSnapList");
  addr[nt] = FileWithSnapList;
  id[nt++] = STRING;

  strcpy(tag[nt], "LastSnapShotNr");
  addr[nt] = &LastSnapShotNr;
  id[nt++] = INT;

  strcpy(tag[nt], "Omega");
  addr[nt] = &Omega;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "OmegaLambda");
  addr[nt] = &OmegaLambda;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "BaryonFrac");
  addr[nt] = &BaryonFrac;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Hubble_h");
  addr[nt] = &Hubble_h;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "PartMass");
  addr[nt] = &PartMass;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "BoxSize");
  addr[nt] = &BoxSize;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
  addr[nt] = &UnitVelocity_in_cm_per_s;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "UnitLength_in_cm");
  addr[nt] = &UnitLength_in_cm;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "UnitMass_in_g");
  addr[nt] = &UnitMass_in_g;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "self_consistent");
  addr[nt] = &self_consistent;
  id[nt++] = INT;

  strcpy(tag[nt], "ReionizationOn");
  addr[nt] = &ReionizationOn;
  id[nt++] = INT;

  strcpy(tag[nt], "SupernovaRecipeOn");
  addr[nt] = &SupernovaRecipeOn;
  id[nt++] = INT;

  strcpy(tag[nt], "DiskInstabilityOn");
  addr[nt] = &DiskInstabilityOn;
  id[nt++] = INT;

  strcpy(tag[nt], "SFprescription");
  addr[nt] = &SFprescription;
  id[nt++] = INT;

  strcpy(tag[nt], "AGNrecipeOn");
  addr[nt] = &AGNrecipeOn;
  id[nt++] = INT;

  strcpy(tag[nt], "QuasarRecipeOn");
  addr[nt] = &QuasarRecipeOn;
  id[nt++] = INT;

  strcpy(tag[nt], "IRA");
  addr[nt] = &IRA;
  id[nt++] = INT;

  strcpy(tag[nt], "TimeResolutionSN");
  addr[nt] = &TimeResolutionSN;
  id[nt++] = INT;

  strcpy(tag[nt], "RescaleSN");
  addr[nt] = &RescaleSN;
  id[nt++] = INT;

  strcpy(tag[nt], "IMF");
  addr[nt] = &IMF;
  id[nt++] = INT;

  strcpy(tag[nt], "SfrEfficiency");
  addr[nt] = &SfrEfficiency;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "FeedbackReheatingEpsilon");
  addr[nt] = &FeedbackReheatingEpsilon;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "FeedbackEjectionEfficiency");
  addr[nt] = &FeedbackEjectionEfficiency;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ReIncorporationFactor");
  addr[nt] = &ReIncorporationFactor;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "RadioModeEfficiency");
  addr[nt] = &RadioModeEfficiency;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "QuasarModeEfficiency");
  addr[nt] = &QuasarModeEfficiency;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "BlackHoleGrowthRate");
  addr[nt] = &BlackHoleGrowthRate;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ThreshMajorMerger");
  addr[nt] = &ThreshMajorMerger;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ThresholdSatDisruption");
  addr[nt] = &ThresholdSatDisruption;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Yield");
  addr[nt] = &Yield;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "RecycleFraction");
  addr[nt] = &RecycleFraction;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Reionization_z0");
  addr[nt] = &Reionization_z0;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Reionization_zr");
  addr[nt] = &Reionization_zr;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "EnergySN");
  addr[nt] = &EnergySN;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "LowSnap");
  addr[nt] = &LowSnap;
  id[nt++] = INT;

  strcpy(tag[nt], "GridSize");
  addr[nt] = &GridSize;
  id[nt++] = INT;

  strcpy(tag[nt], "PhotoionDir");
  addr[nt] = PhotoionDir;
  id[nt++] = STRING;

  strcpy(tag[nt], "PhotoionName");
  addr[nt] = PhotoionName;
  id[nt++] = STRING;

  strcpy(tag[nt], "ReionRedshiftName");
  addr[nt] = ReionRedshiftName; 
  id[nt++] = STRING;

  strcpy(tag[nt], "PhotonPrescription");
  addr[nt] = &PhotonPrescription;
  id[nt++] = INT;

  strcpy(tag[nt], "HaloPartCut");
  addr[nt] = &HaloPartCut;
  id[nt++] = INT;

  strcpy(tag[nt], "TimeResolutionStellar");
  addr[nt] = &TimeResolutionStellar;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "fescPrescription");
  addr[nt] = &fescPrescription;
  id[nt++] = INT;

  strcpy(tag[nt], "alpha");
  addr[nt] = &alpha;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "beta");
  addr[nt] = &beta;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "delta");
  addr[nt] = &delta;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "quasar_baseline");
  addr[nt] = &quasar_baseline;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "quasar_boosted");
  addr[nt] = &quasar_boosted;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "N_dyntime");
  addr[nt] = &N_dyntime;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "MH_low");
  addr[nt] = &MH_low;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "fesc_low");
  addr[nt] = &fesc_low;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "MH_high");
  addr[nt] = &MH_high;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "fesc_high");
  addr[nt] = &fesc_high;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "calcUVmag");
  addr[nt] = &calcUVmag;
  id[nt++] = INT;

  if((fd = fopen(fname, "r")))
  {
    while(!feof(fd))
    {
      *buf = 0;
      fgets(buf, 200, fd);
      if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
        continue;

      if(buf1[0] == '%' || buf1[0] == '-')
        continue;

      for(i = 0, j = -1; i < nt; i++)
        if(strcmp(buf1, tag[i]) == 0)
      {
        j = i;
        tag[i][0] = 0;
        break;
      }

      if(j >= 0)
      {
        switch (id[j])
        {
          case DOUBLE:
          *((double *) addr[j]) = atof(buf2);
          break;
          case STRING:
          strcpy(addr[j], buf2);
          break;
          case INT:
          *((int *) addr[j]) = atoi(buf2);
          break;
        }
      }
      else
      {
        //printf("Error in file %s:   Tag '%s' not allowed or multiple defined.\n", fname, buf1);        
      }
    }
    fclose(fd);

  }
  else
  {
    printf("Parameter file %s not found.\n", fname);
    errorFlag = 1;
  }

  for(i = 0; i < nt; i++)
  {
    if(*tag[i]) 
    {
      if (strcmp(tag[i], "TreeExtension") != 1)
      {
        memset(tag[i], 0, 50); 
      }
      else
      {
        printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
        errorFlag = 1;
      }
    }
  }
	
  if (errorFlag == 1)
  {
    return EXIT_FAILURE;
  }
	
	assert(LastSnapShotNr+1 > 0 && LastSnapShotNr+1 < ABSOLUTEMAXSNAPS);
	MAXSNAPS = LastSnapShotNr + 1;

  ListOutputSnaps[0] = LastSnapShotNr;

  int32_t status = check_ini_parameters();
  if (status != EXIT_SUCCESS)
  {
    return EXIT_FAILURE;
  }
 
  return EXIT_SUCCESS;
}

// Local Functions //

int32_t check_ini_parameters(void)
{
  if (GridSize > 1024)
  {
    fprintf(stderr, "The grid size cannot be greater than 1024.\n");
    return EXIT_FAILURE;
  }

  // First check that the recipe flags are set in acceptable ranges.

  switch (self_consistent)
  {
    case 0:
    case 1:
      break;

    default:
      fprintf(stderr, "The only valid options for `self_consistent` are 0 or 1. You entered %d\n", self_consistent);
      return EXIT_FAILURE;
  }

  switch (ReionizationOn)
  {
    case 0:
    case 1:
    case 2:
    case 3:
      break;

    default:
      fprintf(stderr, "The only valid options for `ReionizationOn` are 0, 1, 2 or 3. You entered %d\n", ReionizationOn);
      return EXIT_FAILURE;
  }

  switch (SupernovaRecipeOn)
  {
    case 0:
    case 1:
      break;

    default:
      fprintf(stderr, "The only valid options for `SupernovaRecipeOn` are 0 or 1. You entered %d\n", SupernovaRecipeOn);
      return EXIT_FAILURE;
  }

  switch (DiskInstabilityOn)
  {
    case 0:
    case 1:
      break;

    default:
      fprintf(stderr, "The only valid options for `DiskInstabilityOn` are 0 or 1. You entered %d\n", DiskInstabilityOn);
      return EXIT_FAILURE;
  }

  switch (SFprescription)
  {
    case 0:
      break;

    default:
      fprintf(stderr, "The only valid options for `SFprescription` are 0 or 1. You entered %d\n", SFprescription);
      return EXIT_FAILURE;
  }

  switch (AGNrecipeOn)
  {
    case 0:
    case 1:
    case 2:
    case 3:
      break;

    default:
      fprintf(stderr, "The only valid options for `AGNrecipeOn` are 0, 1, 2 or 3. You entered %d\n", AGNrecipeOn);
      return EXIT_FAILURE;
  }

  switch (QuasarRecipeOn)
  {
    case 1:
    case 2:
      break;

    default:
      fprintf(stderr, "The only valid options for `QuasarRecipeOn` are 1 or 2. You entered %d\n", QuasarRecipeOn);
      return EXIT_FAILURE;
  }

  switch (IRA)
  {
    case 0:
    case 1:
      break;

    default:
      fprintf(stderr, "The only valid options for `IRA` are 0 or 1. You entered %d\n", IRA);
      return EXIT_FAILURE;
  }

  switch (RescaleSN)
  {
    case 0:
    case 1:
      break;

    default:
      fprintf(stderr, "The only valid options for `RescaleSN` are 0 or 1. You entered %d\n", RescaleSN);
      return EXIT_FAILURE;
  }

  switch (IMF)
  {
    case 0:
    case 1:
      break;

    default:
      fprintf(stderr, "The only valid options for `IMF` are 0 or 1. You entered %d\n", IMF);
      return EXIT_FAILURE;
  }

  switch (fescPrescription)
  {
    case 0:
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
      break;

  default:
    fprintf(stderr, "The only valid options for `fescPrescription` are 0, 1, 2, 3, 4 or 5. You entered %d\n", fescPrescription);
    return EXIT_FAILURE;

  }

  // Then ensure that any cross dependancies are required.

  if (self_consistent == 0 && (ReionizationOn == 2 || ReionizationOn == 3))
  {
    fprintf(stderr, "You have turned off self consistent reionization (`self_consistent` = 0).  However your reionization prescription (`ReionizationOn`) relies on self consistent pipeline. Change either `self_consistent` or `ReionizationOn`.\n");
    return EXIT_FAILURE; 
  }

#ifndef RSAGE

  if (self_consistent == 1 || (ReionizationOn == 2 || ReionizationOn == 3))
  {
    fprintf(stderr, "You have requested self-consistent reionization to be on (`self_consistent = 1` or `ReionizationOn = 2 or 3`)." 
                    "  However you have not compiled with the self-consistent modules.\nEnsure that `BUILD_RSAGE` is set in the Makefile.\n" );
    return EXIT_FAILURE;
  }

#endif

  return EXIT_SUCCESS;
}
