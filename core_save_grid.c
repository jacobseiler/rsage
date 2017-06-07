#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars_grid.h"
#include "core_proto_grid.h"

// keep a static file handle to remove the need to do constant seeking.
FILE* file_HI = NULL;
FILE* file_HeI = NULL;
FILE* file_HeII = NULL;
FILE* save_density = NULL;
FILE* file_redshift = NULL;
FILE* file_count = NULL;

void save_grid(int GridNr)
{
  char name_HI[MAXLEN], name_HeI[MAXLEN], name_HeII[MAXLEN], tag[MAXLEN];
  int i;

  if ((NGrid -1 - GridNr) < 10)
    snprintf(tag, MAXLEN, "%s", "00");
  else if ((NGrid - 1 - GridNr) < 100 && (NGrid - 1 - GridNr) >= 10)
    snprintf(tag, MAXLEN, "%s", "0");
  else
    snprintf(tag, MAXLEN, "%s", "");

  if (PhotonPrescription == 0) // Changes the output tag depending on whether we are using halo or galaxy photon prescription.
          snprintf(name_HI, MAXLEN, "%s/Halos_fgamma_%.0f_%s_nionHI_%s%d", GridOutputDir, SourceEfficiency, FileNameGalaxies, tag, (NGrid-1) - GridNr); 
  else
  {
    if (fescPrescription == 0)
    {
      snprintf(name_HI, MAXLEN, "%s/Galaxies_%s_fesc%.2f_HaloPartCut%d_nionHI_%s%d", GridOutputDir, FileNameGalaxies, fesc, HaloPartCut, tag, (NGrid - 1) - GridNr); 
      snprintf(name_HeI, MAXLEN, "%s/Galaxies_%s_fesc%.2f_nionHeI_%s%d", GridOutputDir, FileNameGalaxies, fesc, tag, (NGrid - 1) - GridNr); 
      snprintf(name_HeII, MAXLEN, "%s/Galaxies_%s_fesc%.2f_nionHeII_%s%d", GridOutputDir, FileNameGalaxies, fesc, tag, (NGrid - 1) - GridNr); 
    } else if (fescPrescription == 1)
    {
      snprintf(name_HI, MAXLEN, "%s/Galaxies_%s_fesc%.2f_KimmFiducial_nionHI_%s%d", GridOutputDir, FileNameGalaxies, fesc, tag, (NGrid - 1) - GridNr); 
      snprintf(name_HeI, MAXLEN, "%s/Galaxies_%s_fesc%.2f_KimmFiducial_nionHeI_%s%d", GridOutputDir, FileNameGalaxies, fesc, tag, (NGrid - 1) - GridNr); 
      snprintf(name_HeII, MAXLEN, "%s/Galaxies_%s_fesc%.2f_KimmFiducial_nionHeII_%s%d", GridOutputDir, FileNameGalaxies, fesc, tag, (NGrid - 1) - GridNr); 
    } else if (fescPrescription == 2)
    {
      snprintf(name_HI, MAXLEN, "%s/Galaxies_%s_MH_alpha%.2fbeta%.2f_nionHI_%s%d", GridOutputDir, FileNameGalaxies, log10(alpha), beta, tag, (NGrid - 1) - GridNr); 
      snprintf(name_HeI, MAXLEN, "%s/Galaxies_%s_MH_alpha%.2fbeta%.2f_nionHeI_%s%d", GridOutputDir, FileNameGalaxies, log10(alpha), beta, tag, (NGrid - 1) - GridNr); 
      snprintf(name_HeII, MAXLEN, "%s/Galaxies_%s_MH_alpha%.2fbeta%.2f_nionHeII_%s%d", GridOutputDir, FileNameGalaxies, log10(alpha), beta, tag, (NGrid - 1) - GridNr); 
    } else if (fescPrescription == 3)
    { 
      snprintf(name_HI, MAXLEN, "%s/Galaxies_%s_Ejected_alpha%.3fbeta%.3f_nionHI_%s%d", GridOutputDir, FileNameGalaxies, alpha, beta, tag, (NGrid - 1) - GridNr); 
      snprintf(name_HeI, MAXLEN, "%s/Galaxies_%s_Ejected_alpha%.3fbeta%.3f_nionHeI_%s%d", GridOutputDir, FileNameGalaxies, alpha, beta, tag, (NGrid - 1) - GridNr); 
      snprintf(name_HeII, MAXLEN, "%s/Galaxies_%s_Ejected_alpha%.3fbeta%.3f_nionHeII_%s%d", GridOutputDir, FileNameGalaxies, alpha, beta, tag, (NGrid - 1) - GridNr); 
    } else if (fescPrescription == 4)
    {
      snprintf(name_HI, MAXLEN, "%s/Galaxies_%s_supersteepslope_nionHI_%s%d", GridOutputDir, FileNameGalaxies, tag, (NGrid - 1) - GridNr); 
      //snprintf(name_HeI, MAXLEN, "%s/Galaxies_%s_slope_nionHeI_%s%d", GridOutputDir, FileNameGalaxies, tag, (NGrid - 1) - GridNr); 
      //snprintf(name_HeII, MAXLEN, "%s/Galaxies_%s_slope_nionHeII_%s%d", GridOutputDir, FileNameGalaxies, tag, (NGrid - 1) - GridNr); 
    }

    // Only calculate Helium ionizing photons if we aren't using Halo-prescription.
    /*
    if (!(file_HeI = fopen(name_HeI,"w")))
    {
      printf("can't open file `%s'\n", name_HeI);
      exit(0);
    }

    if (!(file_HeII = fopen(name_HeII,"w")))
    {
      printf("can't open file `%s'\n", name_HeII);
      exit(0);
    }
    */
  }

  if (!(file_HI = fopen(name_HI,"w")))
  {
    printf("can't open file `%s'\n", name_HI);
    exit(0);
  }

  /*
  snprintf(name_count, MAXLEN, "%s/Galaxies_%s_fesc%.2f_count_%s%d", GridOutputDir, FileNameGalaxies, fesc, tag, (NGrid - 1) - GridNr); 
  if (!(file_count = fopen(name_count,"w")))
  {
    printf("can't open file `%s'\n", name_count);
    exit(0);
  }
  */
  for (i = 0; i < CUBE(GridSize); ++i)
  {
    fwrite(&Grid[i].Nion_HI, sizeof(double), 1, file_HI);
    if (PhotonPrescription != 0)
    {
//      fwrite(&Grid[i].Nion_HeI, sizeof(double), 1, file_HeI);
//      fwrite(&Grid[i].Nion_HeII, sizeof(double), 1, file_HeII); 
    }
//    fwrite(&Grid[i].Count, sizeof(int), 1, file_count);
    //fwrite(&Grid[i].Density, sizeof(double), 1, save_density);
  }

  printf("Saved Snapshot %d to file %s.\n", ListOutputGrid[GridNr], name_HI);
  fclose(file_HI);
 
  if (PhotonPrescription != 0)
  {
//    fclose(file_HeI);
//    fclose(file_HeII);
    file_HeI = NULL;
    file_HeII = NULL;
  }

  file_HI = NULL;
  save_density = NULL;
}

void save_redshift(void)
{

  char name_redshift[MAXLEN];
  int j;

  snprintf(name_redshift, MAXLEN, "%s/redshift_file.dat", GridOutputDir);

  printf("Writing out the saved redshifts to file %s\n", name_redshift);

  if (!(file_redshift = fopen(name_redshift,"w")))
  {
    printf("can't open file `%s'\n", name_redshift);
    exit(0);
  }

  for (j = 0; j < NGrid; ++j)
  {
    fprintf(file_redshift, "%.5f 1\n", ZZ[ListOutputGrid[NGrid - 1 - j]]);
  }

  fprintf(file_redshift, "6.000 0\n");
  fclose(file_redshift);

}
