#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars_grid.h"
#include "core_proto_grid.h"

// keep a static file handle to remove the need to do constant seeking.

int32_t save_grid(struct GRID_STRUCT *save_grid)
{
  char name_HI[MAXLEN], name_HeI[MAXLEN], name_HeII[MAXLEN];

  int32_t grid_num_idx;

  FILE* file_HI; 
  FILE* file_HeI;
  FILE* file_HeII;

  for (grid_num_idx = 0; grid_num_idx < Grid->NumGrids; ++grid_num_idx)  
  {    
 
    if (fescPrescription == 0)
    {
      snprintf(name_HI, MAXLEN, "%s/%s_fesc%.2f_HaloPartCut%d_nionHI_%03d", GridOutputDir, FileNameGalaxies, fesc, HaloPartCut, ListOutputGrid[grid_num_idx]);
      snprintf(name_HeI, MAXLEN, "%s/%s_fesc%.2f_HaloPartCut%d_nionHeI_%03d", GridOutputDir, FileNameGalaxies, fesc, HaloPartCut, ListOutputGrid[grid_num_idx]); 
      snprintf(name_HeII, MAXLEN, "%s/%s_fesc%.2f_HaloPartCut%d_nionHeII_%03d", GridOutputDir, FileNameGalaxies, fesc, HaloPartCut, ListOutputGrid[grid_num_idx]); 
    } else if (fescPrescription == 1)
    {
      snprintf(name_HI, MAXLEN, "%s/%s_fesc%.2f_KimmFiducial_nionHI_%03d", GridOutputDir, FileNameGalaxies, fesc, ListOutputGrid[grid_num_idx]); 
      snprintf(name_HeI, MAXLEN, "%s/%s_fesc%.2f_KimmFiducial_nionHeI_%03d", GridOutputDir, FileNameGalaxies, fesc, ListOutputGrid[grid_num_idx]); 
      snprintf(name_HeII, MAXLEN, "%s/%s_fesc%.2f_KimmFiducial_nionHeII_%03d", GridOutputDir, FileNameGalaxies, fesc, ListOutputGrid[grid_num_idx]); 
    } else if (fescPrescription == 2)
    {
      snprintf(name_HI, MAXLEN, "%s/%s_MH_%.3e_%.2f_%.3e_%.2f_HaloPartCut%d_nionHI_%03d", GridOutputDir, FileNameGalaxies, MH_low, fesc_low, MH_high, fesc_high, HaloPartCut, ListOutputGrid[grid_num_idx]); 
      snprintf(name_HeI, MAXLEN, "%s/%s_MH_%.3e_%.2f_%.3e_%.2f_HaloPartCut%d_nionHeI_%03d", GridOutputDir, FileNameGalaxies, MH_low, fesc_low, MH_high, fesc_high, HaloPartCut, ListOutputGrid[grid_num_idx]); 
      snprintf(name_HeII, MAXLEN, "%s/%s_MH_%.3e_%.2f_%.3e_%.2f_HaloPartCut%d_nionHeII_%03d", GridOutputDir, FileNameGalaxies, MH_low, fesc_low, MH_high, fesc_high, HaloPartCut, ListOutputGrid[grid_num_idx]); 
    } else if (fescPrescription == 3)
    { 
      snprintf(name_HI, MAXLEN, "%s/%s_Ejected_alpha%.3fbeta%.3f_HaloPartCut%d_nionHI_%03d", GridOutputDir, FileNameGalaxies, alpha, beta, HaloPartCut, ListOutputGrid[grid_num_idx]); 
      snprintf(name_HeI, MAXLEN, "%s/%s_Ejected_alpha%.3fbeta%.3f_HaloPartCut%d_nionHeI_%03d", GridOutputDir, FileNameGalaxies, alpha, beta, HaloPartCut, ListOutputGrid[grid_num_idx]); 
      snprintf(name_HeII, MAXLEN, "%s/%s_Ejected_alpha%.3fbeta%.3f_HaloPartCut%d_nionHeII_%03d", GridOutputDir, FileNameGalaxies, alpha, beta, HaloPartCut, ListOutputGrid[grid_num_idx]); 
    } 
    else if (fescPrescription == 4)
    {
      snprintf(name_HI, MAXLEN, "%s/%s_quasar_%.2f_%.2f_%.2f_HaloPartCut%d_nionHI_%03d", GridOutputDir, FileNameGalaxies, quasar_baseline, quasar_boosted, N_dyntime, HaloPartCut, ListOutputGrid[grid_num_idx]); 
      snprintf(name_HeI, MAXLEN, "%s/%s_quasar_%.2f_%.2f_%.2f_HaloPartCut%d_nionHeI_%03d", GridOutputDir, FileNameGalaxies, quasar_baseline, quasar_boosted, N_dyntime, HaloPartCut, ListOutputGrid[grid_num_idx]); 
      snprintf(name_HeII, MAXLEN, "%s/%s_quasar_%.2f_%.2f_%.2f_HaloPartCut%d_nionHeII_%03d", GridOutputDir, FileNameGalaxies, quasar_baseline, quasar_boosted, N_dyntime, HaloPartCut, ListOutputGrid[grid_num_idx]); 
    }
    else if (fescPrescription == 5 || fescPrescription == 6)
    {
      snprintf(name_HI, MAXLEN, "%s/%s_Anne_MH_%.3e_%.2f_%.3e_%.2f_HaloPartCut%d_nionHI_%03d", GridOutputDir, FileNameGalaxies, MH_low, fesc_low, MH_high, fesc_high, HaloPartCut, ListOutputGrid[grid_num_idx]); 
      snprintf(name_HeI, MAXLEN, "%s/%s_Anne_MH_%.3e_%.2f_%.3e_%.2f_HaloPartCut%d_nionHeI_%03d", GridOutputDir, FileNameGalaxies, MH_low, fesc_low, MH_high, fesc_high, HaloPartCut, ListOutputGrid[grid_num_idx]); 
      snprintf(name_HeII, MAXLEN, "%s/%s_Anne_MH_%.3e_%.2f_%.3e_%.2f_HaloPartCut%d_nionHeII_%03d", GridOutputDir, FileNameGalaxies, MH_low, fesc_low, MH_high, fesc_high, HaloPartCut, ListOutputGrid[grid_num_idx]);       
    }

    file_HI = fopen(name_HI, "w"); 
    if (file_HI == NULL)
    {
      fprintf(stderr, "Can't open file %s for writing\n", name_HI);
      return EXIT_FAILURE;
    }
    fwrite(save_grid->GridProperties[grid_num_idx].Nion_HI, sizeof(*(save_grid->GridProperties[grid_num_idx].Nion_HI)) * save_grid->NumCellsTotal, 1, file_HI);
    fclose(file_HI);

    printf("Saved to %s\n", name_HI);
    /*
    file_HeI = fopen(name_HeI, "w");
    if (file_HeI == NULL)
    {
      fprintf(stderr, "Can't open file %s for writing\n", name_HeI);
      return EXIT_FAILURE;
    }
    fwrite(save_grid->GridProperties[grid_num_idx].Nion_HeI, sizeof(*(save_grid->GridProperties[grid_num_idx].Nion_HeI)) * save_grid->NumCellsTotal, 1, file_HeI);
    fclose(file_HeI);

    file_HeII = fopen(name_HeII, "w");
    if (file_HeII == NULL)
    {
      fprintf(stderr, "Can't open file %s for writing\n", name_HeII);
      return EXIT_FAILURE;
    }
    fwrite(save_grid->GridProperties[grid_num_idx].Nion_HeII, sizeof(*(save_grid->GridProperties[grid_num_idx].Nion_HeII)) * save_grid->NumCellsTotal, 1, file_HeII);
    fclose(file_HeII);
    */   

  } // Snapshot loop.

  return EXIT_SUCCESS;    
}


int32_t save_fesc_properties(int filenr, int32_t merged_gal_flag) 
{

  FILE *save_properties;
  char name_properties[MAXLEN];
  int32_t i, j; 

  printf("Saving fesc properties\n");
  for (i = 0; i < Grid->NumGrids; ++i)
  {
    if (fescPrescription == 0)
    {
      snprintf(name_properties, MAXLEN, "%s/%s_fesc%.2f_HaloPartCut%d_fescproperties_%d.txt", GridOutputDir, FileNameGalaxies, fesc, HaloPartCut, i + LowSnap);
    }
    else if (fescPrescription == 3)
    {
      snprintf(name_properties, MAXLEN, "%s/%s_Ejected_alpha%.3fbeta%.3f_HaloPartCut%d_fescproperties_%d.txt", GridOutputDir, FileNameGalaxies, alpha, beta, HaloPartCut, i + LowSnap); 
    }
    else if (fescPrescription == 4)
    {
      snprintf(name_properties, MAXLEN, "%s/%s_quasar_%.2f_%.2f_%.2f_HaloPartCut%d_fescproperties_%d.txt", GridOutputDir, FileNameGalaxies, quasar_baseline, quasar_boosted, N_dyntime, HaloPartCut, i + LowSnap);
    }

    if (filenr == 0 && merged_gal_flag == 0) // For the first grid want to open up a fresh file
    {
      save_properties = fopen(name_properties, "w");
    }
    else
    {  
      save_properties = fopen(name_properties, "a"); // But for the others we want to append.
    }
    if (save_properties == NULL)
    {
      fprintf(stderr, "Cannot open file %s\n", name_properties);
      return EXIT_FAILURE;
    }

    for (j = 0; j < NtotGals; ++j)
    {
      if (Grid->GridProperties[i].SnapshotGalaxy[j] != -1)
      {
        fprintf(save_properties, "%d %.4f %.4f %.4f %.4e %.4e\n", Grid->GridProperties[i].SnapshotGalaxy[j], Grid->GridProperties[i].fescGalaxy[j], Grid->GridProperties[i].MvirGalaxy[j], Grid->GridProperties[i].MstarGalaxy[j], Grid->GridProperties[i].NgammaGalaxy[j], Grid->GridProperties[i].NgammafescGalaxy[j]); 
      }
    }
    fclose(save_properties);
  }
  return EXIT_SUCCESS;

}
/*
void save_redshift(void)
{

  char name_redshift[MAXLEN];
  int j;

  snprintf(name_redshift, MAXLEN, "%s/redshift_file.dat", GridOutputDir);

  printf("Writing out the saved redshifts to file %s\n", name_redshift);

  if (!(file_redshift = fopen(name_redshift,"w")))
  {
    printf("can't open file `%s'\n", name_redshift);
    ABORT(0);
  }

  for (j = 0; j < NGrid; ++j)
  {
    fprintf(file_redshift, "%.5f 1\n", ZZ[ListOutputGrid[NGrid - 1 - j]]);
  }

  fprintf(file_redshift, "6.000 0\n");
  fclose(file_redshift);

}
*/
