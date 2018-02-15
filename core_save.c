#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"

#define TREE_MUL_FAC        (1000000000LL)
#define FILENR_MUL_FAC      (1000000000000000LL)

// keep a static file handle to remove the need to do constant seeking.
FILE* save_fd[ABSOLUTEMAXSNAPS] = { 0 };
FILE* save_fd2 = NULL;

void save_galaxies(int filenr, int tree)
{
  char buf[1000];
  int i, n;
  //struct GALAXY_OUTPUT galaxy_output;

  int OutputGalCount[MAXSNAPS], *OutputGalOrder;
  
  OutputGalOrder = (int*)malloc( NumGals*sizeof(int) );
  assert( OutputGalOrder );

  // reset the output galaxy count and order
  for(i = 0; i < MAXSNAPS; i++)
    OutputGalCount[i] = 0;
  for(i = 0; i < NumGals; i++)
    OutputGalOrder[i] = -1;

  // first update mergeIntoID to point to the correct galaxy in the output
  for(n = 0; n < NOUT; n++)
  {
    for(i = 0; i < NumGals; i++)
    {
      if(HaloGal[i].SnapNum == ListOutputSnaps[n])
      {
        OutputGalOrder[i] = OutputGalCount[n];
        OutputGalCount[n]++;
      }
    }
  }
    
  for(i = 0; i < NumGals; i++)
    if(HaloGal[i].mergeIntoID > -1)
      HaloGal[i].mergeIntoID = OutputGalOrder[HaloGal[i].mergeIntoID];    
  
  // now prepare and write galaxies
  for(n = 0; n < NOUT; n++)
  {

    
    // only open the file if it is not already open.
    if( !save_fd[n] )
    {
      sprintf(buf, "%s/%s_z%1.3f_%d", OutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[n]], filenr);

      if(!(save_fd[n] = fopen(buf, "r+")))
      {
				printf("can't open file `%s'\n", buf);
				ABORT(0);
      }

			// write out placeholders for the header data.
			size_t size = (Ntrees + 2)*sizeof(int);
			int* tmp_buf = (int*)malloc( size );
			memset( tmp_buf, 0, size );
			fwrite( tmp_buf, sizeof(int), Ntrees + 2, save_fd[n] );
			free( tmp_buf );
		}

     
    for(i = 0; i < NumGals; i++)
    {
      if(HaloGal[i].SnapNum == ListOutputSnaps[n]) 
      { 	       
      
        //prepare_galaxy_for_output(filenr, tree, &HaloGal[i], &galaxy_output);
        //myfwrite(&galaxy_output, sizeof(struct GALAXY_OUTPUT), 1, save_fd[n]);
        write_gridarray(&HaloGal[i], save_fd[n]); // Input snapshots are ordered highest -> lowest so it'll be 0th element. 
       	TotGalaxies[n]++;
        TreeNgals[n][tree]++;

      }
    }
  }

  // don't forget to free the workspace.
  free( OutputGalOrder );

}

void prepare_galaxy_for_output(int filenr, int tree, struct GALAXY *g, struct GALAXY_OUTPUT *o)
{
  int j, step;

  o->SnapNum = g->SnapNum;
  o->Type = g->Type;
  
  o->GalaxyIndex = 0;
  o->CentralGalaxyIndex = 0;
  o->SAGETreeIndex = 0;

  if(LastFile>=10000) // Assume that because there are so many files, the trees per file will be less than 100000.  Required for limits of long long.
  {
      assert( g->GalaxyNr < TREE_MUL_FAC ); // breaking tree size assumption
      assert(tree < (FILENR_MUL_FAC/10)/TREE_MUL_FAC);
      o->GalaxyIndex = g->GalaxyNr + TREE_MUL_FAC * tree + (FILENR_MUL_FAC/10) * filenr;
      assert( (o->GalaxyIndex - g->GalaxyNr - TREE_MUL_FAC*tree)/(FILENR_MUL_FAC/10) == filenr );
      assert( (o->GalaxyIndex - g->GalaxyNr -(FILENR_MUL_FAC/10)*filenr) / TREE_MUL_FAC == tree );
      assert( o->GalaxyIndex - TREE_MUL_FAC*tree - (FILENR_MUL_FAC/10)*filenr == g->GalaxyNr );
      o->CentralGalaxyIndex = HaloGal[HaloAux[Halo[g->HaloNr].FirstHaloInFOFgroup].FirstGalaxy].GalaxyNr + TREE_MUL_FAC * tree + (FILENR_MUL_FAC/10) * filenr;
      o->SAGETreeIndex = tree;
  }
  else 
  {
      assert( g->GalaxyNr < TREE_MUL_FAC ); // breaking tree size assumption
      assert(tree < FILENR_MUL_FAC/TREE_MUL_FAC);
      o->GalaxyIndex = g->GalaxyNr + TREE_MUL_FAC * tree + FILENR_MUL_FAC * filenr;
      assert( (o->GalaxyIndex - g->GalaxyNr - TREE_MUL_FAC*tree)/FILENR_MUL_FAC == filenr );
      assert( (o->GalaxyIndex - g->GalaxyNr -FILENR_MUL_FAC*filenr) / TREE_MUL_FAC == tree );
      assert( o->GalaxyIndex - TREE_MUL_FAC*tree - FILENR_MUL_FAC*filenr == g->GalaxyNr );
      o->CentralGalaxyIndex = HaloGal[HaloAux[Halo[g->HaloNr].FirstHaloInFOFgroup].FirstGalaxy].GalaxyNr + TREE_MUL_FAC * tree + FILENR_MUL_FAC * filenr;
      o->SAGETreeIndex = tree;
  }
    
  o->SAGEHaloIndex = g->HaloNr;
 
  o->mergeType = g->mergeType;
  o->mergeIntoID = g->mergeIntoID;
  o->mergeIntoSnapNum = g->mergeIntoSnapNum;
  o->dT = g->dT * UnitTime_in_s / SEC_PER_MEGAYEAR;

  for(j = 0; j < 3; j++)
  {
    o->Pos[j] = g->Pos[j];
    o->Vel[j] = g->Vel[j];
    o->Spin[j] = Halo[g->HaloNr].Spin[j];
  }

  o->Len = g->Len;
  o->Mvir = g->Mvir;
  o->Vmax = g->Vmax;

  o->CentralMvir = get_virial_mass(Halo[g->HaloNr].FirstHaloInFOFgroup);
  o->Rvir = get_virial_radius(g->HaloNr);  // output the actual Rvir, not the maximum Rvir
  o->Vvir = get_virial_velocity(g->HaloNr);  // output the actual Vvir, not the maximum Vvir
  //o->VelDisp = Halo[g->HaloNr].VelDisp;
  //o->SimulationFOFHaloIndex = Halo[g->HaloNr].SubhaloIndex;

  o->VelDisp = 0.0; 
  o->SimulationFOFHaloIndex = -1; 

   
  o->ColdGas = g->ColdGas;
  o->StellarMass = g->StellarMass;
  
  o->BulgeMass = g->BulgeMass;
  o->HotGas = g->HotGas;
  o->EjectedMass = g->EjectedMass;
  o->BlackHoleMass = g->BlackHoleMass;
  o->ICS = g->ICS;

  o->MetalsColdGas = g->MetalsColdGas;
  o->MetalsStellarMass = g->MetalsStellarMass;
  o->MetalsBulgeMass = g->MetalsBulgeMass;
  o->MetalsHotGas = g->MetalsHotGas;
  o->MetalsEjectedMass = g->MetalsEjectedMass;
  o->MetalsICS = g->MetalsICS;
  
  o->SfrDisk = 0.0;
  o->SfrBulge = 0.0;
  o->SfrDiskZ = 0.0;
  o->SfrBulgeZ = 0.0;
 
  // NOTE: in Msun/yr 
  for(step = 0; step < STEPS; step++)
  {
    o->SfrDisk += g->SfrDisk[step] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS / STEPS;
    o->SfrBulge += g->SfrBulge[step] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS / STEPS;

    if(g->SfrDiskColdGas[step] > 0.0)
      o->SfrDiskZ += g->SfrDiskColdGasMetals[step] / g->SfrDiskColdGas[step] / STEPS;

    if(g->SfrBulgeColdGas[step] > 0.0)
      o->SfrBulgeZ += g->SfrBulgeColdGasMetals[step] / g->SfrBulgeColdGas[step] / STEPS;
  }

  o->DiskScaleRadius = g->DiskScaleRadius;

  if (g->Cooling > 0.0)
    o->Cooling = log10(g->Cooling * UnitEnergy_in_cgs / UnitTime_in_s);
  else
    o->Cooling = 0.0;
  if (g->Heating > 0.0)
    o->Heating = log10(g->Heating * UnitEnergy_in_cgs / UnitTime_in_s);
  else
    o->Heating = 0.0;

  o->QuasarModeBHaccretionMass = g->QuasarModeBHaccretionMass;

  o->TimeOfLastMajorMerger = g->TimeOfLastMajorMerger * UnitTime_in_Megayears;
  o->TimeOfLastMinorMerger = g->TimeOfLastMinorMerger * UnitTime_in_Megayears;
	
  o->OutflowRate = g->OutflowRate * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;

  //infall properties
  if(g->Type != 0)
  {
    o->infallMvir = g->infallMvir;
    o->infallVvir = g->infallVvir;
    o->infallVmax = g->infallVmax;
  }
  else
  {
    o->infallMvir = 0.0;
    o->infallVvir = 0.0;
    o->infallVmax = 0.0;
  }

}

void write_gridarray(struct GALAXY *g, FILE *fp)
{

  int j;  
  size_t nwritten;

  ++count;
  
  nwritten = fwrite(&g->HaloNr, sizeof(g->HaloNr), 1, fp);
 
  XASSERT(g->IsMalloced == 1, "We are trying to write out the grid arrays for a galaxies who has already been freed.\n");
 
  float SFR_conversion = UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS / STEPS; 
 
  XASSERT(g->GridHistory != NULL, "GridHistory has a NULL pointer.\n"); 
  nwritten = fwrite(g->GridHistory, sizeof(*(g->GridHistory)), MAXSNAPS, fp);
  XASSERT( nwritten == MAXSNAPS, "Error: While writing GridHistory, expected to write %d times but wrote %zu times instead\n",
	   MAXSNAPS, nwritten);

  XASSERT(g->GridStellarMass != NULL, "GridStellarMass has a NULL pointer.\n"); 
  nwritten = fwrite(g->GridStellarMass, sizeof(*(g->GridStellarMass)), MAXSNAPS, fp);
  XASSERT( nwritten == MAXSNAPS, "Error: While writing GridStellarMass, expected to write %d times but wrote %zu times instead\n",
	   MAXSNAPS, nwritten);
 
  XASSERT(g->GridSFR != NULL, "GridSFR has a NULL pointer.\n"); 
 
  float *SFR_tmp = malloc(sizeof(*(g->GridSFR)) * MAXSNAPS);
  for (j = 0; j < MAXSNAPS; ++j)
  {
    SFR_tmp[j] = g->GridSFR[j]*SFR_conversion;   
  }
 
  nwritten = fwrite(SFR_tmp, sizeof(*(g->GridSFR)), MAXSNAPS, fp);
  XASSERT( nwritten == MAXSNAPS, "Error: While writing GridSFR, expected to write %d times but wrote %zu times instead\n",
	   MAXSNAPS, nwritten);
  free(SFR_tmp);  

  XASSERT(g->GridZ != NULL, "GridZ has a NULL pointer.\n"); 
  nwritten = fwrite(g->GridZ, sizeof(*(g->GridZ)), MAXSNAPS, fp);
  XASSERT( nwritten == MAXSNAPS, "Error: While writing GridZ, expected to write %d times but wrote %zu times instead\n",
	   MAXSNAPS, nwritten);

  XASSERT(g->GridCentralGalaxyMass != NULL, "GridCentralGalaxyMass has a NULL pointer.\n");  
  nwritten = fwrite(g->GridCentralGalaxyMass, sizeof(*(g->GridCentralGalaxyMass)), MAXSNAPS, fp);
  XASSERT( nwritten == MAXSNAPS, "Error: While writing GridCentralGalaxyMass, expected to write %d times but wrote %zu times instead\n",
	   MAXSNAPS, nwritten);

  XASSERT(g->MfiltGnedin != NULL, "MfiltGnedin has a NULL pointer.\n"); 
  nwritten = fwrite(g->MfiltGnedin, sizeof(*(g->MfiltGnedin)), MAXSNAPS, fp);
  XASSERT( nwritten == MAXSNAPS, "Error: While writing MfiltGnedin, expected to write %d times but wrote %zu times instead\n",
	   MAXSNAPS, nwritten);

  XASSERT(g->MfiltSobacchi != NULL, "MfilSobacchi has a NULL pointer.\n"); 
  nwritten = fwrite(g->MfiltSobacchi, sizeof(*(g->MfiltSobacchi)), MAXSNAPS, fp);
  XASSERT( nwritten == MAXSNAPS, "Error: While writing MfiltSobacchi, expected to write %d times but wrote %zu times instead\n",
	   MAXSNAPS, nwritten);

  XASSERT(g->EjectedFraction != NULL, "EjectedFraction has a NULL pointer.\n"); 
  nwritten = fwrite(g->EjectedFraction, sizeof(*(g->EjectedFraction)), MAXSNAPS, fp);
  XASSERT( nwritten == MAXSNAPS, "Error: While writing EjectedFraction, expected to write %d times but wrote %zu times instead\n",
	   MAXSNAPS, nwritten);

  for (j = 0; j < MAXSNAPS; ++j)
  {
    if (j == 80 && g->LenHistory[j] < 1 && g->LenHistory[j] != -1 )
    {
      fprintf(stderr, "%d\n", g->LenHistory[j]);  
      exit(EXIT_FAILURE);
    }
  }

  XASSERT(g->LenHistory != NULL, "LenHistory has a NULL pointer.\n"); 
  nwritten = fwrite(g->LenHistory, sizeof(*(g->LenHistory)), MAXSNAPS, fp);
  XASSERT( nwritten == MAXSNAPS, "Error: While writing LenHistory, expected to write %d times but wrote %zu times instead\n",
	   MAXSNAPS, nwritten);
 
  float *outflow_tmp= malloc(sizeof(*(g->GridOutflowRate)) * MAXSNAPS);
  for (j = 0; j < MAXSNAPS; ++j)
  {
    outflow_tmp[j] = g->GridOutflowRate[j] * SFR_conversion; 
    XASSERT(outflow_tmp[j] == outflow_tmp[j], "Outflow_tmp has a value of %.4e for snapshot %d. The GridOutflowRate value was %.4e\n", outflow_tmp[j], j, g->GridOutflowRate[j]);
  }

  nwritten = fwrite(outflow_tmp, sizeof(*(g->GridOutflowRate)), MAXSNAPS, fp);
  XASSERT( nwritten == MAXSNAPS, "Error: While writing GridOutflowRate, expected to write %d times but wrote %zu times instead\n",
	   MAXSNAPS, nwritten);
  free(outflow_tmp);  

  float *infall_tmp= malloc(sizeof(*(g->GridInfallRate)) * MAXSNAPS);
  XASSERT(g->GridInfallRate != NULL, "GridInfallRate has a NULL pointer.\n");
  for (j = 0; j < MAXSNAPS; ++j)
  {
    infall_tmp[j] = g->GridInfallRate[j] * SFR_conversion; 
    XASSERT(infall_tmp[j] == infall_tmp[j], "infall_tmp has a value of %.4e for snapshot %d. The GridOutflowRate value was %.4e\n", infall_tmp[j], j, g->GridInfallRate[j]);
  }

  nwritten = fwrite(infall_tmp, sizeof(*(g->GridInfallRate)), MAXSNAPS, fp);
  XASSERT( nwritten == MAXSNAPS, "Error: While writing GridInfallRate, expected to write %d times but wrote %zu times instead\n",
	   MAXSNAPS, nwritten);
  free(infall_tmp);

  XASSERT(g->GridEjectedMass != NULL, "GridEjectedMass has a NULL pointer.\n"); 
  nwritten = fwrite(g->GridEjectedMass, sizeof(*(g->GridEjectedMass)), MAXSNAPS, fp);
  XASSERT( nwritten == MAXSNAPS, "Error: While writing GridEjectedMass, expected to write %d times but wrote %zu times instead\n",
	   MAXSNAPS, nwritten);

  XASSERT(g->QuasarActivity != NULL, "QuasarActivity has a NULL pointer.\n"); 
  nwritten = fwrite(g->QuasarActivity, sizeof(*(g->QuasarActivity)), MAXSNAPS, fp);
  XASSERT( nwritten == MAXSNAPS, "Error: While writing QuasarActivity, expected to write %d times but wrote %zu times instead\n",
	   MAXSNAPS, nwritten);

  float *dynamical_tmp = malloc(sizeof(*(g->DynamicalTime)) * MAXSNAPS);
  XASSERT(g->DynamicalTime != NULL, "DynamicalTime has a NULL pointer.\n"); 
  for (j = 0; j < MAXSNAPS; ++j)
  {
    dynamical_tmp[j] = g->DynamicalTime[j] * UnitTime_in_Megayears;
  }
  
  nwritten = fwrite(dynamical_tmp, sizeof(*(g->DynamicalTime)), MAXSNAPS, fp);
  XASSERT( nwritten == MAXSNAPS, "Error: While writing DynamicalTime, expected to write %d times but wrote %zu times instead\n",
	   MAXSNAPS, nwritten);
  free(dynamical_tmp);

  XASSERT(g->QuasarSubstep != NULL, "QuasarSubstep has a NULL pointer.\n"); 
  nwritten = fwrite(g->QuasarSubstep, sizeof(*(g->QuasarSubstep)), MAXSNAPS, fp);
  XASSERT( nwritten == MAXSNAPS, "Error: While writing QuasarSubstep, expected to write %d times but wrote %zu times instead\n",
	   MAXSNAPS, nwritten);

  XASSERT(g->GridColdGas != NULL, "GridColdGas has a NULL pointer.\n"); 
  nwritten = fwrite(g->GridColdGas, sizeof(*(g->GridColdGas)), MAXSNAPS, fp);
  XASSERT( nwritten == MAXSNAPS, "Error: While writing GridColdGas, expected to write %d times but wrote %zu times instead\n",
	   MAXSNAPS, nwritten);
 
  ++times_written;

}

void finalize_galaxy_file(int filenr)
{
  int n;

  for(n = 0; n < NOUT; n++)
  {
    // file must already be open.
    assert( save_fd[n] );

    // seek to the beginning.
    fseek( save_fd[n], 0, SEEK_SET );

    fprintf(stderr, "Ntrees = %d \t TotGalaxies[n] = %d\n", Ntrees, TotGalaxies[n]);
    printf("n = %d\n", n);
    
    myfwrite(&Ntrees, sizeof(int), 1, save_fd[n]); 
    myfwrite(&TotGalaxies[n], sizeof(int), 1, save_fd[n]);
    myfwrite(TreeNgals[n], sizeof(int), Ntrees, save_fd[n]);

//    printf("Ntrees = %d\nTotnGalaxies[0] = %d\nTreeNGals[0] = %d\n", Ntrees, TotGalaxies[0], TreeNgals[0][0]);
    // close the file and clear handle after everything has been written
    fclose( save_fd[n] );
    save_fd[n] = NULL;
  }
  
}


void save_merged_galaxies(int filenr, int tree)
{
  char buf[1000];
  int i;
  //struct GALAXY_OUTPUT galaxy_output;

  int OutputGalCount, *OutputGalOrder;
  
  OutputGalOrder = (int*)malloc( MergedNr*sizeof(int) );
  assert( OutputGalOrder );

  // reset the output galaxy count and order 
  OutputGalCount = 0;
  for(i = 0; i < MergedNr; i++)
  {
    OutputGalOrder[i] = i;
    OutputGalCount++;
  }

  if( !save_fd2 )
  { 
    sprintf(buf, "%s/%s_MergedGalaxies_%d", OutputDir, FileNameGalaxies, filenr);

    if(!(save_fd2 = fopen(buf, "r+")))
    {
      printf("can't open file `%s'\n", buf);		
      ABORT(0);
    }

    // write out placeholders for the header data.
    size_t size = (Ntrees + 2)*sizeof(int);
    int* tmp_buf = (int*)malloc( size );
    memset( tmp_buf, 0, size );
    fwrite( tmp_buf, sizeof(int), Ntrees + 2, save_fd2);
    free( tmp_buf );
  }

  for(i = 0; i < MergedNr; i++)
  {
    //prepare_galaxy_for_output(filenr, tree, &MergedGal[i], &galaxy_output);
    //myfwrite(&galaxy_output, sizeof(struct GALAXY_OUTPUT), 1, save_fd2);
    write_gridarray(&MergedGal[i], save_fd2);

    TotMerged++;
    TreeNMergedgals[tree]++;
  }

  // don't forget to free the workspace.
  free( OutputGalOrder );

}

void finalize_merged_galaxy_file(int filenr)
{

  // file must already be open.
  assert( save_fd2 );

  // seek to the beginning.
  fseek( save_fd2, 0, SEEK_SET );

//  printf("TotMerged = %d.  TotGalaxies = %d\n", TotMerged, TotGalaxies[0]);
  myfwrite(&Ntrees, sizeof(int), 1, save_fd2);
  myfwrite(&TotMerged, sizeof(int), 1, save_fd2);
  myfwrite(TreeNMergedgals, sizeof(int), Ntrees, save_fd2);

  // close the file and clear handle after everything has been written
  fclose( save_fd2 );
  save_fd2 = NULL;
  
}

#undef TREE_MUL_FAC
#undef FILENR_MUL_FAC

