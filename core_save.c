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
  struct GALAXY_OUTPUT galaxy_output;

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
      
        prepare_galaxy_for_output(filenr, tree, &HaloGal[i], &galaxy_output);
        myfwrite(&galaxy_output, sizeof(struct GALAXY_OUTPUT), 1, save_fd[n]);
        if(HaloGal[i].SnapNum == ListOutputSnaps[0]) // We only want to write out the grid properties for the final snapshot.
            write_gridarray(&HaloGal[i], save_fd[n]); // Input snapshots are ordered highest -> lowest so it'll be 0th element. 
        //fprintf(stderr, "On write the GrandSum was %.4e and Stellar Mass was %.4e for Galaxy %d\n", HaloGal[i].GrandSum, HaloGal[i].StellarMass, i);
        /*
        if(HaloGal[i].StellarMass != 0)
	{ 
		//double ratio = HaloGal[i].GrandSum/HaloGal[i].StellarMass;
		//fprintf(stderr, "Ratio: %.4e \t StellarMass = %.4e \t GrandSum = %.4e\n", ratio, HaloGal[i].StellarMass, HaloGal[i].GrandSum); 	
		//fprintf(stderr, "%.4f\n", ratio);

        }
        */
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
  float SFR; 

  ++count;

  XASSERT(g->IsMalloced == 1, "We are trying to write out the grid arrays for a galaxies who has already been freed.\n");
  XASSERT(!(!(fp)), "We are trying to write to a file which has a NULL pointer.\n"); 
  double SFR_conversion = UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS / STEPS; 
 
  XPRINT(!(!(g->GridHistory)), "GridHistory has a NULL pointer.\n"); 
  for (j = 0; j < MAXSNAPS; ++j)
  { 
      myfwrite(&g->GridHistory[j], sizeof(*(g->GridHistory)), 1, fp);
  }

  XPRINT(!(!(g->GridStellarMass)), "GridStellarMass has a NULL pointer.\n"); 
  for (j = 0; j < MAXSNAPS; ++j)
  { 
      myfwrite(&g->GridStellarMass[j], sizeof(*(g->GridStellarMass)), 1, fp);
  }
 
  XPRINT(!(!(g->GridSFR)), "GridSFR has a NULL pointer.\n"); 
  for (j = 0; j < MAXSNAPS; ++j)
  {
      SFR = g->GridSFR[j]*SFR_conversion;    
      myfwrite(&SFR, sizeof(*(g->GridSFR)), 1, fp);
  }

  XPRINT(!(!(g->GridZ)), "GridZ has a NULL pointer.\n"); 
  for (j = 0; j < MAXSNAPS; ++j)
  {
      XPRINT(j < MAXSNAPS, "Trying to save Mfilt Gnedin and have index j = %d \t MAXSNAPS = %d\n", j, MAXSNAPS);


      myfwrite(&g->GridZ[j], sizeof(*(g->GridZ)), 1, fp);
  }

  XPRINT(!(!(g->GridCentralGalaxyMass)), "GridCentralGalaxyMass has a NULL pointer.\n"); 
  for (j = 0; j < MAXSNAPS; ++j)
  {
      myfwrite(&g->GridCentralGalaxyMass[j], sizeof(*(g->GridCentralGalaxyMass)), 1, fp);
  }
 
  for (j = 0; j < MAXSNAPS; ++j)
  {
      float tmp = -999.99;
      myfwrite(&tmp, sizeof(float), 1, fp);
  }

  XPRINT(!(!(g->MfiltGnedin)), "MfiltGnedin has a NULL pointer.\n"); 
  for (j = 0; j < MAXSNAPS; ++j)
  {
      XPRINT(j < MAXSNAPS, "Trying to save Mfilt Gnedin and have index j = %d \t MAXSNAPS = %d\n", j, MAXSNAPS); 
      myfwrite(&g->MfiltGnedin[j], sizeof(*(g->MfiltGnedin)), 1, fp);  
      //XASSERT(g->MfiltGnedin[j] == 1.0, "Somehow got a value for MfiltGned = %.4e.  Snapshot %d\n", g->MfiltGnedin[j], j);
  }

  XPRINT(!(!(g->MfiltSobacchi)), "MfilSobacchi has a NULL pointer.\n"); 
  for (j = 0; j < MAXSNAPS; ++j)
  { 
      myfwrite(&g->MfiltSobacchi[j], sizeof(*(g->MfiltSobacchi)), 1, fp);   
  }

  XPRINT(!(!(g->EjectedFraction)), "EjectedFraction has a NULL pointer.\n"); 
  for (j = 0; j < MAXSNAPS; ++j)
  {  
     myfwrite(&g->EjectedFraction[j], sizeof(*(g->EjectedFraction)), 1, fp); 
  }

  XPRINT(!(!(g->LenHistory)), "LenHistory has a NULL pointer.\n"); 
  for (j = 0; j < MAXSNAPS; ++j)
  {  
     myfwrite(&g->LenHistory[j], sizeof(*(g->LenHistory)), 1, fp);  
  }

  /*
  // Padding // 
  for (j = 0; j < MAXSNAPS; ++j)
  {
      float tmpOutflowRate = g->GridOutflowRate[j] * SFR_conversion * STEPS; 
      if(tmpOutflowRate > 1e4 && j == 63)
	fprintf(stderr, "tmpOutflowRate = %.4e \t g->GridOutflowRate[j] = %.4e \t SFR_conversion * STEPS = %.4e \t Gal[p].OutflowRate = %.4e\n", tmpOutflowRate, g->GridOutflowRate[j], SFR_conversion * STEPS, g->OutflowRate); 
      myfwrite(&tmpOutflowRate, sizeof(float), 1, fp); 
  }
  */

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
  struct GALAXY_OUTPUT galaxy_output;

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
    prepare_galaxy_for_output(filenr, tree, &MergedGal[i], &galaxy_output);
    myfwrite(&galaxy_output, sizeof(struct GALAXY_OUTPUT), 1, save_fd2);
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

