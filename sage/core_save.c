#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"

// keep a static file handle to remove the need to do constant seeking.
FILE* save_fd[ABSOLUTEMAXSNAPS] = { 0 };
FILE* save_fd2 = NULL;

void save_galaxies(int filenr, int tree)
{
  char buf[1024];
  int32_t i, n, max_snap, j;

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
    
    // There are some galaxies that aren't at the root redshift but do not have any descendants.
    // This block of code catches these excepetions.
    // A more detailed comment on this is located in the 'free_galaxies_and_tree' function in 'core_io_tree.c'. 
 
    for(i = 0; i < NumGals; i++)
    {
      max_snap = 0;

      if(HaloGal[i].IsMerged != -1)
        continue;
      for(j = 1; j < MAXSNAPS; ++j)
      {
        if(HaloGal[i].GridHistory[j] != -1)
        {    
          max_snap = j;
        } 
      }

      if(HaloGal[i].SnapNum == max_snap && Halo[HaloGal[i].HaloNr].Descendant == -1)
      {
          write_gridarray(&HaloGal[i], save_fd[n]); // Input snapshots are ordered highest -> lowest so it'll be 0th element. 
          TotGalaxies[n]++;
          TreeNgals[n][tree]++;          
      } 

    } // NumGals loop.
  } // NOUT loop.

  // don't forget to free the workspace.
  free( OutputGalOrder );

}

void write_gridarray(struct GALAXY *g, FILE *fp)
{

#define WRITE_GRID_PROPERTY(name, length)    \
{                                            \
  XASSERT(name != NULL, #name" has a NULL pointer.\n"); \
  nwritten = fwrite(name, sizeof(*(name)), length, fp); \
  XASSERT(nwritten == length, "While writing "#name", we expected to write %d times but wrote %d times instead\n", \
          length, nwritten);                 \
}

#define WRITE_CONVERTED_GRID_PROPERTY(name, length, conversion, type)    \
{                                            \
  XASSERT(name != NULL, #name" has a NULL pointer.\n"); \
  buffer = calloc(length, sizeof(type));\
  XASSERT(buffer != NULL, "Could not allocate memory for a buffer in write_gridarray for "#name".\n"); \
  for (j = 0; j < length; ++j)\
  {\
    ((type*)buffer)[j] = name[j] * conversion;\
  }\
  nwritten = fwrite(buffer, sizeof(type), length, fp); \
  XASSERT(nwritten == length, "While writing"#name", we expected to write %d times but wrote %d times instead\n", \
          length, nwritten);                 \
  free(buffer); \
}

  float SFR_conversion = UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS / STEPS; 
  int j;  
  int32_t nwritten;
  void *buffer;
 
  nwritten = fwrite(&g->HaloNr, sizeof(g->HaloNr), 1, fp);

  fwrite(&g->mergeType, sizeof(g->mergeType), 1, fp);
 
  XASSERT(g->IsMalloced == 1, "We are trying to write out the grid arrays for a galaxies who has already been freed.\n");
 

  WRITE_GRID_PROPERTY(g->GridHistory, MAXSNAPS);
  WRITE_GRID_PROPERTY(g->GridStellarMass, MAXSNAPS);

  WRITE_CONVERTED_GRID_PROPERTY(g->GridSFR, MAXSNAPS, SFR_conversion, typeof(*(g->GridSFR))); 

  WRITE_GRID_PROPERTY(g->GridZ, MAXSNAPS);
  WRITE_GRID_PROPERTY(g->GridCentralGalaxyMass, MAXSNAPS);
  WRITE_GRID_PROPERTY(g->EjectedFraction, MAXSNAPS);
  WRITE_GRID_PROPERTY(g->LenHistory, MAXSNAPS);

  WRITE_CONVERTED_GRID_PROPERTY(g->GridOutflowRate, MAXSNAPS, SFR_conversion, typeof(*(g->GridOutflowRate)));
  WRITE_CONVERTED_GRID_PROPERTY(g->GridInfallRate, MAXSNAPS, SFR_conversion, typeof(*(g->GridInfallRate)));

  WRITE_GRID_PROPERTY(g->GridEjectedMass, MAXSNAPS);
  WRITE_GRID_PROPERTY(g->QuasarActivity, MAXSNAPS);
 
  WRITE_CONVERTED_GRID_PROPERTY(g->DynamicalTime, MAXSNAPS, UnitTime_in_Megayears, typeof(*(g->DynamicalTime)));

  WRITE_GRID_PROPERTY(g->QuasarSubstep, MAXSNAPS);
  WRITE_GRID_PROPERTY(g->GridColdGas, MAXSNAPS);
  WRITE_GRID_PROPERTY(g->LenMergerGal, MAXSNAPS);
  WRITE_GRID_PROPERTY(g->GridBHMass, MAXSNAPS);
  WRITE_GRID_PROPERTY(g->GridReionMod, MAXSNAPS);
  WRITE_GRID_PROPERTY(g->GridDustColdGas, MAXSNAPS);
  WRITE_GRID_PROPERTY(g->GridDustHotGas, MAXSNAPS);
  WRITE_GRID_PROPERTY(g->GridDustEjectedMass, MAXSNAPS);

#undef WRITE_GRID_PROPERTY
#undef WRITE_CONVERTED_GRID_PROPERTY
 
  ++times_written;

}

void finalize_galaxy_file(void)
{
  int n;

  for(n = 0; n < NOUT; n++)
  {
    // file must already be open.
    assert( save_fd[n] );

    // seek to the beginning.
    fseek( save_fd[n], 0, SEEK_SET );
    
    myfwrite(&Ntrees, sizeof(int), 1, save_fd[n]); 
    myfwrite(&TotGalaxies[n], sizeof(int), 1, save_fd[n]);
    myfwrite(TreeNgals[n], sizeof(int), Ntrees, save_fd[n]);

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
    write_gridarray(&MergedGal[i], save_fd2);

    TotMerged++;
    TreeNMergedgals[tree]++;
  }

  // don't forget to free the workspace.
  free( OutputGalOrder );

}

void finalize_merged_galaxy_file(void)
{

  // file must already be open.
  assert( save_fd2 );

  // seek to the beginning.
  fseek( save_fd2, 0, SEEK_SET );

  myfwrite(&Ntrees, sizeof(int), 1, save_fd2);
  myfwrite(&TotMerged, sizeof(int), 1, save_fd2);
  myfwrite(TreeNMergedgals, sizeof(int), Ntrees, save_fd2);

  // close the file and clear handle after everything has been written
  fclose( save_fd2 );
  save_fd2 = NULL;
  
}
