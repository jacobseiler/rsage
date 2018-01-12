#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"


// keep a static file handle to remove the need to do constant seeking
FILE* load_fd = NULL;


void load_tree_table(int filenr)
{
  int i, n, totNHalos;
  char buf[MAXLEN], tag[MAXLEN];
  FILE *fd;

  if(use_tiamat == 1)
  {
    if(filenr < 10)
      snprintf(tag, MAXLEN, "00%d", filenr);
    else if (filenr >= 10 && filenr < 100)
      snprintf(tag, MAXLEN, "0%d", filenr);
    else
      snprintf(tag, MAXLEN, "%d", filenr);
    sprintf(buf, "%s/%s_%s.dat", SimulationDir, TreeName, tag);
  }
  else
  {  
    snprintf(buf, MAXLEN, "%s/%s.%d", SimulationDir, TreeName, filenr);
  }
  fprintf(stderr, "Reading file %s\n", buf);
  if(!(load_fd = fopen(buf, "r")))
  {
    printf("can't open file `%s'\n", buf);
    ABORT(0);
  }

  myfread(&Ntrees, 1, sizeof(int), load_fd);
  myfread(&totNHalos, 1, sizeof(int), load_fd);

  TreeNHalos = mymalloc(sizeof(int) * Ntrees);
  TreeFirstHalo = mymalloc(sizeof(int) * Ntrees);

  for(n = 0; n < NOUT; n++)
    TreeNgals[n] = mymalloc(sizeof(int) * Ntrees);
  TreeNMergedgals = mymalloc(sizeof(int)* Ntrees);
  myfread(TreeNHalos, Ntrees, sizeof(int), load_fd); 

  if(Ntrees)
    TreeFirstHalo[0] = 0;
  for(i = 1; i < Ntrees; i++)
    TreeFirstHalo[i] = TreeFirstHalo[i - 1] + TreeNHalos[i - 1];

  for(n = 0; n < NOUT; n++)
  {
    for(i = 0; i < Ntrees; i++)
    {
      TreeNgals[n][i] = 0;
      TreeNMergedgals[i] = 0;
    }
    sprintf(buf, "%s/%s_z%1.3f_%d", OutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[n]], filenr);

    if(!(fd = fopen(buf, "w")))
    {
      printf("can't open file `%s'\n", buf);
      ABORT(0);
    }
    fclose(fd);
    TotGalaxies[n] = 0;
  }
  TotMerged = 0;
  fprintf(stderr, "Read the table\n");
}



void free_tree_table(void)
{
  int n;

  myfree(TreeNMergedgals); 
  for(n = NOUT - 1; n >= 0; n--)
    myfree(TreeNgals[n]);

  myfree(TreeFirstHalo);
  myfree(TreeNHalos);
	
	// Don't forget to free the open file handle
	if(load_fd) {
		fclose(load_fd);
		load_fd = NULL;
	}
}



void load_tree(int filenr, int nr)
{
  int i;

  // must have an FD
  assert( load_fd );

  Halo = mymalloc(sizeof(struct halo_data) * TreeNHalos[nr]);  
  myfread(Halo, TreeNHalos[nr], sizeof(struct halo_data), load_fd);

  MaxGals = (int)(MAXGALFAC * TreeNHalos[nr]);

  if(MaxGals < 10000)  
    MaxGals = 10000;

  MaxMergedGals = MaxGals;
  FoF_MaxGals = 10000; 

  gal_to_free = malloc(sizeof(int) * MaxMergedGals);
  HaloAux = mymalloc(sizeof(struct halo_aux_data) * TreeNHalos[nr]);
  HaloGal = mymalloc(sizeof(struct GALAXY) * MaxGals);
  Gal = mymalloc(sizeof(struct GALAXY) * FoF_MaxGals);
  MergedGal = mymalloc(sizeof(struct GALAXY) * MaxMergedGals);   

  double Min_Halo = 1e5;
  double Max_Halo = 0.0; 
  for(i = 0; i < TreeNHalos[nr]; i++)
  {
    if(Halo[i].Mvir > Max_Halo)
      Max_Halo = Halo[i].Mvir;
    if(Halo[i].Mvir < Min_Halo)
      Min_Halo = Halo[i].Mvir;
   
    //if(Halo[i].NextHaloInFOFgroup != -1)
     // fprintf(stderr, "%d\n", nr);
  
    //if(Halo[i].FirstProgenitor != -1)
    //XASSERT(Halo[i].SnapNum != Halo[Halo[i].FirstProgenitor].SnapNum, "Halo[%d].SnapNum = %d, Halo[%d].FirstProgenitor = %d, Halo[Halo[%d].FirstProgenitor].SnapNum = %d\n", i, Halo[i].SnapNum, i, Halo[i].FirstProgenitor, i, Halo[Halo[i].FirstProgenitor].SnapNum);
  /* 
    if (Halo[i].FirstProgenitor== -1)
    {
      fprintf(stderr, "Halo: %d. FirstProg = %d, NextProg = %d, FirstHaloInFOFgroup = %d, NextHaloInFOFgroup = %d SnapNum = %d, Mvir = %.4e\n", i, Halo[i].FirstProgenitor, Halo[i].NextProgenitor, Halo[i].FirstHaloInFOFgroup, Halo[i].NextHaloInFOFgroup, Halo[i].SnapNum, Halo[i].Mvir);
    }
    */
    HaloAux[i].DoneFlag = 0;
    HaloAux[i].HaloFlag = 0;
    HaloAux[i].NGalaxies = 0;
//    if (Halo[i].SnapNum != LastSnapShotNr && Halo[i].Descendant == -1)
//fprintf(stderr, "STUPID TREES! Halo[%d].SnapNum = %d Halo[%d].Descendant = %d\n", i, Halo[i].SnapNum, i, Halo[i].Descendant); 
  }
//fprintf(stderr, "Min_Halo = %.4f, Max_Halo = %.4f\n", Min_Halo, Max_Halo); 
//  exit(EXIT_FAILURE);
}



void free_galaxies_and_tree(void)
{
  int i, j, max_snap, count_frees = 0;

  // This block is quite important (and took me a ridiculously long time to figure out) so I'll explain it. 
  // The issue is that Greg's tree building code allows so-called 'stray' or 'dangling' halos.  These are halos which do NOT have descendants but are not at the root redshift.
  // Because a progenitor line all points to the same block of memory for arrays such as GridHistory, we must only free this memory block once. 
  // While normally we could just ask 'is this Halo at the root redshift?' and free it if the answer is yes, the stray halos are not caught by this if condition.
  //
  // To correctly account for the stray halos, we execute the following code for each galaxy:
  // First we find out what snapshot is the last snapshot that this progenitor line is alive. Since the pointers all point to the same memory block, the galaxy at snapshot 20 knows if it will be alive at snapshot 50 and will hence know this maximum snapshot.
  // Once we find the final snapshot of the progenitor line we ask 'is this galaxy THE galaxy that is alive at the final snapshot?'.  
  // If this condition is fulfilled, we add the galaxy index to an array.  WE DO NOT FREE THE GALAXY IMMEDIATELY because we will still be checking the other galaxies in this progenitor line.
  // Once we have checked through all galaxies and constructed the indices that should be freed, we then finally do the freeing.

  for(i = 0; i < NumGals; ++i)
  {
    max_snap = 0;
    XASSERT(HaloGal[i].IsMalloced == 1, "HaloGal %d doesn't have grids mallocced but we're trying to free it.\n", i);

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
      XPRINT(HaloGal[i].IsMalloced == 1, "HaloGal %d doesn't have grids mallocced but we're trying to free it.\n", i); 
      gal_to_free[count_frees] = i;
      ++count_frees;
    } 
  }

  for(i = 0; i < count_frees; ++i)
  {

    free_grid_arrays(&HaloGal[gal_to_free[i]]);
    ++gal_frees; 
  }
  free(gal_to_free);


  // Now we just have to free the arrays for the galaxies that have merged.

  for(i = 0; i < MergedNr; ++i)
  {
      XPRINT(MergedGal[i].IsMalloced == 1, "MergedGal %d doesn't have grids mallocced but we're trying to free it.\n", i); 
    free_grid_arrays(&MergedGal[i]); // These are the Gal[xxx] entries that were copied over to  
    ++mergedgal_frees;
  } 

  // All the inside pointers have now been freed, lets free the structs themselves now.
  myfree(MergedGal);
  myfree(Gal);
  myfree(HaloGal);
  myfree(HaloAux);
  myfree(Halo);
}

void free_grid_arrays(struct GALAXY *g)
{
  free(g->GridHistory);
  free(g->GridStellarMass);
  free(g->GridSFR);
  free(g->GridZ);
  free(g->GridCentralGalaxyMass);
  free(g->MfiltGnedin);
  free(g->MfiltSobacchi);
  free(g->EjectedFraction);
  free(g->LenHistory);
  free(g->Stars);
  free(g->GridOutflowRate);
  free(g->GridInfallRate);
  free(g->GridEjectedMass);
  free(g->QuasarActivity);
  free(g->DynamicalTime);

  g->IsMalloced = 0;
}

void malloc_grid_arrays(struct GALAXY *g)
{
  g->GridHistory = malloc(sizeof(*(g->GridHistory)) * (MAXSNAPS));
  if(g->GridHistory == NULL)
  {
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridHistory.\n", sizeof(*(g->GridHistory))*MAXSNAPS); 
    exit(EXIT_FAILURE);
  }

  g->GridStellarMass = malloc(sizeof(*(g->GridStellarMass)) * (MAXSNAPS));
  if(g->GridStellarMass == NULL)
  {
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridStellarMass.\n", sizeof(*(g->GridStellarMass))*MAXSNAPS); 
    exit(EXIT_FAILURE);
  } 

  g->GridSFR = malloc(sizeof(*(g->GridSFR)) * (MAXSNAPS));
  if(g->GridSFR == NULL)
  {
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridSFR.\n", sizeof(*(g->GridSFR))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  g->GridZ = malloc(sizeof(*(g->GridZ)) * (MAXSNAPS));
  if (g->GridZ == NULL)
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridSFR.\n", sizeof(*(g->GridZ))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }
 
  g->GridCentralGalaxyMass = malloc(sizeof(*(g->GridCentralGalaxyMass)) * (MAXSNAPS));
  if (g->GridCentralGalaxyMass == NULL)
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridCentralGalaxyMass.\n", sizeof(*(g->GridCentralGalaxyMass))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  g->MfiltGnedin = malloc(sizeof(*(g->MfiltGnedin)) * (MAXSNAPS));
  if (g->MfiltGnedin == NULL)
  {
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate MfiltGnedin.\n", sizeof(*(g->MfiltGnedin))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  g->MfiltSobacchi = malloc(sizeof(*(g->MfiltSobacchi)) * (MAXSNAPS));
  if (g->MfiltSobacchi == NULL)
  {   
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate MfiltSobacchi.\n", sizeof(*(g->MfiltSobacchi))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  g->EjectedFraction = malloc(sizeof(*(g->EjectedFraction)) * (MAXSNAPS));
  if (g->EjectedFraction == NULL)
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate EjectedFraction.\n", sizeof(*(g->EjectedFraction))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  g->LenHistory = malloc(sizeof(*(g->LenHistory)) * (MAXSNAPS));
  if (g->LenHistory == NULL)
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate LenHistory.\n", sizeof(*(g->LenHistory))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  g->Stars = malloc(sizeof(*(g->Stars)) * SN_Array_Len);
  if (g->Stars == NULL)
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate Stars.\n", sizeof(*(g->Stars))*SN_Array_Len);
    exit(EXIT_FAILURE);
  }

  g->GridOutflowRate = malloc(sizeof(*(g->GridOutflowRate)) * (MAXSNAPS)); 
  if (g->GridOutflowRate == NULL)
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate OutflowRate.\n", sizeof(*(g->GridOutflowRate))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  g->GridInfallRate = malloc(sizeof(*(g->GridInfallRate)) * (MAXSNAPS)); 
  if (g->GridInfallRate == NULL)
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridInfallRate.\n", sizeof(*(g->GridInfallRate))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  g->GridEjectedMass = malloc(sizeof(*(g->GridEjectedMass)) * (MAXSNAPS)); 
  if (g->GridEjectedMass == NULL)
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridEjectedMass.\n", sizeof(*(g->GridEjectedMass))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  g->QuasarActivity = malloc(sizeof(*(g->QuasarActivity)) * (MAXSNAPS));
  if (g->QuasarActivity == NULL)
  {
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate QuasarActivity.\n", sizeof(*(g->QuasarActivity)) * MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  g->DynamicalTime = malloc(sizeof(*(g->DynamicalTime)) * (MAXSNAPS));
  if (g->DynamicalTime == NULL)
  {
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate DynamicalTime.\n", sizeof(*(g->DynamicalTime)) * MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  g->IsMalloced = 1;
}

int32_t free_grid()
{

  int32_t i;

  for (i = 0; i < Grid->NumGrids; ++i)
  {
    free(Grid->PhotoGrid[i].PhotoRate);
  }

  free(Grid->PhotoGrid);
  free(Grid->ReionRedshift);
  free(Grid);

  return EXIT_SUCCESS;

} 

size_t myfread(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  return fread(ptr, size, nmemb, stream);
}

size_t myfwrite(void *ptr, size_t size, size_t nmemb, FILE * stream)
{
  return fwrite(ptr, size, nmemb, stream);
}

int myfseek(FILE * stream, long offset, int whence)
{
  return fseek(stream, offset, whence);
}
