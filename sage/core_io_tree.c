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
  char buf[MAXLEN];
  FILE *fd;
  
  snprintf(buf, MAXLEN, "%s/%s_%03d%s", SimulationDir, TreeName, filenr, TreeExtension);

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

}



void free_tree_table(void)
{
  int n;

  myfree(TreeNMergedgals, sizeof(*(TreeNMergedgals)) * Ntrees); 
  for(n = NOUT - 1; n >= 0; n--)
    myfree(TreeNgals[n], sizeof(*(TreeNgals[n])) * Ntrees);

  myfree(TreeFirstHalo, sizeof(*(TreeFirstHalo)) * Ntrees);
  myfree(TreeNHalos, sizeof(*(TreeNHalos)) * Ntrees);
	
	// Don't forget to free the open file handle
	if(load_fd) {
		fclose(load_fd);
		load_fd = NULL;
	}
}

void load_tree(int nr)
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

  gal_to_free = mymalloc(sizeof(int) * MaxMergedGals);
  HaloAux = mymalloc(sizeof(struct halo_aux_data) * TreeNHalos[nr]);
  HaloGal = mymalloc(sizeof(struct GALAXY) * MaxGals);
  Gal = mymalloc(sizeof(struct GALAXY) * FoF_MaxGals);
  MergedGal = mymalloc(sizeof(struct GALAXY) * MaxMergedGals);   

  double Min_Halo = 1e5;
  double Max_Halo = 0.0; 
  for(i = 0; i < TreeNHalos[nr]; i++)
  {
    if (Halo[i].FirstHaloInFOFgroup == -1)
    exit(0);

    if(Halo[i].Mvir > Max_Halo)
      Max_Halo = Halo[i].Mvir;
    if(Halo[i].Mvir < Min_Halo)
      Min_Halo = Halo[i].Mvir;
 
#ifdef BRITTON_SIM     
    Halo[i].Pos[0] = Halo[i].Pos[0] - 775.0;
    Halo[i].Pos[1] = Halo[i].Pos[1] - 775.0;
    Halo[i].Pos[2] = Halo[i].Pos[2] - 775.0;
#endif
  
    HaloAux[i].DoneFlag = 0;
    HaloAux[i].HaloFlag = 0;
    HaloAux[i].NGalaxies = 0;
  }
  
}

void free_galaxies_and_tree(int32_t treenr)
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
  myfree(gal_to_free, sizeof(*(gal_to_free)) * MaxMergedGals);


  // Now we just have to free the arrays for the galaxies that have merged.

  for(i = 0; i < MergedNr; ++i)
  {
      XPRINT(MergedGal[i].IsMalloced == 1, "MergedGal %d doesn't have grids mallocced but we're trying to free it.\n", i); 
    free_grid_arrays(&MergedGal[i]); // These are the Gal[xxx] entries that were copied over to  
    ++mergedgal_frees;
  } 

  // All the inside pointers have now been freed, lets free the structs themselves now.
  myfree(MergedGal, sizeof(*(MergedGal)) * MaxMergedGals);
  myfree(Gal, sizeof(*(Gal)) * FoF_MaxGals);
  myfree(HaloGal, sizeof(*(HaloGal)) * MaxGals);
  myfree(HaloAux, sizeof(*(HaloAux)) * TreeNHalos[treenr]);
  myfree(Halo, sizeof(*(Halo)) * TreeNHalos[treenr]);

}

void free_grid_arrays(struct GALAXY *g)
{
  free(g->GridHistory);
  free(g->GridStellarMass);
  free(g->GridSFR);
  free(g->GridZ);
  free(g->GridCentralGalaxyMass);
  free(g->EjectedFraction);
  free(g->LenHistory);
  free(g->Stars);
  free(g->GridOutflowRate);
  free(g->GridInfallRate);
  free(g->GridEjectedMass);
  free(g->QuasarActivity);
  free(g->DynamicalTime);
  free(g->QuasarSubstep);
  free(g->GridColdGas);
  free(g->LenMergerGal);
  free(g->GridBHMass);
  free(g->GridReionMod);
  free(g->GridDustColdGas);
  free(g->GridDustHotGas);
  free(g->GridDustEjectedMass);

  g->IsMalloced = 0;
}

int32_t malloc_grid_arrays(struct GALAXY *g)
{

#define ALLOCATE_GRID_MEMORY(name, length) \
{                                          \
  name = calloc(length, sizeof(*(name)));  \
  if (name == NULL)                        \
  {                                        \
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate"#name".\n", sizeof(*(name)* length)); \
    return EXIT_FAILURE;                   \
  }                                        \
}

  ALLOCATE_GRID_MEMORY(g->GridHistory, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->GridStellarMass, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->GridSFR, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->GridZ, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->GridCentralGalaxyMass, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->EjectedFraction, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->LenHistory, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->Stars, SN_Array_Len);
  ALLOCATE_GRID_MEMORY(g->GridOutflowRate, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->GridInfallRate, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->GridEjectedMass, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->QuasarActivity, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->DynamicalTime, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->QuasarSubstep, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->GridColdGas, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->LenMergerGal, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->GridBHMass, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->GridReionMod, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->GridDustColdGas, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->GridDustHotGas, MAXSNAPS);
  ALLOCATE_GRID_MEMORY(g->GridDustEjectedMass, MAXSNAPS);

  g->IsMalloced = 1; // This way we can check that we're not freeing memory that hasn't been allocated.

  return EXIT_SUCCESS;

}

#undef ALLOCATE_GRID_MEMORY

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

int32_t free_reion_lists(int32_t filenr)
{

  int32_t SnapNum;

  if (ReionSnap == LowSnap)
  {
    return EXIT_SUCCESS;
  }

  for (SnapNum = 0; SnapNum < ReionList->NumLists; ++SnapNum)
  {
    if (ReionList->ReionMod_List[SnapNum].NHalos_Ionized == 0) // No lists were read for this snapshot to move along.
    {
      continue;
    }

    if (ReionList->ReionMod_List[SnapNum].NHalos_Found != ReionList->ReionMod_List[SnapNum].NHalos_Ionized)
    {
      fprintf(stderr, "After processing file %d we only matched %d Halos to the reionization list.  The list contained %d Halos; these Halos MUST be in the tree file.\n", filenr, ReionList->ReionMod_List[SnapNum].NHalos_Found, ReionList->ReionMod_List[SnapNum].NHalos_Ionized);
      return EXIT_FAILURE;
    }

    free(ReionList->ReionMod_List[SnapNum].ReionMod);
    free(ReionList->ReionMod_List[SnapNum].HaloID);
  }

  free(ReionList->ReionMod_List);
  free(ReionList);

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
