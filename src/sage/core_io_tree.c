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
#include "temporal_array.h"
#include "io/tree_binary.h"

// keep a static file handle to remove the need to do constant seeking

void load_tree_table(int filenr, int32_t treestyle)
{
  int32_t i;
  char buf[MAXLEN];
   
  if (treestyle == 0)
    snprintf(buf, MAXLEN, "%s/%s_%03d%s", SimulationDir, TreeName, filenr, TreeExtension);
  else
    snprintf(buf, MAXLEN, "%s/%s.%d%s", SimulationDir, TreeName, filenr, TreeExtension);

  load_tree_table_binary(buf); // Gets the total number of trees, halos and halos per tree.

  TreeNMergedgals = mycalloc(Ntrees, sizeof(*(TreeNMergedgals)));
  if (TreeNMergedgals == NULL)
  {
    fprintf(stderr, "Could not allocate memory for `TreeNMergedgals`.\n");
    ABORT(EXIT_FAILURE);
  }

  TreeNgals = mycalloc(Ntrees, sizeof(*(TreeNgals))); 
  if (TreeNgals == NULL)
  {
    fprintf(stderr, "Could not allocate memory for `TreeNgals`.\n");
    ABORT(EXIT_FAILURE);
  }
   
  for(i = 0; i < Ntrees; i++)
  {
    TreeNgals[i] = 0;
    TreeNMergedgals[i] = 0;
  }

  TotGalaxies = 0;  
  TotMerged = 0;

}



void free_tree_table(void)
{

  myfree(TreeNMergedgals, sizeof(*(TreeNMergedgals)) * Ntrees); 

  myfree(TreeNgals, sizeof(*(TreeNgals)) * Ntrees);

  myfree(TreeFirstHalo, sizeof(*(TreeFirstHalo)) * Ntrees);
  myfree(TreeNHalos, sizeof(*(TreeNHalos)) * Ntrees);
	
	// Don't forget to free the open file handle

  close_binary_file();

}

void load_tree(int nr)
{
  int i;

  // must have an FD
  
  load_tree_binary(nr);

  MaxGals = (int)(MAXGALFAC * TreeNHalos[nr]);

  if(MaxGals < 10000)  
    MaxGals = 10000;

  MaxMergedGals = MaxGals;
  FoF_MaxGals = 10000; 

  // Previously used mycalloc instead of malloc but there was immense slowdowns (~3x).
  // I profiled this using Vtune and it was the implementation of the function itself that was causing the slowdown.
  // This may be because this loop is hit so often that initializing all the memory to 0 took a long time.

  gal_to_free = mymalloc(MaxMergedGals * sizeof(int));
  if (gal_to_free == NULL)
  {
    fprintf(stderr, "Could not allocate memory for `gal_to_free`.\n");
    ABORT(EXIT_FAILURE);
  }

  HaloAux = mymalloc(TreeNHalos[nr] * sizeof(struct halo_aux_data));
  if (HaloAux == NULL)
  {
    fprintf(stderr, "Could not allocate memory for `HaloAux`.\n");
    ABORT(EXIT_FAILURE);
  }
 
  HaloGal = mymalloc(MaxGals * sizeof(struct GALAXY));
  if (HaloGal == NULL)
  {
    fprintf(stderr, "Could not allocate memory for `HaloGal`.\n");
    ABORT(EXIT_FAILURE);
  }

  Gal = mymalloc(FoF_MaxGals * sizeof(struct GALAXY));
  if (Gal == NULL)
  {
    fprintf(stderr, "Could not allocate memory for `Gal`.\n");
    ABORT(EXIT_FAILURE);
  }

  MergedGal = mymalloc(MaxMergedGals * sizeof(struct GALAXY));
  if (MergedGal == NULL)
  {
    fprintf(stderr, "Could not allocate memory for `MergedGal`.\n");
    ABORT(EXIT_FAILURE);
  }

  double Min_Halo = 1e5;
  double Max_Halo = 0.0; 
  for(i = 0; i < TreeNHalos[nr]; i++)
  {

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

    free_temporal_arrays(&HaloGal[gal_to_free[i]]);
    ++gal_frees; 
  }
  myfree(gal_to_free, sizeof(*(gal_to_free)) * MaxMergedGals);


  // Now we just have to free the arrays for the galaxies that have merged.

  for(i = 0; i < MergedNr; ++i)
  {
      XPRINT(MergedGal[i].IsMalloced == 1, "MergedGal %d doesn't have grids mallocced but we're trying to free it.\n", i); 
    free_temporal_arrays(&MergedGal[i]); // These are the Gal[xxx] entries that were copied over to  
    ++mergedgal_frees;
  } 

  // All the inside pointers have now been freed, lets free the structs themselves now.
  myfree(MergedGal, sizeof(*(MergedGal)) * MaxMergedGals);
  myfree(Gal, sizeof(*(Gal)) * FoF_MaxGals);
  myfree(HaloGal, sizeof(*(HaloGal)) * MaxGals);
  myfree(HaloAux, sizeof(*(HaloAux)) * TreeNHalos[treenr]);
  myfree(Halo, sizeof(*(Halo)) * TreeNHalos[treenr]);

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
