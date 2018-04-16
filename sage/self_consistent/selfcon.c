#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>

#include "../core_allvars.h"
#include "../core_proto.h"
#include "selfcon.h"

// Local Variables //

// Local Proto-Types //

// External Functions //

int32_t init_reion_lists(int32_t filenr)
{

  FILE *ListFile;
  char ListFile_name[MAXLEN];
  int32_t SnapNum, SnapNum_Read;
 
  if (ReionSnap == LowSnap) // This is the first iteration of the self-consistent run so there will be no ionization yet.
  {
    return EXIT_SUCCESS;
  }

  snprintf(ListFile_name, MAXLEN, "%s/reionization_modifiers/treefile_%03d", PhotoionDir, filenr);    

  ListFile = fopen(ListFile_name, "rb");
  if (ListFile == NULL)
  {
    fprintf(stderr, "Cannot open file %s\n", ListFile_name);
    return EXIT_FAILURE;
  }
 
  printf("Reading in the reionization modifier lists.\n");
 
  ReionList = malloc(sizeof(struct REIONMOD_STRUCT));
  if (ReionList == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for Reionization List struct\n");
    return EXIT_FAILURE;
  }

  ReionList->NumLists = ReionSnap; 

  ReionList->ReionMod_List = malloc(sizeof(struct REIONMOD_LIST) * ReionList->NumLists);

  for (SnapNum = 0; SnapNum < ReionList->NumLists; ++SnapNum)
  {

    fread(&SnapNum_Read, sizeof(int32_t), 1, ListFile);

    fread(&ReionList->ReionMod_List[SnapNum].NHalos_Ionized, sizeof(int32_t), 1, ListFile);
    //printf("Snapshot %d has %d Halos in the list.\n", SnapNum_Read, ReionList->ReionMod_List[SnapNum].NHalos_Ionized);

    if (SnapNum_Read != SnapNum)
    { 
      fprintf(stderr, "When attempting to read the reionization modifier lists, the read file had a snapshot number %d when we expected a number %d\n", SnapNum_Read, SnapNum);
      return EXIT_FAILURE;
    }

    ReionList->ReionMod_List[SnapNum].NHalos_Found = 0;
 
    if (ReionList->ReionMod_List[SnapNum].NHalos_Ionized == 0) // There were no halos within ionized regions for this snapshot, reionization hasn't started or is in the beginning.
    {
      continue;
    }
   
    ReionList->ReionMod_List[SnapNum].HaloID = mycalloc(ReionList->ReionMod_List[SnapNum].NHalos_Ionized, sizeof(*(ReionList->ReionMod_List[SnapNum].HaloID))); 

    fread(ReionList->ReionMod_List[SnapNum].HaloID, sizeof(*(ReionList->ReionMod_List[SnapNum].HaloID)), ReionList->ReionMod_List[SnapNum].NHalos_Ionized, ListFile);

#ifdef DEBUG_SELFCON    
    int32_t i;
    for (i = 0; i < ReionList->ReionMod_List[SnapNum].NHalos_Ionized; ++i)
    {
      int64_t ID;
      ID = ReionList->ReionMod_List[SnapNum].HaloID[i];
      printf("File %d: HaloID %ld is in the list, corresponding to tree %d and Halo number %d\n", filenr, ID, (int32_t)(ID >> 32), (int32_t)ID); 
    }
#endif

    ReionList->ReionMod_List[SnapNum].ReionMod = mycalloc(ReionList->ReionMod_List[SnapNum].NHalos_Ionized, sizeof(*(ReionList->ReionMod_List[SnapNum].ReionMod))); 

    fread(ReionList->ReionMod_List[SnapNum].ReionMod, sizeof(*(ReionList->ReionMod_List[SnapNum].ReionMod)), ReionList->ReionMod_List[SnapNum].NHalos_Ionized, ListFile);
  
    
  }
  fclose(ListFile);
 
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

    myfree(ReionList->ReionMod_List[SnapNum].ReionMod, sizeof(*(ReionList->ReionMod_List[SnapNum].ReionMod)) * ReionList->ReionMod_List[SnapNum].NHalos_Ionized);
    myfree(ReionList->ReionMod_List[SnapNum].HaloID, sizeof(*(ReionList->ReionMod_List[SnapNum].HaloID)) * ReionList->ReionMod_List[SnapNum].NHalos_Ionized);
  }

  free(ReionList->ReionMod_List);
  free(ReionList);

  return EXIT_SUCCESS;

}

int32_t do_self_consistent_reionization(int32_t gal, int32_t halonr, int32_t increment_counter, double *reionization_modifier)
{

  int32_t treenr; 
  int64_t HaloID, found_idx = -1;

  if (ReionSnap == LowSnap)
  {
    *reionization_modifier = 1.0;
    return EXIT_SUCCESS; 
  }

  if ((Halo[halonr].SnapNum > ReionSnap) || (ReionList->ReionMod_List[Halo[halonr].SnapNum].NHalos_Ionized == 0)) // We have yet to do reionization for this snapshot or if there are no halos within ionized regions for this snapshot.
  {
    *reionization_modifier = 1.0; // Reionization hasn't happened yet for this halo.
    return EXIT_SUCCESS; 
  } 
  
  ReionList->NumLists = ReionSnap; 

  treenr = Gal[gal].TreeNr;

  HaloID = ((int64_t)treenr << 32) | (int64_t)halonr; // Generates the unique ID for each halo within this file. 

  *reionization_modifier = search_for_modifier(HaloID, Halo[halonr].SnapNum, increment_counter, &found_idx); 

  if (found_idx == -1)
  {
    *reionization_modifier = 1.0;
    return EXIT_SUCCESS;
  }

  // If ReionizationOn is 3 we want to do the proper self-consistent treatment.
  // However if ReionizationOn is 4 we want to use the analytic formula to determine reionization_modifier for those halos in ionized regions. 
  if (ReionizationOn == 4)
  {
    *reionization_modifier = do_reionization(gal, ZZ[Halo[halonr].SnapNum], 0);
  }

  if (ReionList->ReionMod_List[Halo[halonr].SnapNum].HaloID[found_idx] != HaloID)
  {
    fprintf(stderr, "After searching for HaloID %ld (corresponding to tree number %d and halo number %d) within the reionization list, the found_idx was %ld.  However queuring the reionization list at snapshot %d, the HaloID at %ld index is %ld\n", (long)HaloID, treenr, halonr, (long)found_idx, Halo[halonr].SnapNum, (long)found_idx, ReionList->ReionMod_List[Halo[halonr].SnapNum].HaloID[found_idx]);

    return EXIT_FAILURE; 
  } 

 
  return EXIT_SUCCESS; 

}

// Local Functions //

double search_for_modifier(int64_t match_HaloID, int32_t SnapNum, int32_t increment_counter, int64_t *found_idx)
{

  int32_t is_found;
  int64_t count, search_idx, search_HaloID;
  int32_t number_search_IDs;
  double reionization_modifier;

  is_found = 0;
  count = 0;
  number_search_IDs = ReionList->ReionMod_List[SnapNum].NHalos_Ionized; // This is the numbers of IDs that we are searching through. 

  search_idx = ceil(number_search_IDs / 2.0) - 1; 

  if (match_HaloID == 38654705745)
  {
    printf("search_idx = %ld\t number_search_IDs = %d\n", (long)search_idx, number_search_IDs);
  }


  while (is_found == 0)
  {
    ++count;
    search_HaloID = ReionList->ReionMod_List[SnapNum].HaloID[search_idx];

    if (match_HaloID == search_HaloID) 
    {
      is_found = 1;
    }
    else if (number_search_IDs / pow(2, count) < 1.0) // The smallest index movement is less than 1.  The HaloID isn't in the list. 
    {
      break;
    }
    else if (match_HaloID > search_HaloID) // The HaloID we are trying to match is larger than the ID we're currently comparing against. 
    {
      search_idx = search_idx + ceil(number_search_IDs / pow(2, count + 1)); // Move up the list.
    }
    else // Otherwise the HaloID we are trying to match is smaller than the ID we're currently comparing against.
    {
      search_idx = search_idx - ceil(number_search_IDs / pow(2, count + 1)); // Move down the list. 
    }

    // When the HaloID we're trying to match is at idx 0 it can be tough to
    // find.  For example, if we're at idx 1 then we try and move down the
    // list, we will end up in negatives!  In this case, fix up the edge case.
    if (search_idx < 0) 
    {
      search_idx = 0;
    }

    if (search_idx >= number_search_IDs) // Fix edge case.
    {
      search_idx = number_search_IDs -1;
    }


    if (match_HaloID == 38654705745)
    {
      printf("search_idx = %ld\t search_HaloID = %ld\tcount = %ld\n", (long)search_idx, (long)search_HaloID, (long)count);
    }

  }

  if (is_found == 1)
  {
    *found_idx = search_idx;
    if (increment_counter == 1) // Only want to increment our counters once as for satellite galaxies with a subhalo, strip_from_satellite is called multiple times. 
    {
      ++ReionList->ReionMod_List[SnapNum].NHalos_Found;
#ifdef DEBUG_SELFCON
      printf("Filenr %d: Found unique HaloID %ld at index %ld with modifier %.4f\n", filenr, match_HaloID, *found_idx, reionization_modifier);
#endif 
    }
  }
  else
  {
    reionization_modifier = 1.0;
    *found_idx = -1;
  }

  return reionization_modifier;

}



