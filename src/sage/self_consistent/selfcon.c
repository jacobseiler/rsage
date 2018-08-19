/*
This file contains functions that control the reionization feedback for the self-consistent RSAGE model.

In the RSAGE model, all halos with a `ReionizationModifier` less than 1 are saved in a file with a unique 64 bit ID (left-most 32 bits are the `treenr` and the right-most 32 bits are the tree-local `HaloID`).
Whenever `ReionizationModifier` is needed we search through this file to see if the specified halo has a reionization modifier less than 1. 

Since we need to remember the history of `ReionizationModifier`, we load in all the `ReionizationModifier` lists for the current snapshot iteration in addition to all previous snapshots. 

Author: Jacob Seiler
Version: 0.0.1
*/
 
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

/*
Initializes the array (`ReionMod_List`) to hold the halos with `ReionizationModifier` less than 1.  
Also reads the file that holds the `ReionizationModifier` values and their unique 64 bit ID.

Parameters
----------

filenr: Integer. 
  The current tree file number. 

Returns
----------

EXIT_SUCCESS or EXIT_FAILURE.
  If this is the first iteration of the self-consistent model (as RSAGE is run for every snapshot), EXIT_SUCCESS is returned (no reionization has occurred). 

  If the `ReionizationModifier` file cannot be opened, EXIT_FAILURE is returned. 
  If memory cannot be allocated for the outer struct (REIONMOD_STRUCT), EXIT_FAILURE is returned. 
  If memory cannot be allocated for the inner struct (REIONMOD_LIST), EXIT_FAILURE is returned. 
  If memory cannot be allocated for the `HaloID` and `ReionMod` arrays for each snapshot, EXIT_FAILURE is returned.

  The data structure of the `ReionizationModifier` file is:
  <Snapshot Number (int32_t)> 
  <Number of Reionized Halos in this Snapshot (int32_t, NHalos)>
  <NHalos HaloIDs (int64_t)>
  <NHalos ReionizationModifier Values (32 bit float)>
  
  <...Next Snapshot Number...> 

  If this data structure is not correct, EXIT_FAILURE is returned.

  Otherwise EXIT_SUCCESS is returned.

Pointer Updates
----------

`ReionList` is allocated memory.
`ReionList.ReionMod_List` is allocated enough memoery to hold lists for the current snapshot and all previous ones.
For each snapshot, `ReionList->ReionMod_List.HaloID` and `ReionList->ReionMod_List.ReionMod` are allocated memory to hold all of the `HaloIDs` and `ReionizationModifier` values for all halos in ionized regions. 

Units  
----------

The `HaloID` variable is unique and defined as `int64_t ((int32_t) treenr << 32) | (int32_t) ForestSpecificHaloID`.
The `ReionizationModifier` variable is unitless and specifies the suppression of baryonic infall. 
*/

int32_t init_reion_lists(int32_t filenr)
{

  FILE *ListFile;
  char ListFile_name[MAXLEN];
  int32_t SnapNum, SnapNum_Read;
 
  if (ReionSnap == LowSnap) // This is the first iteration of the self-consistent run so there will be no ionization yet.
  {
    return EXIT_SUCCESS;
  }

  snprintf(ListFile_name, MAXLEN, "%s/reionization_modifiers/%s_treefile_%03d", PhotoionDir, FileNameGalaxies, filenr);    

  ListFile = fopen(ListFile_name, "rb");
  if (ListFile == NULL)
  {
    fprintf(stderr, "Cannot open file %s\n", ListFile_name);
    return EXIT_FAILURE;
  }

  ReionList = malloc(sizeof(struct REIONMOD_STRUCT));
  if (ReionList == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for Reionization List struct\n");
    return EXIT_FAILURE;
  }

  ReionList->NumLists = ReionSnap; 

  ReionList->ReionMod_List = malloc(sizeof(struct REIONMOD_LIST) * ReionList->NumLists);
  if (ReionList == NULL)
  {
    fprintf(stderr, "Cannot allocate memory for the ReionMod_List struct\n");
    return EXIT_FAILURE;
  }

  for (SnapNum = 0; SnapNum < ReionList->NumLists; ++SnapNum)
  {

    fread(&SnapNum_Read, sizeof(SnapNum_Read), 1, ListFile);

    fread(&ReionList->ReionMod_List[SnapNum].NHalos_Ionized, sizeof(SnapNum_Read), 1, ListFile);
    //printf("Snapshot %d has %d Halos in the list.\n", SnapNum_Read, ReionList->ReionMod_List[SnapNum].NHalos_Ionized);

    if (SnapNum_Read != SnapNum)
    { 
      fprintf(stderr, "When attempting to read the reionization modifier lists, the read file had a snapshot number %d when we expected a number %d\n", SnapNum_Read, SnapNum);
      fprintf(stderr, "Check that `filter_mass` was properly compiled with MPI.\n");
      return EXIT_FAILURE;
    }

    ReionList->ReionMod_List[SnapNum].NHalos_Found = 0;
 
    if (ReionList->ReionMod_List[SnapNum].NHalos_Ionized == 0) // There were no halos within ionized regions for this snapshot, reionization hasn't started or is in the beginning.
    {
      continue;
    }
   
    ReionList->ReionMod_List[SnapNum].HaloID = mycalloc(ReionList->ReionMod_List[SnapNum].NHalos_Ionized, sizeof(*(ReionList->ReionMod_List[SnapNum].HaloID))); 
    if (ReionList->ReionMod_List[SnapNum].HaloID == NULL)
    {
      fprintf(stderr, "Cannot allocate memoery for ReionMod_List.HaloID for Snapshot Number %d\n", SnapNum);
      return EXIT_FAILURE;
    }

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
    if (ReionList->ReionMod_List[SnapNum].ReionMod == NULL)
    {
      fprintf(stderr, "Cannot allocate memoery for ReionMod_List.ReionMod for Snapshot Number %d\n", SnapNum);
      return EXIT_FAILURE;
    }

    fread(ReionList->ReionMod_List[SnapNum].ReionMod, sizeof(*(ReionList->ReionMod_List[SnapNum].ReionMod)), ReionList->ReionMod_List[SnapNum].NHalos_Ionized, ListFile);      
  }
  fclose(ListFile);
 
  return EXIT_SUCCESS;
}

/*
Frees the unique `HaloID` and `ReionizationModifier` arrays. 

Parameters
----------

filenr: Integer. 
  The current tree file number. 

Returns
----------

EXIT_SUCCESS or EXIT_FAILURE.
  If this is the first iteration of the self-consistent model (as RSAGE is run for every snapshot), EXIT_SUCCESS is returned (no reionization has occurred). 

  As this function is called after all galaxies have been fully evolved, all halos in the `ReionMod_List` struct MUST have been found in the main code.
  If there were halos that were not found in the main code, EXIT_FAILURE is returned. 
   
  Otherwise EXIT_SUCCESS is returned.

Pointer Updates
----------

For each snapshot, `ReionList->ReionMod_List.HaloID` and `ReionList->ReionMod_List.ReionMod` are freed. 
`ReionList.ReionMod_List` is freed for the current snapshot and all previous ones.
`ReionList` is freed.

Units  
----------

None.
*/

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

/*
For a given halo, determines if the halo requires checking for a non-one `ReionizationModifier` value.  As RSAGE is iterated over snapshots, halos at snapshots later than the one being processed are ignored.
E.g., if we are processing Snapshot 40, a halo at Snapshot 45 will have a `ReionizationModifier` value of 1.

If the halo needs to be checked, the `search_for_modifier()` function is called to check the `ReionList` struct for a matching unique 64 bit `HaloID`.
If the halo is found in the list, the corresponding value for `ReionizationModifier` is used, otherwise a value of 1 is used. 

**NOTE** If ``ReionizationOn`` = 4, we do not use the `ReionizationModifier` value in the list.  Instead we use the value from the Gnedin Analytic formulat (see `do_reionization()` in `model_infall.c`).

Parameters
----------

gal: Integer.
  The local galaxy index.  Used to reference the galaxy property arrays.

halonr: Integer. 
  The background FoF Halo. 
  **IMPORTANT**: This is the BACKGROUND FOF HALO. This has important consequences for satellite galaxies/halos (see https://github.com/jacobseiler/self_consistent_SAGE/issues/6).

increment_counter: Integer.  
  This function is called multiple times for an individual halo.  Since we keep track of the number of halos successfully found in the list, we only want to increment this counter once.
  This flag denotes whether we have accounted for this halo being found yet.

*reionization modifier: Double pointer.
  Pointer that stores the value for the reionization modifier for this galaxy.

Returns
----------

EXIT_SUCCESS or EXIT_FAILURE.
  If this is the first iteration of the self-consistent model (as RSAGE is run for every snapshot), EXIT_SUCCESS is returned (no reionization has occurred). 
  If the halo is at a snapshot greater than the current snapshot being processed, EXIT_SUCCESS is returned. 
  If the halo is found in the `ReionMod_List` struct, we do a check to ensure all the indices line up.  If they don't, EXIT_FAILURE is returned.
  
  Otherwise EXIT_SUCCESS is returned.

Pointer Updates
----------

The `*reionization_modifier` pointer is updated with the suppressed value for this halo.

Units  
----------

None.
*/

int32_t do_self_consistent_reionization(int32_t gal, int32_t halonr, int32_t increment_counter, double *reionization_modifier)
{

  int32_t treenr; 
  int64_t HaloID, found_idx = -1;

  if (ReionSnap == LowSnap)
  {
    *reionization_modifier = 1.0;
    return EXIT_SUCCESS; 
  }

  if (Halo[halonr].SnapNum >= ReionSnap) // We have yet to do reionization for this snapshot.
  {
    *reionization_modifier = 1.0; // Reionization hasn't happened yet for this halo.
    return EXIT_SUCCESS; 
  } 
  
  if (ReionList->ReionMod_List[Halo[halonr].SnapNum].NHalos_Ionized == 0) // There are no halos within ionized regions for this snapshot.
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

/*
Given a unique HaloID, does a binary search for a matching HaloID in the `ReionMod_List` struct.
If the halo is found, uses the `ReionMod` value in the list, otherwise a value of 1 will be used.


Parameters
----------

match_HaloID: 64 bit Integer. 
  The unique HaloID we're searching for.

SnapNum: Integer. 
  Snapshot the halo is at.

increment_counter: Integer.  
  This function is called multiple times for an individual halo.  Since we keep track of the number of halos successfully found in the list, we only want to increment this counter once.
  This flag denotes whether we have accounted for this halo being found yet.

*found_idx: 64 bit Integer pointer. 
  Pointer that stores the index in the `ReionMod_List` array that corresponds to the found halo.
  If the halo is not in the list, this will be -1.  

Returns
----------

reionization_modifier: Double.
  If the Halo is found, the associated reionization_modifier in the list. 
  Otherwise, the halo is not in an ionized region and reionization_modifier = 1. 

Pointer Updates
----------

The `*found_idx` pointer is updated with the `ReionMod_List` index that corresponds to the found halo.
If the halo is not found, this will be -1.

Units  
----------

None.
*/

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
  }

  if (is_found == 1)
  {
    reionization_modifier = ReionList->ReionMod_List[SnapNum].ReionMod[search_idx];
    *found_idx = search_idx;
    if (increment_counter == 1) // Only want to increment our counters once as for satellite galaxies with a subhalo, strip_from_satellite is called multiple times. 
    {
      ++ReionList->ReionMod_List[SnapNum].NHalos_Found;
#ifdef DEBUG_SELFCON
      printf("Found unique HaloID %ld at index %ld\n", match_HaloID, *found_idx);
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
