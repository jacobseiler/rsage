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
  char buf[1000];
  FILE *fd;

	// open the file each time this function is called
  sprintf(buf, "%s/%s.%d", SimulationDir, TreeName, filenr);
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
  
  MaxGals = 100000;

  MaxMergedGals = MaxGals;
  FoF_MaxGals = 100000;

  HaloAux = mymalloc(sizeof(struct halo_aux_data) * TreeNHalos[nr]);
  HaloGal = mymalloc(sizeof(struct GALAXY) * MaxGals);
  Gal = mymalloc(sizeof(struct GALAXY) * FoF_MaxGals);
  MergedGal = mymalloc(sizeof(struct GALAXY) * MaxMergedGals);   
 
  for(i = 0; i < TreeNHalos[nr]; i++)
  {
    HaloAux[i].DoneFlag = 0;
    HaloAux[i].HaloFlag = 0;
    HaloAux[i].NGalaxies = 0;
    if (Halo[i].SnapNum == 99)
      count += Halo[i].Mvir; 
  }

}



void free_galaxies_and_tree(void)
{
  int i;

  for(i = 0; i < NumGals; ++i)
  { 
    if(HaloGal[i].SnapNum == ListOutputSnaps[0]) // The pointed to memory that we malloced is NOT copied over when we generate a new MergedGal entry. 
      free_grid_arrays(&HaloGal[i]);   
  }

  for(i = 0; i < MergedNr; ++i)
  {
    free_grid_arrays(&MergedGal[i]); // These are the Gal[xxx] entries that were copied over to 
  } 

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
}

void malloc_grid_arrays(struct GALAXY *g)
{
  if(NULL == (g->GridHistory = malloc(sizeof(*(g->GridHistory)) * MAXSNAPS)))
  {
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridHistory.", sizeof(*(g->GridHistory))*MAXSNAPS); 
    exit(EXIT_FAILURE);
  }

  if(NULL == (g->GridStellarMass = malloc(sizeof(*(g->GridStellarMass)) * MAXSNAPS)))
  {
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridStellarMass.", sizeof(*(g->GridStellarMass))*MAXSNAPS); 
    exit(EXIT_FAILURE);
  } 

  if(NULL == (g->GridSFR = malloc(sizeof(*(g->GridSFR)) * MAXSNAPS)))
  {
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridSFR.", sizeof(*(g->GridSFR))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  if (NULL == (g->GridZ = malloc(sizeof(*(g->GridZ)) * MAXSNAPS)))
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridSFR.", sizeof(*(g->GridZ))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }
 
  if (NULL == (g->GridCentralGalaxyMass = malloc(sizeof(*(g->GridCentralGalaxyMass)) * MAXSNAPS)))
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridCentralGalaxyMass.", sizeof(*(g->GridCentralGalaxyMass))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  if (NULL == (g->MfiltGnedin = malloc(sizeof(*(g->MfiltGnedin)) * MAXSNAPS)))
  {
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate MfiltGnedin.", sizeof(*(g->MfiltGnedin))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  if (NULL == (g->MfiltSobacchi = malloc(sizeof(*(g->MfiltSobacchi)) * MAXSNAPS)))
  {   
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate MfiltSobacchi.", sizeof(*(g->MfiltSobacchi))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }
 
  if (NULL == (g->EjectedFraction = malloc(sizeof(*(g->EjectedFraction)) * MAXSNAPS)))
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate EjectedFraction.", sizeof(*(g->EjectedFraction))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }
 
  if (NULL == (g->LenHistory = malloc(sizeof(*(g->LenHistory)) * MAXSNAPS)))
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate LenHistory.", sizeof(*(g->LenHistory))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }

  if (NULL == (g->Stars = malloc(sizeof(*(g->Stars)) * SN_Array_Len)))
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate Stars.", sizeof(*(g->Stars))*SN_Array_Len);
    exit(EXIT_FAILURE);
  }

  if (NULL == (g->GridOutflowRate = malloc(sizeof(*(g->GridOutflowRate)) * MAXSNAPS)))
  { 
    fprintf(stderr, "Out of memory allocating %ld bytes, could not allocate GridOutflowRate.", sizeof(*(g->GridOutflowRate))*MAXSNAPS);
    exit(EXIT_FAILURE);
  }


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
