#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>

#include "core_allvars_grid.h"
#include "core_proto_grid.h"

// keep a static file handle to remove the need to do constant seeking
FILE* load_fd = NULL;

void load_halos(int filenr)
{
  char buf[1000];

  // open the file each time this function is called
  sprintf(buf, "%s/%s.%d", SimulationDir, TreeName, filenr);
  if(!(load_fd = fopen(buf, "r")))
  {
    printf("can't open file `%s'\n", buf);
    exit(0); 
  }

#ifdef DEBUG_HALOS
  printf("Loading Halos from file %s\n", buf);
#endif
  fread(&Ntrees, 1, sizeof(int), load_fd); // Necessary to seek to the correct line.
  fread(&totNHalos, 1, sizeof(int), load_fd);

  TreeNHalos = malloc(sizeof(int) * Ntrees); // Seeking again.
  fread(TreeNHalos, Ntrees, sizeof(int), load_fd); // Seeking again.

  Halo = mymalloc(sizeof(struct halo_data) * totNHalos);
  fread(Halo, totNHalos, sizeof(struct halo_data), load_fd);

  fclose(load_fd);
//  if (Verbose == 1)
//    printf("Read in a total of %d Halos for file %d.\n", totNHalos, filenr);
}

void load_gals(char *fname)
{

  int i;

  if(!(load_fd = fopen(fname, "r")))
  {
    printf("can't open file `%s'\n", fname);
    exit(0);
  }

#ifdef DEBUG_GALS 
  printf("Loading Galaxies from file %s\n", fname);
#endif
  fread(&Ntrees, 1, sizeof(int), load_fd);
  fread(&NtotGals, sizeof(int), 1, load_fd);

  GalsForTree = malloc(Ntrees * sizeof(int));
  fread(GalsForTree, Ntrees, sizeof(int), load_fd);

  Gal = mymalloc(NtotGals * sizeof(struct GALAXY_INPUT));
  GalGrid = mymalloc(NtotGals * sizeof(struct GALAXY_GRID));
  //estimate_gal_memory(NtotGals);

  for (i = 0; i < NtotGals; ++i)
  { 
    fread(Gal, sizeof(struct GALAXY_INPUT), 1, load_fd);
    
    GalGrid[i].History = malloc(sizeof(int) * MAXSNAPS);
    fread(GalGrid[i].History, sizeof(int), MAXSNAPS, load_fd);

    GalGrid[i].StellarMass = malloc(sizeof(float) * MAXSNAPS);
    fread(GalGrid[i].StellarMass, sizeof(float), MAXSNAPS, load_fd);
    
    GalGrid[i].SFR = malloc(sizeof(float) * MAXSNAPS); 
    fread(GalGrid[i].SFR, sizeof(float), MAXSNAPS, load_fd);

    GalGrid[i].Z = malloc(sizeof(float) * MAXSNAPS);
    fread(GalGrid[i].Z, sizeof(float), MAXSNAPS, load_fd);

    GalGrid[i].CentralGalaxyMass = malloc(sizeof(float) * MAXSNAPS);
    fread(GalGrid[i].CentralGalaxyMass, sizeof(float), MAXSNAPS, load_fd);

    GalGrid[i].Photons_HI = malloc(sizeof(float) * MAXSNAPS);
    fread(GalGrid[i].Photons_HI, sizeof(float), MAXSNAPS, load_fd);

    GalGrid[i].Photons_HeI = malloc(sizeof(float) * MAXSNAPS);
    fread(GalGrid[i].Photons_HeI, sizeof(float), MAXSNAPS, load_fd);

    GalGrid[i].Photons_HeII = malloc(sizeof(float) * MAXSNAPS);
    fread(GalGrid[i].Photons_HeII, sizeof(float), MAXSNAPS, load_fd);

  }

  fclose(load_fd);
  free(GalsForTree);
}

void free_gals(void)
{

    int i;
    for (i = 0; i < NtotGals; ++i)
    { 
      free(GalGrid[i].History);
      free(GalGrid[i].StellarMass);
      free(GalGrid[i].SFR);
      free(GalGrid[i].Z);
      free(GalGrid[i].CentralGalaxyMass);
      free(GalGrid[i].Photons_HI);
      free(GalGrid[i].Photons_HeI);
      free(GalGrid[i].Photons_HeII);
    }

    myfree(GalGrid);
    myfree(Gal);

}

void load_merged_gals(char *fname)
{

  if(!(load_fd = fopen(fname, "r")))
  {
    printf("can't open file `%s'\n", fname);
    exit(0);
  }

  fread(&NtotGals, 1, sizeof(int), load_fd);
  Gal = malloc(NtotGals * sizeof(struct GALAXY_INPUT));
  fread(Gal, NtotGals, sizeof(struct GALAXY_INPUT), load_fd); 

  fclose(load_fd);

}
