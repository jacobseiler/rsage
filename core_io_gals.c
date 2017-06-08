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
  char buf[MAXLEN];

  // open the file each time this function is called
  snprintf(buf, MAXLEN, "%s/%s.%d", SimulationDir, TreeName, filenr);
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
    
    GalGrid[i].History = malloc(sizeof(*(GalGrid[i].History)) * MAXSNAPS);
    fread(GalGrid[i].History, sizeof(*(GalGrid[i].History)), MAXSNAPS, load_fd);

    GalGrid[i].StellarMass = malloc(sizeof(*(GalGrid[i].StellarMass)) * MAXSNAPS);
    fread(GalGrid[i].StellarMass, sizeof(*(GalGrid[i].StellarMass)), MAXSNAPS, load_fd);
    
    GalGrid[i].SFR = malloc(sizeof(*(GalGrid[i].SFR)) * MAXSNAPS); 
    fread(GalGrid[i].SFR, sizeof(*(GalGrid[i].SFR)), MAXSNAPS, load_fd);

    GalGrid[i].Z = malloc(sizeof(*(GalGrid[i].Z)) * MAXSNAPS);
    fread(GalGrid[i].Z, sizeof(*(GalGrid[i].Z)), MAXSNAPS, load_fd);

    GalGrid[i].CentralGalaxyMass = malloc(sizeof(*(GalGrid[i].CentralGalaxyMass)) * MAXSNAPS);
    fread(GalGrid[i].CentralGalaxyMass, sizeof(*(GalGrid[i].CentralGalaxyMass)), MAXSNAPS, load_fd);

    GalGrid[i].Pad = malloc(sizeof(*(GalGrid[i].Pad)) * MAXSNAPS);
    fread(GalGrid[i].Pad, sizeof(*(GalGrid[i].Pad)), MAXSNAPS, load_fd);

    GalGrid[i].MfiltGnedin = malloc(sizeof(*(GalGrid[i].MfiltGnedin)) * MAXSNAPS);
    fread(GalGrid[i].MfiltGnedin, sizeof(*(GalGrid[i].MfiltGnedin)), MAXSNAPS, load_fd);

    GalGrid[i].MfiltSobacchi = malloc(sizeof(*(GalGrid[i].MfiltSobacchi)) * MAXSNAPS);
    fread(GalGrid[i].MfiltSobacchi, sizeof(*(GalGrid[i].MfiltSobacchi)), MAXSNAPS, load_fd);
 
    GalGrid[i].EjectedFraction = malloc(sizeof(*(GalGrid[i].EjectedFraction)) * MAXSNAPS);
    fread(GalGrid[i].EjectedFraction, sizeof(*(GalGrid[i].EjectedFraction)), MAXSNAPS, load_fd);

    GalGrid[i].LenHistory = malloc(sizeof(*(GalGrid[i].LenHistory)) * MAXSNAPS);
    fread(GalGrid[i].LenHistory, sizeof(*(GalGrid[i].LenHistory)), MAXSNAPS, load_fd);

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
      free(GalGrid[i].MfiltGnedin);
      free(GalGrid[i].MfiltSobacchi); 
      free(GalGrid[i].EjectedFraction);
      free(GalGrid[i].LenHistory);
      free(GalGrid[i].Pad);
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
