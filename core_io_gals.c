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

void free_meraxes_halos(void)
{
  myfree(meraxes_Halo);
}

int load_meraxes_halos(int snapnum)
{

  char buf[MAXLEN];
  char tag[MAXLEN];
  char line[MAXLEN];

  int num_halos = 0;

  XASSERT(snapnum < 1000, "Help, we have more than 1000 snapshots!  The code was written for only three digits.\n");

  double SFR;
  double StellarMass;
  double Mvir;
  double ColdGas;
  double HotGas;
  double EjectedGas;
  double MetalsColdGas;
  double Pos[3]; 

  if(snapnum < 10)
  {
    snprintf(tag, MAXLEN, "00%d", snapnum);
  }
  else if (snapnum < 100)
  {	
    snprintf(tag, MAXLEN, "0%d", snapnum);
  }
  else
  {
    snprintf(tag, MAXLEN, "%d", snapnum);
  }

  snprintf(buf, MAXLEN, "%s/halolist_%s.Sfr.txt", GalaxiesInputDir, tag); 

  if(!(load_fd = fopen(buf, "r")))
  {
    printf("can't open file `%s'\n", buf);
    exit(0); 
  }

  while (fgets(line, MAXLEN, load_fd))
  {
    if(*line == '#')
      continue;    
    num_halos++;
  }
  fclose(load_fd);
  fprintf(stderr, "For snapshot %d there is %d halos.\n", snapnum, num_halos);  
  meraxes_Halo = mymalloc(sizeof(struct meraxes_halo_data) * (num_halos)); 

  int counter = 0; 

  // SFR //

  snprintf(buf, MAXLEN, "%s/halolist_%s.Sfr.txt", GalaxiesInputDir, tag);  
  if(!(load_fd = fopen(buf, "r")))
  {
    printf("can't open file `%s'\n", buf);
    exit(0); 
  }
 
  while (fgets(line, MAXLEN, load_fd))
  {
    if(*line == '#')
	continue; 
    sscanf(line, "%lf", &SFR);
    meraxes_Halo[counter].SFR = SFR;
    counter++;
  } 
  fclose(load_fd);
  fprintf(stderr, "Read in the SFR data.\n");
  // Stellar Mass //

  counter = 0;

  snprintf(buf, MAXLEN, "%s/halolist_%s.StellarMass.txt", GalaxiesInputDir, tag);  
  if(!(load_fd = fopen(buf, "r")))
  {
    printf("can't open file `%s'\n", buf);
    exit(0); 
  }
 
  while (fgets(line, MAXLEN, load_fd))
  {
    if(*line == '#')
	continue;
    sscanf(line, "%lf", &StellarMass);
    meraxes_Halo[counter].StellarMass = StellarMass;
    counter++;
  } 
  fclose(load_fd);

  fprintf(stderr, "Read in the StellarMass data.\n");
  // Mvir //  

  counter = 0;

  snprintf(buf, MAXLEN, "%s/halolist_%s.Mvir.txt", GalaxiesInputDir, tag);  
  if(!(load_fd = fopen(buf, "r")))
  {
    printf("can't open file `%s'\n", buf);
    exit(0); 
  }
 
  while (fgets(line, MAXLEN, load_fd))
  {
    if(*line == '#')
	continue;
    sscanf(line, "%lf", &Mvir);
    meraxes_Halo[counter].Mvir = Mvir;
    counter++;
  } 
  fclose(load_fd);

  fprintf(stderr, "Read in the Mvir data.\n");
  // ColdGas //  

  counter = 0;

  snprintf(buf, MAXLEN, "%s/halolist_%s.ColdGas.txt", GalaxiesInputDir, tag);  
  if(!(load_fd = fopen(buf, "r")))
  {
    printf("can't open file `%s'\n", buf);
    exit(0); 
  }
 
  while (fgets(line, MAXLEN, load_fd))
  {
    if(*line == '#')
	continue;
    sscanf(line, "%lf", &ColdGas);
    meraxes_Halo[counter].ColdGas = ColdGas;
    counter++;
  } 
  fclose(load_fd);

  fprintf(stderr, "Read in the ColdGas data.\n");
  // HotGas //  

  counter = 0;

  snprintf(buf, MAXLEN, "%s/halolist_%s.HotGas.txt", GalaxiesInputDir, tag);  
  if(!(load_fd = fopen(buf, "r")))
  {
    printf("can't open file `%s'\n", buf);
    exit(0); 
  }
 
  while (fgets(line, MAXLEN, load_fd))
  {
    if(*line == '#')
	continue;
    sscanf(line, "%lf", &HotGas);
    meraxes_Halo[counter].HotGas = HotGas;
    counter++;
  } 
  fclose(load_fd);

  fprintf(stderr, "Read in the HotGas data.\n");
  // EjectedGas // 

  counter = 0;

  snprintf(buf, MAXLEN, "%s/halolist_%s.EjectedGas.txt", GalaxiesInputDir, tag);  
  if(!(load_fd = fopen(buf, "r")))
  {
    printf("can't open file `%s'\n", buf);
    exit(0); 
  }
 
  while (fgets(line, MAXLEN, load_fd))
  {
    if(*line == '#')
	continue;
    sscanf(line, "%lf", &EjectedGas);
    meraxes_Halo[counter].EjectedGas = EjectedGas;
    counter++;
  } 
  fclose(load_fd);

  fprintf(stderr, "Read in the EjectedGas data.\n");
  // Pos // 

  counter = 0;

  snprintf(buf, MAXLEN, "%s/halolist_%s.Pos.txt", GalaxiesInputDir, tag);  
  if(!(load_fd = fopen(buf, "r")))
  {
    printf("can't open file `%s'\n", buf);
    exit(0); 
  }
 
  while (fgets(line, MAXLEN, load_fd))
  {
    if(*line == '#')
	continue;
    sscanf(line, "%lf %lf %lf", &Pos[0], &Pos[1], &Pos[2]);
    meraxes_Halo[counter].Pos[0] = Pos[0]; // Note to Future Jacob, if Manodeep ever says anything about putting this in a loop tell him that I considered it at 6:13pm Monday 19th June, 2017.
    meraxes_Halo[counter].Pos[1] = Pos[1]; // I was unsure if it was worth though. 
    meraxes_Halo[counter].Pos[2] = Pos[2];

    counter++;
  } 
  fclose(load_fd);

  fprintf(stderr, "Read in the Pos data.\n");

  // Pos // 

  counter = 0;

  snprintf(buf, MAXLEN, "/lustre/projects/p004_swin/jseiler/meraxes-tiamat/halolist_%s.MetalsColdGas.txt", tag);  
  if(!(load_fd = fopen(buf, "r")))
  {
    printf("can't open file `%s'\n", buf);
    exit(0); 
  }
 
  while (fgets(line, MAXLEN, load_fd))
  {
    if(*line == '#')
	continue;
    sscanf(line, "%lf", &MetalsColdGas); 
    meraxes_Halo[counter].MetalsColdGas = MetalsColdGas;

    counter++;
  } 
  fclose(load_fd);

  fprintf(stderr, "Read in the MetalsColdGas data.\n");
  return num_halos;
}
