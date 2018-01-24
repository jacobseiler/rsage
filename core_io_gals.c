#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
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

int32_t load_gals(char *fname)
{

  int i;
  FILE *infile;

  infile = fopen(fname, "rb");
  if (infile == NULL) 
  {
    printf("can't open file `%s'\n", fname);
    exit(0);
  }

#ifdef DEBUG_GALS 
  printf("Loading Galaxies from file %s\n", fname);
#endif
  fread(&Ntrees, 1, sizeof(int), infile);
  fread(&NtotGals, sizeof(int), 1, infile);

  GalsForTree = malloc(Ntrees * sizeof(int));
  fread(GalsForTree, Ntrees, sizeof(int), infile);

  GalGrid = mymalloc(NtotGals * sizeof(struct GALAXY_GRID));
  
  for (i = 0; i < NtotGals; ++i)
  {

    fread(&GalGrid[i].HaloNr, sizeof(int), 1, infile); 
 
    GalGrid[i].History = malloc(sizeof(*(GalGrid[i].History)) * MAXSNAPS);
    if (GalGrid[i].History == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for GalGrid.History for galaxy %d\n", i);
      exit(EXIT_FAILURE);
    }
    fread(GalGrid[i].History, sizeof(*(GalGrid[i].History)), MAXSNAPS, infile);

    GalGrid[i].StellarMass = malloc(sizeof(*(GalGrid[i].StellarMass)) * MAXSNAPS);
    if (GalGrid[i].StellarMass == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for GalGrid.StellarMass for galaxy %d\n", i);
      exit(EXIT_FAILURE);
    }
    fread(GalGrid[i].StellarMass, sizeof(*(GalGrid[i].StellarMass)), MAXSNAPS, infile);
    
    GalGrid[i].SFR = malloc(sizeof(*(GalGrid[i].SFR)) * MAXSNAPS); 
    if (GalGrid[i].SFR == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for GalGrid.SFR for galaxy %d\n", i);
      exit(EXIT_FAILURE);
    }
    fread(GalGrid[i].SFR, sizeof(*(GalGrid[i].SFR)), MAXSNAPS, infile);

    GalGrid[i].Z = malloc(sizeof(*(GalGrid[i].Z)) * MAXSNAPS);
    if (GalGrid[i].Z == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for GalGrid.Z for galaxy %d\n", i);
      exit(EXIT_FAILURE);
    }
    fread(GalGrid[i].Z, sizeof(*(GalGrid[i].Z)), MAXSNAPS, infile);

    GalGrid[i].CentralGalaxyMass = malloc(sizeof(*(GalGrid[i].CentralGalaxyMass)) * MAXSNAPS);
    if (GalGrid[i].CentralGalaxyMass == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for GalGrid.CentralGalaxyMass for galaxy %d\n", i);
      exit(EXIT_FAILURE);
    }
    fread(GalGrid[i].CentralGalaxyMass, sizeof(*(GalGrid[i].CentralGalaxyMass)), MAXSNAPS, infile);

    GalGrid[i].MfiltGnedin = malloc(sizeof(*(GalGrid[i].MfiltGnedin)) * MAXSNAPS);
    if (GalGrid[i].MfiltGnedin == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for GalGrid.MfiltGnedin for galaxy %d\n", i);
      exit(EXIT_FAILURE);
    }
    fread(GalGrid[i].MfiltGnedin, sizeof(*(GalGrid[i].MfiltGnedin)), MAXSNAPS, infile);

    GalGrid[i].MfiltSobacchi = malloc(sizeof(*(GalGrid[i].MfiltSobacchi)) * MAXSNAPS);
    if (GalGrid[i].MfiltSobacchi == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for GalGrid.MfiltSobacchi for galaxy %d\n", i);
      exit(EXIT_FAILURE);
    }
    fread(GalGrid[i].MfiltSobacchi, sizeof(*(GalGrid[i].MfiltSobacchi)), MAXSNAPS, infile);
 
    GalGrid[i].EjectedFraction = malloc(sizeof(*(GalGrid[i].EjectedFraction)) * MAXSNAPS);
    if (GalGrid[i].EjectedFraction == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for GalGrid.EjectedFraction for galaxy %d\n", i);
      exit(EXIT_FAILURE);
    }
    fread(GalGrid[i].EjectedFraction, sizeof(*(GalGrid[i].EjectedFraction)), MAXSNAPS, infile);

    GalGrid[i].LenHistory = malloc(sizeof(*(GalGrid[i].LenHistory)) * MAXSNAPS);
    if (GalGrid[i].LenHistory == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for GalGrid.LenHistory for galaxy %d\n", i);
      exit(EXIT_FAILURE);
    }
    fread(GalGrid[i].LenHistory, sizeof(*(GalGrid[i].LenHistory)), MAXSNAPS, infile);

    GalGrid[i].OutflowRate= malloc(sizeof(*(GalGrid[i].OutflowRate)) * MAXSNAPS);
    if (GalGrid[i].OutflowRate == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for GalGrid.OutflowRate for galaxy %d\n", i);
      exit(EXIT_FAILURE);
    }
    fread(GalGrid[i].OutflowRate, sizeof(*(GalGrid[i].OutflowRate)), MAXSNAPS, infile);

    GalGrid[i].InfallRate= malloc(sizeof(*(GalGrid[i].InfallRate)) * MAXSNAPS);
    if (GalGrid[i].InfallRate == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for GalGrid.InfallRate for galaxy %d\n", i);
      exit(EXIT_FAILURE);
    }
    fread(GalGrid[i].InfallRate, sizeof(*(GalGrid[i].InfallRate)), MAXSNAPS, infile);

    GalGrid[i].EjectedMass= malloc(sizeof(*(GalGrid[i].EjectedMass)) * MAXSNAPS);
    if (GalGrid[i].EjectedMass == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for GalGrid.EjectedMass for galaxy %d\n", i);
      exit(EXIT_FAILURE);
    }
    fread(GalGrid[i].EjectedMass, sizeof(*(GalGrid[i].EjectedMass)), MAXSNAPS, infile);

    GalGrid[i].QuasarActivity= malloc(sizeof(*(GalGrid[i].QuasarActivity)) * MAXSNAPS);
    if (GalGrid[i].QuasarActivity == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for GalGrid.QuasarActivity for galaxy %d\n", i);
      exit(EXIT_FAILURE);
    }
    fread(GalGrid[i].QuasarActivity, sizeof(*(GalGrid[i].QuasarActivity)), MAXSNAPS, infile);

    GalGrid[i].DynamicalTime= malloc(sizeof(*(GalGrid[i].DynamicalTime)) * MAXSNAPS);
    if (GalGrid[i].DynamicalTime == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for GalGrid.DynamicalTime for galaxy %d\n", i);
      exit(EXIT_FAILURE);
    }
    fread(GalGrid[i].DynamicalTime, sizeof(*(GalGrid[i].DynamicalTime)), MAXSNAPS, infile);

    GalGrid[i].QuasarSubstep= malloc(sizeof(*(GalGrid[i].QuasarSubstep)) * MAXSNAPS);
    if (GalGrid[i].QuasarSubstep == NULL)
    {
      fprintf(stderr, "Cannot allocate memory for GalGrid.QuasarSubstep for galaxy %d\n", i);
      exit(EXIT_FAILURE);
    }
    fread(GalGrid[i].QuasarSubstep, sizeof(*(GalGrid[i].QuasarSubstep)), MAXSNAPS, infile);

#ifdef DEBUG_GALS
    if (i == 0)
    {
      int tmp = MAXSNAPS - 1;
      printf("History[Final] = %d\tStellarMass = %.4f\tSFR = %.4f\tZ = %.4f\tCentralGalaxyMass = %.4f\tMfiltGnedin = %.4f\tMfiltSobacchi = %.4f\tEjectedFraction = %.4f\tLenHistory = %d\tOutflowRate = %.4f\tInfallRate = %.4f\tEjectedMass = %.4f\tQuasarActivity = %d\tDynamicalTime = %.4f\tQuasarSubstep = %d\n", GalGrid[i].History[tmp], GalGrid[i].StellarMass[tmp], GalGrid[i].SFR[tmp], GalGrid[i].Z[tmp], GalGrid[i].CentralGalaxyMass[tmp], GalGrid[i].MfiltGnedin[tmp], GalGrid[i].MfiltSobacchi[tmp], GalGrid[i].EjectedFraction[tmp], GalGrid[i].LenHistory[tmp], GalGrid[i].OutflowRate[tmp], GalGrid[i].InfallRate[tmp], GalGrid[i].EjectedMass[tmp], GalGrid[i].QuasarActivity[tmp], GalGrid[i].DynamicalTime[tmp], GalGrid[i].QuasarSubstep[tmp]);
      fclose(infile);
      return EXIT_FAILURE; 
    }
#endif

  } // Galaxy Loop

  printf("Read in %ld total galaxies.\n", (long)NtotGals);

  fclose(infile);
  free(GalsForTree);

  return EXIT_SUCCESS;

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
      free(GalGrid[i].OutflowRate); 
      free(GalGrid[i].InfallRate); 
      free(GalGrid[i].EjectedMass); 
      free(GalGrid[i].QuasarActivity); 
      free(GalGrid[i].DynamicalTime);
      free(GalGrid[i].QuasarSubstep); 
    }

    myfree(GalGrid);
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

// NOTE: BEYOND DEPRECATED

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
    XPRINT(Pos[0] < 80, "Galaxy has a position of %.4f\n", Pos[0]);
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

