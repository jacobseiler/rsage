#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "main.h"

#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

#define MAXLEN 1024

int32_t read_parameter_file(char *fname, SAGE_params *params) 
{
  FILE *fd;
  char buf[400], buf1[400], buf2[400], buf3[400];
  int i, j, nt = 0;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int errorFlag = 0; 

  char TreeDir_tmp[MAXLEN], TreeName_tmp[MAXLEN], PhotoionDir_tmp[MAXLEN], PhotoionName_tmp[MAXLEN], ReionRedshiftName_tmp[MAXLEN], SnapListFile_tmp[MAXLEN];
  int32_t FirstFile_tmp, LastFile_tmp, LastSnapshotNr_tmp, GridSize_tmp;
  double BoxSize_tmp, Hubble_h_tmp;

#ifdef MPI
  if(ThisTask == 0)
#endif
  printf("\nreading parameter file:\n\n");
  
  strcpy(tag[nt], "TreeName");
  addr[nt] = TreeName_tmp;
  id[nt++] = STRING;

  strcpy(tag[nt], "SimulationDir");
  addr[nt] = TreeDir_tmp;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileWithSnapList");
  addr[nt] = SnapListFile_tmp;
  id[nt++] = STRING;

  strcpy(tag[nt], "PhotoionDir");
  addr[nt] = PhotoionDir_tmp;
  id[nt++] = STRING;

  strcpy(tag[nt], "PhotoionName");
  addr[nt] = PhotoionName_tmp;
  id[nt++] = STRING;

  strcpy(tag[nt], "FirstFile");
  addr[nt] = &FirstFile_tmp;
  id[nt++] = INT;

  strcpy(tag[nt], "LastFile");
  addr[nt] = &LastFile_tmp;
  id[nt++] = INT;

  strcpy(tag[nt], "LastSnapShotNr");
  addr[nt] = &LastSnapshotNr_tmp;
  id[nt++] = INT;

  strcpy(tag[nt], "GridSize");
  addr[nt] = &GridSize_tmp;
  id[nt++] = INT;

  strcpy(tag[nt], "BoxSize");
  addr[nt] = &BoxSize_tmp;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "Hubble_h");
  addr[nt] = &Hubble_h_tmp;
  id[nt++] = DOUBLE;

  strcpy(tag[nt], "ReionRedshiftName");
  addr[nt] = ReionRedshiftName_tmp;
  id[nt++] = STRING;

  if((fd = fopen(fname, "r")))
  {
    while(!feof(fd))
    {
      *buf = 0;
      fgets(buf, 200, fd);
      if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
        continue;

      if(buf1[0] == '%' || buf1[0] == '-')
        continue;

      for(i = 0, j = -1; i < nt; i++)
        if(strcmp(buf1, tag[i]) == 0)
      {
        j = i;
        tag[i][0] = 0;
        break;
      }

      if(j >= 0)
      {
#ifdef MPI
        if(ThisTask == 0)
#endif  
        switch (id[j])
        {
          case DOUBLE:
          *((double *) addr[j]) = atof(buf2);
          break;
          case STRING:
          strcpy(addr[j], buf2);
          break;
          case INT:
          *((int *) addr[j]) = atoi(buf2);
          break;
        }
      }
      else
      {
        //printf("Error in file %s:   Tag '%s' not allowed or multiple defined.\n", fname, buf1);        
      }
    }
    fclose(fd);          
  }
  else
  {
    printf("Parameter file %s not found.\n", fname);
    errorFlag = 1;
  }

  for(i = 0; i < nt; i++)
  {
    if(*tag[i])
    {
      printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
      errorFlag = 1;
    }
  }
	
  if (errorFlag == 1)
  {
    return EXIT_FAILURE;
  }
	printf("\n");

  *params = malloc(sizeof(struct SAGE_PARAMETERS));

  (*params)->TreeDir = TreeDir_tmp;
  (*params)->TreeName = TreeName_tmp;
  (*params)->PhotoionDir = PhotoionDir_tmp;
  (*params)->SnapListFile = SnapListFile_tmp;
  (*params)->PhotoionName = PhotoionName_tmp;
  (*params)->ReionRedshiftName = ReionRedshiftName_tmp;

  (*params)->FirstFile = FirstFile_tmp;
  (*params)->LastFile = LastFile_tmp;
  (*params)->LastSnapshotNr = LastSnapshotNr_tmp;
  (*params)->TreeDir = TreeDir_tmp;

  (*params)->GridSize = GridSize_tmp;
  (*params)->BoxSize = BoxSize_tmp;
  (*params)->Hubble_h = Hubble_h_tmp;

  return EXIT_SUCCESS;

}
