#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

#define MAXLEN 1024

int32_t read_parameter_file(char *fname, char **treedir, char **treename, char **photoiondir, char **photoionname, char **reionredshiftname, int32_t *numfiles)
{
  FILE *fd;
  char buf[400], buf1[400], buf2[400], buf3[400];
  int i, j, nt = 0;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int errorFlag = 0; 

  char treedir_tmp[MAXLEN], treename_tmp[MAXLEN], photoiondir_tmp[MAXLEN], photoionname_tmp[MAXLEN], reionredshiftname_tmp[MAXLEN];
  int32_t FirstFile, LastFile;

#ifdef MPI
  if(ThisTask == 0)
#endif
  printf("\nreading parameter file:\n\n");
  
  strcpy(tag[nt], "TreeName");
  addr[nt] = treename_tmp;
  id[nt++] = STRING;

  strcpy(tag[nt], "SimulationDir");
  addr[nt] = treedir_tmp;
  id[nt++] = STRING;

  strcpy(tag[nt], "PhotoionDir");
  addr[nt] = photoiondir_tmp;
  id[nt++] = STRING;

  strcpy(tag[nt], "PhotoionName");
  addr[nt] = photoionname_tmp;
  id[nt++] = STRING;

  strcpy(tag[nt], "FirstFile");
  addr[nt] = &FirstFile;
  id[nt++] = INT;

  strcpy(tag[nt], "LastFile");
  addr[nt] = &LastFile;
  id[nt++] = INT;

  strcpy(tag[nt], "ReionRedshiftName");
  addr[nt] = reionredshiftname_tmp;
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

  *treedir = treedir_tmp;
  *treename = treename_tmp;
  *photoiondir = photoiondir_tmp; 
  *photoionname = photoionname_tmp; 
  *reionredshiftname = reionredshiftname_tmp; 
  *numfiles = LastFile - FirstFile + 1;

  return EXIT_SUCCESS;

}
