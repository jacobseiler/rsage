#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <assert.h>
#include <time.h>

#include "main.h"

#include "../grid-model/src/confObj.h"

/*
Within the SAGE and cifog .ini files, some directory paths can be set to the
value of `None`.  In these instances, we want to use the runtime parameters to
determine the correct paths.

This function checks these fields and updates them if necessary.

Parameters
----------


Returns
----------


Pointer Updates
----------

Units
----------

Notes
----------

strncmp returns 0 if the two strings are equal.
*/

int32_t check_ini_files(char *OutputDir, char *GalaxyOutputDir, char *GridOutputDir,
                        char **argv, int32_t ThisTask)
{

  struct stat st = {0};
  char buf[MAX_STRING_LEN];

  // Check the Galaxy output directory.

  if (strncmp(GalaxyOutputDir, "None", 4) == 0)
  {
    snprintf(GalaxyOutputDir, MAX_STRING_LEN - 1, "%s/galaxies", OutputDir);
    printf("The GalaxyOutputDir has been updated to %s\n", GalaxyOutputDir);
  }

  // Check the Grid output directory.
  if (strncmp(GridOutputDir, "None", 4) == 0)
  {
    snprintf(GridOutputDir, MAX_STRING_LEN - 1, "%s/grids", OutputDir);
    printf("The GridOutputDir has been updated to %s\n", GalaxyOutputDir);
  }

  // Now check that the sub-directory structure exists.
  // Outer directory.
  if (ThisTask == 0)
  {
    if (stat(OutputDir, &st) == - 1)
    {
      mkdir(OutputDir, 0700);
      printf("Made directory %s.\n", OutputDir);
    }

    // Galaxy directory.
    if (stat(GalaxyOutputDir, &st) == - 1)
    {
      mkdir(GalaxyOutputDir, 0700);
      fprintf(stderr, "Made directory %s.\n", GalaxyOutputDir);
    }

    // Grid Directory.
    if (stat(GridOutputDir, &st) == - 1)
    {
      mkdir(GridOutputDir, 0700);
      fprintf(stderr, "Made directory %s.\n", GridOutputDir);
    }

    // Directory for the Nion files.
    snprintf(buf, MAX_STRING_LEN - 1, "%s/nion", GridOutputDir);
    if (stat(buf, &st) == - 1)
    {
      mkdir(buf, 0700);
      fprintf(stderr, "Made directory %s.\n", buf);
    }

    // Directory for the cifog files.
    snprintf(buf, MAX_STRING_LEN - 1, "%s/cifog", GridOutputDir);
    if (stat(buf, &st) == - 1)
    {
      mkdir(buf, 0700);
      fprintf(stderr, "Made directory %s.\n", buf);
    }

    // Directory for the reionization files.
    snprintf(buf, MAX_STRING_LEN - 1, "%s/cifog/reionization_modifiers", GridOutputDir);
    if (stat(buf, &st) == - 1)
    {
      mkdir(buf, 0700);
      fprintf(stderr, "Made directory %s.\n", buf);
    }

    // Directory for the ini files.
    snprintf(buf, MAX_STRING_LEN - 1, "%s/ini_files", OutputDir);
    if (stat(buf, &st) == - 1)
    {
      mkdir(buf, 0700);
      fprintf(stderr, "Made directory %s.\n", buf);
    }

    // Directory for the slurm files.
    snprintf(buf, MAX_STRING_LEN - 1, "%s/slurm_files", OutputDir);
    if (stat(buf, &st) == - 1)
    {
      mkdir(buf, 0700);
      fprintf(stderr, "Made directory %s.\n", buf);
    }

    // Directory for the log files.
    snprintf(buf, MAX_STRING_LEN - 1, "%s/log_files", OutputDir);
    if (stat(buf, &st) == - 1)
    {
      mkdir(buf, 0700);
      fprintf(stderr, "Made directory %s.\n", buf);
    }

    // For log keeping purposes, copy the parameter files to the ini_files
    // directory. Furthermore, let's add some extra information to the top
    // of these ini files denoting when the code was executed and the git
    // version.

    // Get the current time. Snippet largely adapted from:
    // https://stackoverflow.com/questions/1442116/how-to-get-the-date-and-time-values-in-a-c-program
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    char time_string[MAX_STRING_LEN];
    snprintf(time_string, MAX_STRING_LEN - 1, "Code was executed on: %d-%d-%d %d:%d:%d",
             tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);


    // RSAGE is called with the SAGE ini file as the first input. This will be
    // the full path to the ini file, so grab just the file name. Taken from:
    // https://bytes.com/topic/c/answers/214979-how-get-filename-using-c-language
    char *SAGE_ini_file = strrchr(argv[1], '/') + 1;
 
    // Copy the SAGE ini file.
    snprintf(buf, MAX_STRING_LEN - 1, "cp %s %s/ini_files/%s", argv[1], OutputDir, SAGE_ini_file); 
    system(buf);

    // At the top of the file, append the time of run and the Git version.
    // https://stackoverflow.com/questions/20543289/how-do-i-append-a-line-to-the-beginning-of-a-very-large-file-in-linux
    snprintf(buf, MAX_STRING_LEN - 1, "sed -i '1s/^/%% %s\\n%% Git Version: %s\\n/' %s/ini_files/%s", time_string, VERSION, OutputDir, SAGE_ini_file);
    printf("%s\n", buf);
    system(buf);

    // The cifog will only be passed if we're using the fully coupled RSAGE.
#ifdef RSAGE
    char *cifog_ini_file = strrchr(argv[2], '/') + 1;
 
    snprintf(buf, MAX_STRING_LEN - 1, "cp %s %s/ini_files/%s", argv[2], OutputDir, cifog_ini_file); 
    system(buf);

    snprintf(buf, MAX_STRING_LEN - 1, "sed -i '1s/^/%% %s\\n%% Git Version: %s\\n/' %s/ini_files/%s", time_string, VERSION, OutputDir, cifog_ini_file);
    printf("%s\n", buf);
    system(buf);

#endif

  }


  return EXIT_FAILURE;

  // Next adjust the nion_file variable using the parameters from the SAGE ini file.
  

  //return EXIT_SUCCESS;
}
