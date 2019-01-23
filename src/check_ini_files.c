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
#include "sage/self_consistent/selfcon_grid.h"

#ifdef RSAGE

#include "../grid-model/src/confObj.h"

#endif

int32_t check_directory(char *path)
{
  struct stat st = {0};

  if (stat(path, &st) == - 1)
  {
    mkdir(path, 0700);
    printf("Made directory %s.\n", path);
  }

  return EXIT_SUCCESS;
}

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

int32_t check_sage_ini(char *RunPrefix, char *OutputDir, char *GalaxyOutputDir,
                       char *GridOutputDir, char *PhotoionDir, char *PhotoionName,
                       char *ReionRedshiftName, char **argv, int32_t ThisTask)
{

  char buf[MAX_STRING_LEN];

  // Check the Galaxy output directory.
  if (strncmp(GalaxyOutputDir, "None", 4) == 0)
  {
    snprintf(GalaxyOutputDir, MAX_STRING_LEN - 1, "%s/galaxies", OutputDir);
    printf("The `GalaxyOutputDir` variable has been updated to %s\n", GalaxyOutputDir);
  }

  // Grid output directory.
  if (strncmp(GridOutputDir, "None", 4) == 0)
  {
    snprintf(GridOutputDir, MAX_STRING_LEN - 1, "%s/grids", OutputDir);
    printf("The `GridOutputDir` variable has been updated to %s\n", GridOutputDir);
  }

  // Photoionization rate directory. 
  if (strncmp(PhotoionDir, "None", 4) == 0)
  {
    snprintf(PhotoionDir, MAX_STRING_LEN - 1, "%s/grids/cifog", OutputDir);
    printf("The `PhotoionDir` variable has been updated to %s\n", PhotoionDir);
  }

  // Photoionization rate file names. 
  if (strncmp(PhotoionName, "None", 4) == 0)
  {
    snprintf(PhotoionName, MAX_STRING_LEN - 1, "%s_photHI", RunPrefix);
    printf("The `PhotoionName` variable has been updated to %s\n", PhotoionName);
  }

  // Reionization rdshift rate file names. 
  if (strncmp(ReionRedshiftName, "None", 4) == 0)
  {
    snprintf(ReionRedshiftName, MAX_STRING_LEN - 1, "%s_reionization_redshift", RunPrefix);
    printf("The `PhotoionName` variable has been updated to %s\n", PhotoionName);
  }

  // Now check that the sub-directory structure exists.
  if (ThisTask == 0)
  {
    // Outer level directory.
    check_directory(OutputDir);

    // Directory for the galaxies.
    check_directory(GalaxyOutputDir);

    // Directories for all the grids.
    check_directory(GridOutputDir);

    // Directory for the ionizing photon grids.
    snprintf(buf, MAX_STRING_LEN - 1, "%s/nion", GridOutputDir);
    check_directory(buf);

    // Directory for the reionization fields.
    snprintf(buf, MAX_STRING_LEN - 1, "%s/cifog", GridOutputDir);
    check_directory(buf);

    // Directory for the reionization baryonic modifier files.
    snprintf(buf, MAX_STRING_LEN - 1, "%s/cifog/reionization_modifiers", GridOutputDir);
    check_directory(buf);

    // Directory for the ini files.
    snprintf(buf, MAX_STRING_LEN - 1, "%s/ini_files", OutputDir);
    check_directory(buf);

    // Directory for the slurm files.
    snprintf(buf, MAX_STRING_LEN - 1, "%s/slurm_files", OutputDir);
    check_directory(buf);

    // Directory for the log files.
    snprintf(buf, MAX_STRING_LEN - 1, "%s/log_files", OutputDir);
    check_directory(buf);

    // For log keeping purposes, copy the parameter files to the ini_files
    // directory. Furthermore, let's add some extra information to the top
    // of these ini files denoting when the code was executed and the git
    // version.

    // Get the current time. Snippet largely adapted from:
    // https://stackoverflow.com/questions/1442116/how-to-get-the-date-and-time-values-in-a-c-program
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    char time_string[MAX_STRING_LEN];
    snprintf(time_string, MAX_STRING_LEN - 1, "Code was executed on: %d-%02d-%02d %02d:%02d:%02d",
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
    system(buf);
  }

  return EXIT_SUCCESS;
}

#ifdef RSAGE

int32_t check_cifog_ini(confObj_t *simParam, char *OutputDir, char *RunPrefix, char **argv, int32_t ThisTask)
{

  char buf[MAX_STRING_LEN], nion_prefix[MAX_STRING_LEN];

  // We have ensured all the directory structure is correct in `check_sage_ini`.

  if (ThisTask == 0)
  {
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);
    char time_string[MAX_STRING_LEN];
    snprintf(time_string, MAX_STRING_LEN - 1, "Code was executed on: %d-%02d-%02d %02d:%02d:%02d",
             tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);

    // First copy over the cifog ini file and append the time we ran the code
    // and the Git version.
    char *cifog_ini_file = strrchr(argv[2], '/') + 1;

    snprintf(buf, MAX_STRING_LEN - 1, "cp %s %s/ini_files/%s", argv[2], OutputDir, cifog_ini_file); 
    system(buf);

    snprintf(buf, MAX_STRING_LEN - 1, "sed -i '1s/^/%% %s\\n%% Git Version: %s\\n/' %s/ini_files/%s", time_string, VERSION, OutputDir, cifog_ini_file);
    system(buf);
  }

  // The prefix for the ionizing photon file depends upon the fesc prescription
  // chose in addition to the constants used.
  get_nion_prefix(nion_prefix);

  if (strncmp((*simParam)->nion_file, "None", 4) == 0)
  {
    snprintf((*simParam)->nion_file, MAX_STRING_LEN - 1, "%s/grids/nion/%s_%s_nionHI", OutputDir, RunPrefix, nion_prefix);
    printf("The `inputNionFile` variable has been updated to %s\n", (*simParam)->nion_file);
  }

  // Reionization maps files.
  if (strncmp((*simParam)->out_XHII_file, "None", 4) == 0)
  {
    snprintf((*simParam)->out_XHII_file, MAX_STRING_LEN - 1, "%s/grids/nion/%s_XHII", OutputDir, RunPrefix);
    printf("The `output_XHII_file` variable has been updated to %s\n", (*simParam)->out_XHII_file);
  }

  // Photoionization rate files.
  if (strncmp((*simParam)->out_photHI_file, "None", 4) == 0)
  {
    snprintf((*simParam)->out_photHI_file, MAX_STRING_LEN - 1, "%s/grids/nion/%s_photHI", OutputDir, RunPrefix);
    printf("The `output_photHI_file` variable has been updated to %s\n", (*simParam)->out_photHI_file);
  }

  // Restart files.
  if (strncmp((*simParam)->out_restart_file, "None", 4) == 0)
  {
    snprintf((*simParam)->out_restart_file, MAX_STRING_LEN - 1, "%s/grids/nion/%s_restart", OutputDir, RunPrefix);
    printf("The `output_restart_file` variable has been updated to %s\n", (*simParam)->out_restart_file);
  }

  return EXIT_SUCCESS;
}

#endif
