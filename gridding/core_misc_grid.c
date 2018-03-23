#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>

#include "core_allvars_grid.h"
#include "core_proto_grid.h"

void count_halomass() 
{
  int i;
  double halomass = 0;
  for(i = 0; i < CUBE(GridSize); ++i)
  {
    halomass += Grid[i].HaloMass;
  }

  printf("Total Halo mass at redshift %.3f is %.3e\n", ZZ[ListOutputGrid[0]], halomass);
}
