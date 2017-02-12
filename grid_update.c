#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"

void grid_update_mass_metals(double mass, double metallicity, int GridPos, int SnapNum)
{

  int i;

  for (i = SnapNum; i < MAXSNAPS; ++i)
  {
    Grid[GridPos].EjectedMass[i] += mass;
    Grid[GridPos].MetalsEjectedMass[i] += mass*metallicity;
  }


}

void grid_update_mass_metals_mass(double mass, double metals_mass, int GridPos, int SnapNum)
{

  int i;

  for (i = SnapNum; i < MAXSNAPS; ++i)
  {
    Grid[GridPos].EjectedMass[i] += mass;
    Grid[GridPos].MetalsEjectedMass[i] += metals_mass;
  }

}

