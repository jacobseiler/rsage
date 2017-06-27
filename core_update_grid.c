#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>
#include <assert.h>

#ifdef MPI
#include <mpi.h>
#endif

#include "core_allvars_grid.h"
#include "core_proto_grid.h"

// Takes index of a galaxy (p) and maps the properties (SFR, StellarMass etc) onto the grid.

void update_grid_properties(int p, int merged, int GridNr, int filenr)
{

  int i, grid_position;
  double fesc_local;
  float Ngamma_HI, Ngamma_HeI, Ngamma_HeII; // Number of ionizing photons (in s^-1). 
  for(i = 0; i < MAXSNAPS; ++i)
  {
    grid_position = GalGrid[p].History[i];
 
    if (grid_position == -1) continue; // If the galaxy hasn't formed yet, move on to next snapshot. 
    if (i > ListOutputGrid[GridNr]) break; // If the current snapshot redshift is lower than an output redshift, proceed to next one.
    XASSERT(grid_position >= 0 && grid_position < CUBE(GridSize), "The grid index must be between 0 and the GridSize cubed (%d cubed = %d).  The index for Galaxy %d is %d.\n", GridSize, CUBE(GridSize), p, grid_position);

    //These are instantaneous properties measured at a specific redshift. //
    if((i == ListOutputGrid[GridNr]) && (GalGrid[p].StellarMass[i] > 0.0) & (GalGrid[p].SFR[i] > 0.0) & (GalGrid[p].CentralGalaxyMass[i] > 0.0) & (GalGrid[p].LenHistory[i] > HaloPartCut)) // Only want to take those galaxies with a non-zero stellar mass (they actually exist and are evolving).
    {
	
      Grid[grid_position].Sfr += GalGrid[p].SFR[i];
      Grid[grid_position].Count += 1;
      Grid[grid_position].StellarMass += GalGrid[p].StellarMass[i];

      if (PhotonPrescription == 0)
      { 
        update_grid_nion_halo(GridNr); // If we are using a halo-based photon prescription, calculate the number of photons based on halo mass.
      } else if (PhotonPrescription == 1)  // If our photon prescription is using the galaxy photon calculated from the STARBURST spectra.
      {
	calculate_photons(GalGrid[p].SFR[i], GalGrid[p].Z[i], &Ngamma_HI, &Ngamma_HeI, &Ngamma_HeII);
	
        fesc_local = calculate_fesc(p, i, filenr);

        Grid[grid_position].Nion_HI += pow(10, Ngamma_HI)*fesc_local;
	
	XASSERT(Grid[grid_position].Nion_HI >= 0.0 && Grid[grid_position].Nion_HI < 1e100, "Somehow have a grid cell with less than zero HI ionizing photons.  Grid cell = %d, Galaxy Number = %d, Snapshot = %d, Ngamma_HI = %.4e, fesc_local = %.4e\n", grid_position, p, i, Ngamma_HI, fesc_local);
     
      
        Grid[grid_position].Nion_HeI += pow(10, Ngamma_HeI)*fesc_local; 
        Grid[grid_position].Nion_HeII += pow(10, Ngamma_HeII)*fesc_local;
 
      }
    }
      
    // Tracking the Merging Types. //
  }
}

void update_meraxes_grid_properties(int p, int GridNr)
{
  if (meraxes_Halo[p].SFR > 1e-10)
  {
    int grid_position;
    double fesc_local;
    float Ngamma_HI, Ngamma_HeI, Ngamma_HeII; // Number of ionizing photons (in s^-1). 
  
    int x_grid, y_grid, z_grid;
    
    x_grid = meraxes_Halo[p].Pos[0]*GridSize/BoxSize; // Convert the (x,y,z) position to a grid (x,y,z).
    y_grid = meraxes_Halo[p].Pos[1]*GridSize/BoxSize;
    z_grid = meraxes_Halo[p].Pos[2]*GridSize/BoxSize; 
 
    XPRINT(meraxes_Halo[p].Pos[0] < 80, "Halo %d has an x position of %.4f Mpc/h", p, meraxes_Halo[p].Pos[0]); 
    grid_position = (x_grid*GridSize+y_grid)*GridSize+z_grid; // Convert the grid (x,y,z) to a 1D value.
   
    XASSERT(grid_position >= 0 && grid_position < CUBE(GridSize), "The grid index must be between 0 and the GridSize cubed (%d cubed = %d).  The index for Halo %d is %d.\n", GridSize, CUBE(GridSize), p, grid_position);
  
    Grid[grid_position].Sfr += meraxes_Halo[p].SFR; 
    Grid[grid_position].Count += 1;
    Grid[grid_position].StellarMass += meraxes_Halo[p].StellarMass;
  
    XASSERT(PhotonPrescription == 1, "We're using the MERAXES galaxies so we have to be using a galaxy based photon prescription. The current variable (PhotonPrescription) is set to a Halo based model.\n");
  
    double Z = get_metallicity(meraxes_Halo[p].ColdGas, meraxes_Halo[p].MetalsColdGas); 
    
    calculate_photons(meraxes_Halo[p].SFR, Z, &Ngamma_HI, &Ngamma_HeI, &Ngamma_HeII);
    fesc_local = calculate_fesc(p, 999, 999);

    Grid[grid_position].Nion_HI += pow(10, Ngamma_HI)*fesc_local;

    XASSERT(Grid[grid_position].Nion_HI >= 0.0 && Grid[grid_position].Nion_HI < 1e100, "Somehow have a grid cell with less than zero HI ionizing photons.  Grid cell = %d, Galaxy Number = %d, Snapshot = %d, Ngamma_HI = %.4e, fesc_local = %.4e\n", grid_position, p, GridNr, Ngamma_HI, fesc_local);
       
    Grid[grid_position].Nion_HeI += pow(10, Ngamma_HeI)*fesc_local; 
    Grid[grid_position].Nion_HeII += pow(10, Ngamma_HeII)*fesc_local;
  }  
}


// Takes the number of halos loaded in memory (totNHalos) and maps the properties onto the grid. //

void update_grid_halo(int totNHalos, int GridNr)
{
  int i;
  int x_grid, y_grid, z_grid, grid_position; 

  for(i = 0; i < totNHalos; ++i)
  {
    if(Halo[i].SnapNum == ListOutputGrid[GridNr] && Halo[i].Mvir > 0.0) // Checks to make sure the halo existed at the output redshift.
    {

      x_grid = Halo[i].Pos[0]*GridSize/BoxSize; // Convert the (x,y,z) position to a grid (x,y,z).
      y_grid = Halo[i].Pos[1]*GridSize/BoxSize;
      z_grid = Halo[i].Pos[2]*GridSize/BoxSize; 

      grid_position = (x_grid*GridSize+y_grid)*GridSize+z_grid; // Convert the grid (x,y,z) to a 1D value.

      Grid[grid_position].HaloMass += Halo[i].Mvir;
      ++Grid[grid_position].HaloCount; 
    
    }
  }
}

// 

void update_grid_diffuse(int GridNr)
{
 
  char bufz0[MAXLEN];
  FILE* load_fd = NULL;


  snprintf(bufz0, MAXLEN, "%s/%d_Unbound.dat", DiffuseDir, ListOutputGrid[GridNr]);
 
  if(!(load_fd = fopen(bufz0, "r")))
  {
    printf("Can't open file `%s'\n", bufz0);
    exit(0);
  }
   
  int i;
  printf("Reading in diffuse particles for Snapshot %d.\n", ListOutputGrid[GridNr]); 
  for (i = 0; i < CUBE(GridSize); ++i)
    fread(&Grid[i].Diffuse, sizeof(int), 1, load_fd);

//  fread(&GridDiffuse, sizeof(int), CUBE(GridSize), load_fd); 

  fclose(load_fd);
}

void update_grid_density(int GridNr)
{ 
  int i, count = 0;
  double total_density = 0, max_density = 0;
  double tot_StellarMass = 0, tot_DiffuseMass = 0, tot_HaloMass = 0;

  for(i = 0; i < CUBE(GridSize); ++i)
  {
    tot_StellarMass += Grid[i].StellarMass;
    tot_DiffuseMass += Grid[i].Diffuse*PartMass;
    tot_HaloMass += Grid[i].HaloMass;
    total_density += Grid[i].StellarMass + Grid[i].Diffuse*PartMass + Grid[i].HaloMass;  
  }
    
  printf("=======================================\n\n");
  printf("Calculating density for Snapshot %d (z = %.3f).\n", ListOutputGrid[GridNr], ZZ[ListOutputGrid[GridNr]]);
  printf("The total mass in the grid is %.4e (10^10 MSun/h)\n", total_density);
  printf("With %.4e Stellar, %.4e Diffuse and %.4e Halo Mass\n", tot_StellarMass, tot_DiffuseMass, tot_HaloMass);

  total_density /= CUBE(GridSize);

  printf("Average density for a grid cell is %.4f (10^10 MSun/h)\n", total_density);
  printf("One particle has a mass of %.4f (10^10MSun/h)\n", PartMass);
  printf("Hence we make a cell with no particles have a density of %.4f\n", PartMass/total_density);
    
  // Convert to overdensity. //
  for(i = 0; i < CUBE(GridSize); ++i)
  {
    Grid[i].Density = (Grid[i].StellarMass + Grid[i].Diffuse*PartMass + Grid[i].HaloMass)/total_density;
    if (Grid[i].Density > max_density) // Tracks maximum overdensity.
	    max_density = Grid[i].Density;
    if(Grid[i].Density == 0.0) // We don't want a grid cell to have zero density so if it does, just populate it with a single dark matter particle.
    {
      ++count;
      Grid[i].Density = PartMass/total_density;
    }
      
  }

  printf("%d cells with zero Density\n", count);
  printf("The maximum density was %.4f\n", max_density); 
  printf("=======================================\n\n");

}

void update_grid_nion_halo(int GridNr) // Calculates number of ionizing photons using the halos. 
{
  int i;
  double evolve_time;

  printf("Calculating number of photons using halo-based prescription.\n");
  
 
  if (GridNr == NGrid - 1)
	  evolve_time = (Age[ListOutputGrid[GridNr]-1] - Age[ListOutputGrid[GridNr]]) * UnitTime_in_Megayears;
  else
	  evolve_time = (Age[ListOutputGrid[GridNr+1]] - Age[ListOutputGrid[GridNr]]) * UnitTime_in_Megayears;

  for (i = 0; i < CUBE(GridSize); ++i)
  {
    // Using Illiev 2012. 
//    UnitConversion = SOLAR_MASS/0.53/PROTONMASS/SEC_PER_MEGAYEAR;
    Grid[i].Nion_HI = Grid[i].HaloMass*1e10/Hubble_h*SOLAR_MASS*SourceEfficiency*BaryonFrac/0.53/PROTONMASS/evolve_time/SEC_PER_MEGAYEAR;
  }
  
}

void count_grid_properties(int GridNr) // Count number of galaxies/halos in the grid.
{
  int i;

  int GalCount = 0, HaloCount = 0, DiffuseCount = 0, SourcesCount = 0;
  double totPhotons_HI = 0, totPhotons_HeI = 0, totPhotons_HeII = 0, HaloMass = 0;
  for (i = 0; i < CUBE(GridSize); ++i)
  {
    GalCount += Grid[i].Count;
    HaloCount += Grid[i].HaloCount;
    totPhotons_HI += Grid[i].Nion_HI;
    totPhotons_HeI += Grid[i].Nion_HeI;
    totPhotons_HeII += Grid[i].Nion_HeII;
    HaloMass += Grid[i].HaloMass;
    DiffuseCount += Grid[i].Diffuse;
    if(Grid[i].Nion_HI > 0.0)
	SourcesCount++;
  }

  printf("At redshift %.3f (Snapshot %d) there was %d galaxies, %d halos and %d diffuse particles with [%.4e, %.4e, %.4e] {HI, HeI, HeII}  ionizing Photons emitted per second ([%.4e %.4e %.4e] s^-1 Mpc^-3), spread across %d cells (%.4f of the total cells) and %.4e (Msun) halo mass.\n", ZZ[ListOutputGrid[GridNr]], ListOutputGrid[GridNr], GalCount, HaloCount, DiffuseCount, totPhotons_HI, totPhotons_HeI, totPhotons_HeII, totPhotons_HI / pow(BoxSize/Hubble_h,3), totPhotons_HeI / pow(BoxSize/Hubble_h, 3), totPhotons_HeII / pow(BoxSize/Hubble_h,3), SourcesCount, (double)SourcesCount / (CUBE(GridSize)), HaloMass*1e10/Hubble_h);

}

void normalize_photon(int GridNr)
{

  int i;
  double totPhotons_HI_normgrid = 0, totPhotons_HI = 0;
  char buf[MAXLEN], tag[MAXLEN];
  double *NormGrid;
  FILE *load_fd; 

  if ((NGrid - 1 - GridNr) < 10)
    snprintf(tag, MAXLEN, "%s", "00");
  else if ((NGrid - 1 - GridNr) < 100 && (NGrid - 1 - GridNr) >= 10)
    snprintf(tag, MAXLEN, "%s", "0");
  else
    snprintf(tag, MAXLEN, "%s", "");


  NormGrid = malloc(CUBE(GridSize) * sizeof(double));

  snprintf(buf, MAXLEN, "/lustre/projects/p004_swin/jseiler/SAGE_output/1024/grid/January_input/Galaxies_noreion_z5.000_fesc0.15_nionHI_%s%d", tag, (NGrid-1) - GridNr);
  if(!(load_fd = fopen(buf, "r")))
  {
    printf("can't open file `%s'\n", buf);
    exit(0); 
  }

  for (i = 0; i < CUBE(GridSize); ++i)
  {
    fread(&NormGrid[i], 1, sizeof(double), load_fd);
    totPhotons_HI_normgrid += NormGrid[i];
    totPhotons_HI += Grid[i].Nion_HI; 
  }

  double ratio = totPhotons_HI_normgrid/totPhotons_HI;

  for (i = 0; i < CUBE(GridSize); ++i)
  {  
    Grid[i].Nion_HI = Grid[i].Nion_HI * ratio;
  }

  printf("Total photons from normalization grid is %.4e compared to the %.4e photons from fescPrescription %d giving a ratio of %.4e.\n", totPhotons_HI_normgrid, totPhotons_HI, fescPrescription, ratio); 

}

void normalize_slope_photons(int GridNr)
{
  
   int i;
   double totPhotons = 0, targetPhotons = pow(10,-0.8*ZZ[ListOutputGrid[GridNr]] + 64.3); 


   for (i = 0; i < CUBE(GridSize); ++i)
   {
      totPhotons += Grid[i].Nion_HI;
   }

   double ratio = targetPhotons/totPhotons;

   for (i = 0; i < CUBE(GridSize); ++i)
   {
     Grid[i].Nion_HI = Grid[i].Nion_HI * ratio;
   }

   printf("The total number of photons in the grid is %.4e.  We only want %.4e photons for redshift %.4e giving a ratio of %.4e.\n", totPhotons, targetPhotons, ZZ[ListOutputGrid[GridNr]], ratio);
}

// INPUT:
// Star formation rate of the galaxy (in units of Msun/yr) (SFR).
// Metallicity (NOT Solar Units) (Z).
//
// OUTPUT/USE:
// Returns the number of HI ionizing photons for the galaxy.  
//
// NOTE: These relationships have been fit manually from STARBURST99 results.  
// DOUBLE NOTE: These relationships assume a constant starformation scenario; a Starburst scenario is completely different.
void calculate_photons(float SFR, float Z, float *Ngamma_HI, float *Ngamma_HeI, float *Ngamma_HeII)
{

  assert(Ngamma_HI);
  assert(Ngamma_HeI);
  assert(Ngamma_HeII);

  if (SFR == 0)
  {
    *Ngamma_HI = 0;
    *Ngamma_HeI = 0;
    *Ngamma_HeII = 0; 
  }   
  else if (Z < 0.0025) // 11
  { 
    *Ngamma_HI = log10(SFR) + 53.354;
    *Ngamma_HeI = log10(SFR) + 52.727;
    *Ngamma_HeII = log10(SFR) + 48.941;
  }
  else if (Z >= 0.0025 && Z < 0.006) // 12
  {
    *Ngamma_HI = log10(SFR) + 53.290;
    *Ngamma_HeI = log10(SFR) + 52.583;
    *Ngamma_HeII = log10(SFR) + 49.411;
  }
  else if (Z>= 0.006 && Z < 0.014) // 13
  {
    *Ngamma_HI = log10(SFR) + 53.248;
    *Ngamma_HeI = log10(SFR) + 52.481;
    *Ngamma_HeII = log10(SFR) + 49.254;
  }
  else if (Z >= 0.014 && Z < 0.030) // 14
  {
    *Ngamma_HI = log10(SFR) + 53.166;
    *Ngamma_HeI = log10(SFR) + 52.319;
    *Ngamma_HeII = log10(SFR) + 48.596;
  }
  else // 15
  {
    *Ngamma_HI = log10(SFR) + 53.041;
    *Ngamma_HeI = log10(SFR) + 52.052;
    *Ngamma_HeII = log10(SFR) + 47.939;
  }

  if (SFR != 0)
  {
    assert(*Ngamma_HI > 0.0);
    assert(*Ngamma_HeI > 0.0);
    assert(*Ngamma_HeII > 0.0);
  }
}

// Here we calculate the escape fraction for a specific galaxy. 
// The user defines a value (fescPrescription) that determines how to calculate fesc.
// 0: Constant Escape fraction.  User defines the exact value.
// 1: Scaling with Halo Mass using the functional form defined by Kimm et al (2016).
// 2: Power Law as a function of halo mass. The user defines the smallest and largest halo mass in the simulation in addition to the escape fractions at these masses.
// 3: Linear relationship as a function of the ejected fraction of a galaxy.  User defines the escape fraction for ejected fractions of 0 and 1.	
//
// The values of alpha/beta are determined within the 'core_init.c' module.
//
// INPUT: The galaxy index (p). 
// 	: The snapshot index (i)
// 	: The filenumber of the galaxy (filenr); useful for debugging. 
//
// OUTPUT: The escape fraction for the galaxy. 

double calculate_fesc(int p, int i, int filenr)
{

  double fesc_local, halomass, ejectedfraction, SFR;

  if (use_sage == 1)
  {
    halomass = GalGrid[p].CentralGalaxyMass[i];
    ejectedfraction = GalGrid[p].EjectedFraction[i];
    SFR = GalGrid[p].SFR[i];
  }
  else
  { 
    halomass = meraxes_Halo[p].Mvir;
    ejectedfraction = meraxes_Halo[p].EjectedGas / (meraxes_Halo[p].EjectedGas + meraxes_Halo[p].ColdGas + meraxes_Halo[p].HotGas);
    SFR = meraxes_Halo[p].SFR;
  }

  if (fescPrescription == 0) 
    fesc_local = fesc;
  else if (fescPrescription == 1)
    fesc_local = pow(10,1.0 - 0.2*log10(halomass * 1.0e10 / Hubble_h));
  else if (fescPrescription == 2)
    fesc_local = pow(10, log10(alpha) + beta * log10(halomass * 1.0e10 / Hubble_h));
  else if (fescPrescription == 3)	
    fesc_local = beta * ejectedfraction + alpha;
	
  if (fesc_local > 1.0)
  {
    fprintf(stderr, "Had fesc_local = %.4f for galaxy %d in file %d with halo mass %.4e (log Msun), Stellar Mass %.4e (log Msun), SFR %.4e (log Msun yr^-1) and Ejected Fraction %.4e\n", fesc_local, p, filenr, log10(halomass * 1.0e10 / Hubble_h), log10(halomass * 1.0e10 / Hubble_h), log10(SFR), ejectedfraction);
    fesc_local = 1.0;
  
  }
 
  if (fesc_local < 0.0)
  {
    fprintf(stderr, "Had fesc_local = %.4f for galaxy %d in file %d with halo mass %.4e (log Msun), Stellar Mass %.4e (log Msun), SFR %.4e (log Msun yr^-1) and Ejected Fraction %.4e\n", fesc_local, p, filenr, log10(GalGrid[p].CentralGalaxyMass[i] * 1.0e10 / Hubble_h), log10(GalGrid[p].StellarMass[i] * 1.0e10 / Hubble_h), log10(GalGrid[p].SFR[i]), GalGrid[p].EjectedFraction[i]);
    fesc_local = 0.0; 
  }

  return fesc_local;

}

double get_metallicity(double gas, double metals)
{
  double metallicity;

  if(gas > 0.0 && metals > 0.0)
  {
    metallicity = metals / gas;
    if(metallicity < 1.0)
      return metallicity;
    else
      return 1.0;
  }
  else
    return 0.0;

}
