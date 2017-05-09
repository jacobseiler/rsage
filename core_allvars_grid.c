#include "core_allvars_grid.h"

// galaxy data 
struct GRID
  *Grid;

struct GALAXY_INPUT
  *Gal;

struct halo_data *Halo;

double UnitLength_in_cm = 3.08578e+24,
  UnitTime_in_s,
  UnitVelocity_in_cm_per_s = 100000,
  UnitMass_in_g = 1.989e+43,
  RhoCrit,
  UnitPressure_in_cgs,
  UnitDensity_in_cgs,
  UnitCoolingRate_in_cgs,
  UnitEnergy_in_cgs,
  UnitTime_in_Megayears, 
  G,
  Hubble,
  a0, ar;

gsl_rng *random_generator;

int Ntrees;
int NtotGals; // Number of galaxies per file.
int *GalsForTree;

double ZZ[ABSOLUTEMAXSNAPS];
double AA[ABSOLUTEMAXSNAPS];
double Age[ABSOLUTEMAXSNAPS];

// Parameter file variables //

int FirstFile;    // first and last file for processing 
int LastFile;

int NumGals;     // Total number of galaxies stored for current tree 
int MaxGals;     // Maximum number of galaxies allowed for current tree  

int LastSnapShotNr;

char GalaxiesInputDir[512];
char GridOutputDir[512];
char FileNameGalaxies[512];
char FileNameMergedGalaxies[512];
char TreeName[512];
char SimulationDir[512];
char DiffuseDir[512];
char FileWithSnapList[512];

int TotHalos;
int totNHalos;
int TotGalaxies[ABSOLUTEMAXSNAPS];
int *TreeNgals[ABSOLUTEMAXSNAPS];

int *FirstHaloInSnap;

int *TreeNHalos;
int *TreeFirstHalo;

#ifdef MPI
int ThisTask, NTask, nodeNameLen;
char *ThisNode;
#endif

double Omega;
double OmegaLambda;
double PartMass;
double BoxSize;
int GridSize;
int GridSnap;
int NGrid;
int LastOutputSnap;

double Hubble_h;
double EnergySNcode, EnergySN;
double EtaSNcode, EtaSN;

// recipe flags 
int ReionizationOn;
int SupernovaRecipeOn;
int DiskInstabilityOn;
int AGNrecipeOn;
int SFprescription;

// recipe parameters 
double RecycleFraction;
double Yield;
double FracZleaveDisk;
double ReIncorporationFactor;
double ThreshMajorMerger;
double BaryonFrac;
double SfrEfficiency;
double FeedbackReheatingEpsilon;
double FeedbackEjectionEfficiency;
double RadioModeEfficiency;
double QuasarModeEfficiency;
double BlackHoleGrowthRate;
double Reionization_z0;
double Reionization_zr;
double ThresholdSatDisruption;

int MAXSNAPS;
int NOUT;
int ListOutputSnaps[ABSOLUTEMAXSNAPS];
int ListOutputGrid[ABSOLUTEMAXSNAPS];
int Snaplistlen;

// Variables for tracking Merger Types //
int *GalaxyCount;
int *MajorMergerCount;
int *MinorMergerCount;
int *DiskCount;
int *ICSCount;

double SourceEfficiency;
int PhotonPrescription;
int Verbose;
int fescPrescription;
double fesc;

double tot_Halomass;

double MH_min;
double MH_max;
double alpha;
double beta;
double kappa;
double delta;
int HaloPartCut;
