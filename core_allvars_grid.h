#ifndef ALLVARSGRID_H
#define ALLVARSGRID_H

#include <stdio.h>
#include <gsl/gsl_rng.h>

#define ABORT(sigterm)                                                  \
do {                                                                \
  printf("Error in file: %s\tfunc: %s\tline: %i\n", __FILE__, __FUNCTION__, __LINE__); \
  myexit(sigterm);                                                \
} while(0)

#define  STEPS 10         // Number of integration intervals between two snapshots
#define  MAXGRIDOUT 30
#define  MAXGALFAC 1
#define  ALLOCPARAMETER 10.0
#define  MAX_NODE_NAME_LEN 50
#define  ABSOLUTEMAXSNAPS 1000


#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  HUBBLE      3.2407789e-18   /* in h/sec */

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#define	 CUBE(x) (x*x*x)

#define MAXLEN 1024

#ifdef NDEBUG
#define XASSERT(EXP, ...)                                do{} while(0)
#else
#define XASSERT(EXP, ...)                                              \
    do { if (!(EXP)) {                                                  \
            printf("Error in file: %s\tfunc: %s\tline: %d with expression `"#EXP"'\n", __FILE__, __FUNCTION__, __LINE__); \
            printf(__VA_ARGS__);                                        \
            fflush(stdout);                                             \
            exit(EXIT_FAILURE);                                         \
        } \
    } while (0)
#endif

struct GALAXY_INPUT
{
  int   SnapNum;
  int   Type;

  long long   GalaxyIndex;
  long long   CentralGalaxyIndex;
  int   SAGEHaloIndex;
  int   SAGETreeIndex;
  int   SimulationFOFHaloIndex;
  
  int   mergeType;  //0=none; 1=minor merger; 2=major merger; 3=disk instability; 4=disrupt to ICS
  int   mergeIntoID;
  int   mergeIntoSnapNum;
  float dT;

  // (sub)halo properties
  float Pos[3];
  float Vel[3];
  float Spin[3];
  int   Len;   
  float Mvir;
  float CentralMvir;
  float Rvir;
  float Vvir;
  float Vmax;
  float VelDisp;

  // baryonic reservoirs 
  float ColdGas;
  float StellarMass;
  float BulgeMass;
  float HotGas;
  float EjectedMass;
  float BlackHoleMass;
  float ICS;

  // metals
  float MetalsColdGas;
  float MetalsStellarMass;
  float MetalsBulgeMass;
  float MetalsHotGas;
  float MetalsEjectedMass;
  float MetalsICS;

  // to calculate magnitudes
  float SfrDisk;
  float SfrBulge;
  float SfrDiskZ;
  float SfrBulgeZ;

  // misc 
  float DiskScaleRadius;
  float Cooling;
  float Heating;
  float QuasarModeBHaccretionMass;
  float TimeOfLastMajorMerger;
  float TimeOfLastMinorMerger;
  float OutflowRate;

  // infall properties
  float infallMvir;
  float infallVvir;
  float infallVmax;

  // NOTE: Because the code has moved to pointers for the galaxy-grid information
  // They are no longer stored explicitly in the galaxy input struct.
 
}*Gal;

struct GALAXY_GRID
{

  int *History; // Integers that describe the grid position at every redshift.
  float *StellarMass; // Units of 1.0e10 Msun/h.
  float *SFR; // Units of Msun yr^-1
  float *Z; // NOT solar metallicity, actual metallicity.
  float *CentralGalaxyMass; // Units of 1.0e10 Msun/h.
  float *Pad; // Used just to pad things out.
  double *MfiltGnedin; // Units of 1.0e10 Msun/h.
  double *MfiltSobacchi; // Units of 1.0e10 Msun/h.
  float *EjectedFraction; // Unitless (fraction).
  int *LenHistory; // Number of particles in FoF Halo.

}*GalGrid;

struct GRID
{

  double Sfr;
  double StellarMass;

  double Density; // Overdensity, rho/<rho>.
  double Nion_HI;
  double Nion_HeI;
  double Nion_HeII;

  int Diffuse;
  int Count;

//  double *ActiveTime;

  double HaloMass; // Halo Mvir mass.
  int HaloCount;
  
}
*Grid;

struct halo_data
{
  // merger tree pointers 
  int Descendant;
  int FirstProgenitor;
  int NextProgenitor;
  int FirstHaloInFOFgroup;
  int NextHaloInFOFgroup;

  // properties of halo 
  int Len;
  float M_Mean200, Mvir, M_TopHat;  // for Millennium, Mvir=M_Crit200
  float Pos[3];
  float Vel[3];
  float VelDisp;
  float Vmax;
  float Spin[3];
  long long MostBoundID;

  // original position in simulation tree files
  int SnapNum;
  int FileNr;
  int SubhaloIndex;
  float SubHalfMass;
}
*Halo;

struct meraxes_halo_data
{
  double SFR;
  double StellarMass;
  double Mvir;
  double ColdGas;
  double HotGas;
  double EjectedGas;
  double Pos[3]; 

}
*meraxes_Halo;

extern double UnitLength_in_cm,
  UnitTime_in_s,
  UnitVelocity_in_cm_per_s,
  UnitMass_in_g,
  RhoCrit,
  UnitPressure_in_cgs,
  UnitDensity_in_cgs,
  UnitCoolingRate_in_cgs,
  UnitEnergy_in_cgs,
  UnitTime_in_Megayears, 
  G,
  Hubble,
  a0, ar;

extern int Ntrees; // Number of trees per file.
extern int NtotGals; // Number of galaxies per file.
extern int *GalsForTree;

extern int Snaplistlen;
extern double ZZ[ABSOLUTEMAXSNAPS];
extern double AA[ABSOLUTEMAXSNAPS];
extern double Age[ABSOLUTEMAXSNAPS];

extern gsl_rng *random_generator;

extern int MAXSNAPS;

// Parameter File Variables //

extern int FirstFile;    // first and last file for processing 
extern int LastFile;

extern int NumGals;     // Total number of galaxies stored for current tree 
extern int MaxGals;     // Maximum number of galaxies allowed for current tree  

extern int LastSnapShotNr;

extern char GridOutputDir[512];
extern char GalaxiesInputDir[512];
extern char FileNameGalaxies[512];
extern char FileNameMergedGalaxies[512];
extern char SimulationDir[512];
extern char TreeName[512];
extern char DiffuseDir[512];
extern char FileWithSnapList[512];

extern int TotHalos;
extern int totNHalos; // Total number of halos in a file.
extern int TotGalaxies[ABSOLUTEMAXSNAPS];
extern int *TreeNgals[ABSOLUTEMAXSNAPS];

extern int *FirstHaloInSnap;

extern int *TreeNHalos;
extern int *TreeFirstHalo;

#ifdef MPI
extern int ThisTask, NTask, nodeNameLen;
extern char *ThisNode;
#endif

extern double Omega;
extern double OmegaLambda;
extern double PartMass;
extern double BoxSize;
extern int GridSize;
extern int GridSnap;
extern int NGrid;
extern int LastOutputSnap;


extern double Hubble_h;
extern double EnergySNcode, EnergySN;
extern double EtaSNcode, EtaSN;

// recipe flags 
extern int ReionizationOn;
extern int SupernovaRecipeOn;
extern int DiskInstabilityOn;
extern int AGNrecipeOn;
extern int SFprescription;

// recipe parameters 
extern double RecycleFraction;
extern double Yield;
extern double FracZleaveDisk;
extern double ReIncorporationFactor;
extern double ThreshMajorMerger;
extern double BaryonFrac;
extern double SfrEfficiency;
extern double FeedbackReheatingEpsilon;
extern double FeedbackEjectionEfficiency;
extern double RadioModeEfficiency;
extern double QuasarModeEfficiency;
extern double BlackHoleGrowthRate;
extern double Reionization_z0;
extern double Reionization_zr;
extern double ThresholdSatDisruption;

extern int OutputMode;
extern int LowSnap;
extern int HighSnap;
extern int ListOutputSnaps[ABSOLUTEMAXSNAPS];
extern int ListOutputGrid[ABSOLUTEMAXSNAPS];

extern int Verbose;

// Input parameters for the photon prescription // 
extern double SourceEfficiency;
extern int PhotonPrescription;
extern int fescPrescription;
extern double fesc; 

extern double MH_min;
extern double MH_max;
extern double alpha;
extern double beta;
extern double kappa;
extern double delta;
extern int HaloPartCut;

extern int use_sage;
#endif
