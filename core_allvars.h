#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <stdbool.h>
#include <gsl/gsl_rng.h>
#include "core_simulation.h"

#define ABORT(sigterm)                                                  \
do {                                                                \
  printf("Error in file: %s\tfunc: %s\tline: %i\n", __FILE__, __FUNCTION__, __LINE__); \
  myexit(sigterm);                                                \
} while(0)

#define  STEPS 10        // Number of integration intervals between two snapshots
#define  MAXGALFAC 1
#define  ALLOCPARAMETER 10.0
#define  MAX_NODE_NAME_LEN 50
#define  ABSOLUTEMAXSNAPS 1000

#define  MAXLEN      1024
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

#ifdef NDEBUG
#define XPRINT(EXP, ...)                                do{} while(0)
#else
#define XPRINT(EXP, ...)                                              \
    do { if (!(EXP)) {                                                  \
            printf("Warning in file: %s\tfunc: %s\tline: %d with expression `"#EXP"'\n", __FILE__, __FUNCTION__, __LINE__); \
            printf(__VA_ARGS__);                                        \
            fflush(stdout);                                             \
        } \
    } while (0)
#endif


struct GALAXY_OUTPUT
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

  // NOTE: Because we have moved to pointers for the galaxy-grid parameters
  // they are no longered stored explicitly in the galaxy output struct.
  // Instead we have to write them out individually (otherwise we would just be writing out pointers and not their values).

};

// This structure contains the properties used within the code
struct GALAXY
{
  int   SnapNum;
  int   Type;

  int   GalaxyNr;
  int   CentralGal;
  int   HaloNr;
  long long  MostBoundID;

  int   mergeType;  //0=none; 1=minor merger; 2=major merger; 3=disk instability; 4=disrupt to ICS
  int   mergeIntoID;
  int   mergeIntoSnapNum;
  float   dT;
  float Age;
  // (sub)halo properties
  float Pos[3];
  float Vel[3];
  int   Len;   
  float Mvir;
  float deltaMvir;
  float CentralMvir;
  float Rvir;
  float Vvir;
  float Vmax;

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
  float SfrDisk[STEPS];
  float SfrBulge[STEPS];
  float SfrDiskColdGas[STEPS];
  float SfrDiskColdGasMetals[STEPS];
  float SfrBulgeColdGas[STEPS];
  float SfrBulgeColdGasMetals[STEPS];

  // misc 
  float DiskScaleRadius;
  float MergTime;
  double Cooling;
  double Heating;
  float r_heat;
  float QuasarModeBHaccretionMass;
  float TimeOfLastMajorMerger;
  float TimeOfLastMinorMerger;
  float OutflowRate;
	float TotalSatelliteBaryons;

  // infall properties
  float infallMvir;
  float infallVvir;
  float infallVmax;

  int GridPos;
  int *GridHistory;
  float *GridStellarMass;
  float *GridSFR;
  float *GridZ;
  float *GridCentralGalaxyMass;
  double *MfiltGnedin;
  double *MfiltSobacchi;
  float *EjectedFraction;
  int *LenHistory;
  double *Stars;
  int *MergerHistory;

  double StellarAge_Numerator;
  double StellarAge_Denominator;
  double Total_SF_Time;
  double Total_Stars;

  int IsMerged;
  double GrandSum;
  int IsMalloced;
}
*Gal, *HaloGal, *MergedGal;


// auxiliary halo data
struct halo_aux_data   
{
  int DoneFlag;
  int HaloFlag;
  int NGalaxies;
  int FirstGalaxy;
}
*HaloAux;


struct PHOTO_GRID
{
  double *PhotoRate;
  double ReionRedshift;
}
*PhotoGrid;

struct GRID
{
  double EjectedMass[64]; 
  double MetalsEjectedMass[64];
  double StellarMass[64];
  double MetalsStellarMass[64]; 
  double Diffuse[64];
  int GalCount[64];
  double StellarAge[64];
  float Ionization_State[64];
  int ID;
}
*Grid;

struct GRID_OUTPUT
{
  double EjectedMass;
  double MetalsEjectedMass;
  double StellarMass;
  double Diffuse;
  double MfiltGnedin;

};

extern int    FirstFile;    // first and last file for processing 
extern int    LastFile;

extern int count_onehalo;

extern int    Ntrees;      // number of trees in current file 
extern int    NumGals;     // Total number of galaxies stored for current tree 
extern int    MaxGals;     // Maximum number of galaxies allowed for current tree  
extern int    FoF_MaxGals;

extern int    GalaxyCounter;     // unique galaxy ID for main progenitor line in tree

extern int    LastSnapShotNr;

extern char   OutputDir[512];
extern char   FileNameGalaxies[512];
extern char   TreeName[512];
extern char   SimulationDir[512];
extern char   IonizationDir[512];
extern char   FileWithSnapList[512];
extern char   PhotoionDir[512];
extern char   PhotoionName[512];
extern char   ReionRedshiftName[512];

extern int    TotHalos;
extern int    TotGalaxies[ABSOLUTEMAXSNAPS];
extern int    *TreeNgals[ABSOLUTEMAXSNAPS];

extern int    *FirstHaloInSnap;

extern int    *TreeNHalos;
extern int    *TreeFirstHalo;

#ifdef MPI
extern int ThisTask, NTask, nodeNameLen;
extern char *ThisNode;
#endif

extern int IMF;

extern double Omega;
extern double OmegaLambda;
extern double PartMass;
extern double BoxSize;
extern int GridSize;
extern double Hubble_h;
extern double EnergySNcode, EnergySN;
extern double EtaSNcode, EtaSN;
extern int use_tiamat;

// recipe flags 
extern int    ReionizationOn;
extern int    SupernovaRecipeOn;
extern int    DiskInstabilityOn;
extern int    AGNrecipeOn;
extern int    SFprescription;

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

extern int    ListOutputSnaps[ABSOLUTEMAXSNAPS];

extern double ZZ[ABSOLUTEMAXSNAPS];
extern double AA[ABSOLUTEMAXSNAPS];
extern double Age[ABSOLUTEMAXSNAPS];

extern int    MAXSNAPS;
extern int    NOUT;
extern int    Snaplistlen;

extern gsl_rng *random_generator;

extern int TreeID;
extern int FileNum;
extern int neg_cell;
extern int pos_cell;

extern double ejectedmass_total;
extern double metalsejectedmass_total;

extern int MergedNr;
extern int TotMerged; 
extern int *TreeNMergedgals;
extern int MaxMergedGals;     // Maximum number of galaxies allowed for current tree  

extern int zeromass_count;
extern int suppression_count;
extern int previous_tree;
extern int lowmass_halo;

extern double smallest_mass;

extern double count;

extern double IMF_norm;
extern double IMF_slope;
extern double Eta_SNII;
extern double m_SNII;
extern int IRA; 
extern int TimeResolutionSN;
extern int SN_Array_Len; 
extern int Time_SFH;

extern double alpha_energy;
extern double V_energy;
extern double beta_energy;
extern double alpha_mass;
extern double V_mass;
extern double beta_mass; 
extern double epsilon_mass_max;

extern int count_firstSF;
extern int count_notfirstSF;

extern int RescaleSN;
extern int mergedgal_mallocs;
extern int gal_mallocs;

extern int mergedgal_frees;
extern int gal_frees;
#endif  // #ifndef ALLVARS_H
