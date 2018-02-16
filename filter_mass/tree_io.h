#ifndef TREE_IO_H 
#define TREE_IO_H 

#include <stdint.h>

#define MAXLEN 1024

// Structs //

struct HALO_STRUCT 
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
};

typedef struct HALO_STRUCT *Halo_t;

// Proto-Types //

int32_t load_tree_table(int32_t filenr, char *treedir, char *treename, int32_t *Ntrees, int32_t *totNHalos, int32_t **TreeNHalos);
int32_t allocate_array_memory(int32_t totNHalos, int64_t **HaloID, float **Mfilt);   
int32_t free_memory(int32_t **TreeNHalos, int64_t **HaloID, float **Mfilt);   
int32_t load_halos(int32_t treenr, int32_t NHalos_ThisTree, Halo_t *Halos, int32_t *NHalos_ThisSnap);      
#endif
