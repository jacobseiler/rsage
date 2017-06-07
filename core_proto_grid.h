#include "core_allvars_grid.h"

void init(void);
void init_grid(int GridNr, int ThisTask_GridNr);
void set_units(void);
void read_snap_list(void);

void read_parameter_file(char *fname);

void load_gals(char *fname);
void load_merged_gals(char *fname);
void free_gals(void);

double time_to_present(double z);
double integrand_time_to_present(double a, void *param);
void estimate_grid_memory(void);
void estimate_gal_memory(int NtotGals);

void update_grid_properties(int p, int merged, int GridNr);
void update_grid_halo(int totNHalos, int GridNr);
void update_grid_diffuse(int GridNr);
void update_grid_density(int GridNr);
void update_grid_nion_halo(int GridNr);
void count_grid_properties(int GridNr);
void normalize_photon(int GridNr); 
void normalize_slope_photons(int GridNr);
void calculate_photons(float SFR, float Z, float *Ngamma_HI, float *Ngamma_HeI, float *Ngamma_HeII);
double calculate_fesc(int p, int i);

const char* getfield(char* line, int num);

void save_grid(int GridNr);
void save_redshift(void);

void load_halos(int filenr);

void calculate_halomass(void);

void *mymalloc(size_t n);
void myfree(void *p);
