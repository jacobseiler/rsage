#include "core_allvars_grid.h"

int32_t init(void);
int32_t init_grid(struct GRID_STRUCT *grid);
int32_t free_grid(void);
void set_units(void);
void read_snap_list(void);
void determine_fesc_constants(void);

void read_parameter_file(char *fname);

int32_t load_gals(char *fname);
int load_meraxes_halos(int snapnum);
void load_merged_gals(char *fname);
void free_gals(void);
void free_meraxes_halos(void);

double time_to_present(double z);
double integrand_time_to_present(double a, void *param);
void estimate_grid_memory(void);
void estimate_gal_memory(int NtotGals);

int32_t update_grid_properties(int32_t filenr);
void update_meraxes_grid_properties(int p, int GridNr);
int32_t update_quasar_tracking(int64_t gal_idx, int32_t snapshot_idx);
void count_grid_properties(struct GRID_STRUCT *count_grid);
void normalize_photon(int GridNr); 
void normalize_slope_photons(int GridNr);
void calculate_photons(float SFR, float Z, float *Ngamma_HI, float *Ngamma_HeI, float *Ngamma_HeII);
float calculate_fesc(int p, int i, int filenr);
double get_metallicity(double gas, double metals);
#ifdef MPI
struct GRID_STRUCT *MPI_sum_grids(void);
#endif

const char* getfield(char* line, int num);

int32_t save_grid(struct GRID_STRUCT *save_grid);
void save_redshift(void);

void load_halos(int filenr);

void calculate_halomass(void);

void *mymalloc(size_t n);
void myfree(void *p);
