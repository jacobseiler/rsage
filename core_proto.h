#include "core_allvars.h"

size_t myfread(void  *ptr,  size_t  size,  size_t  nmemb,  FILE *stream);
size_t myfwrite(void  *ptr,  size_t  size,  size_t  nmemb,  FILE *stream);
int myfseek(FILE *stream, long offset, int whence);

void construct_galaxies(int halonr, int tree, int filenr);
void evolve_galaxies(int halonr, int ngal, int tree, int filenr);
int  join_galaxies_of_progenitors(int halonr, int nstart);
void init(void);
int32_t init_grid(void);
void set_units(void);

void load_tree_table(int filenr);
void load_tree(int filenr, int nr);
void save_galaxies(int filenr, int tree);
void save_merged_galaxies(int MergedNr, int filenr);

void prepare_galaxy_for_output(int filenr, int tree, struct GALAXY *g, struct GALAXY_OUTPUT *o);

void free_galaxies_and_tree(void);
void free_grid_arrays(struct GALAXY *g);
void free_tree_table(void);
int32_t free_grid(void);
void malloc_grid_arrays(struct GALAXY *g);
void print_allocated(void);

void read_parameter_file(char *fname);
void *mymalloc(size_t n);
void *myrealloc(void *p, size_t n); 
void myfree(void *p);
void myexit(int signum);

void finalize_galaxy_file(int filenr);
void write_gridarray(struct GALAXY *g, FILE *fp);
void finalize_merged_galaxy_file(int filenr);

void starformation_and_feedback(int p, int centralgal, double time, double dt, int halonr, int step, int tree, int ngal);
void add_galaxies_together(int t, int p);
void init_galaxy(int p, int halonr);
double infall_recipe(int centralgal, int ngal, double Zcurr, int halonr);
void add_infall_to_hot(int centralgal, double infallingGas, double dt);
double cooling_recipe(int centralgal, double dt);
void cool_gas_onto_galaxy(int centralgal, double coolingGas);
void reincorporate_gas(int centralgal, double dt);
double estimate_merging_time(int prog, int mother_halo, int ngal);
void deal_with_galaxy_merger(int p, int merger_centralgal, int centralgal, double time, double dt, int halonr, int step, int tree, int ngal);
void add_galaxy_to_merger_list(int p);
double dmax(double x, double y);
double do_reionization(int centralgal, double Zcurr, int ReturnMfilt);
double do_myreionization(int centralgal, double Zcurr, double *Mfilt);
double do_AGN_heating(double coolingGas, int centralgal, double dt, double x, double rcool);
void collisional_starburst_recipe(double mass_ratio, int merger_centralgal, int centralgal, double time, double dt, int halonr, int mode, int step, int tree, int ngal);
void update_from_star_formation(int p, double stars, double dt, int step, bool ismerger, int tree, int ngal);
void update_from_feedback(int p, int centralgal, double reheated_mass, double ejected_mass, double metallicity);
void make_bulge_from_burst(int p);
void grow_black_hole(int merger_centralgal, double mass_ratio, int32_t step);
void check_disk_instability(int p, int centralgal, int halonr, double time, double dt, int step, int tree, int ngal);

void strip_from_satellite(int halonr, int centralgal, int gal);
void disrupt_satellite_to_ICS(int centralgal, int gal, int tree);
void quasar_mode_wind(int gal, float BHaccrete, int32_t step);

double get_metallicity(double gas, double metals);
double get_virial_velocity(int halonr);
double get_virial_radius(int halonr);
double get_virial_mass(int halonr);
double get_disk_radius(int halonr, int p);

void read_output_snaps(void);
void read_snap_list(void);
void read_cooling_functions(void);
double get_metaldependent_cooling_rate(double logTemp, double logZ);
double get_rate(int tab, double logTemp);

double time_to_present(double z);
double integrand_time_to_present(double a, void *param);

double metallicity_dependent_star_formation(int p);
double Z_dependent_SF(float lower_limit, float upper_limit, float Sigma_c0, float Xi, float gamma);
double integrand_Z_dependent_SF(double q, void *p);

void grid_update_mass_metals(double mass, double metallicity, int GridPos, int SnapNum);
void grid_update_mass_metals_mass(double mass, double mass_metals, int GridPos, int SnapNum);

void update_grid_array(int p, int halonr, int steps_completed, int centralgal);
void calculate_photons(float SFR, float Z, float *Ngamma_HI, float *Ngamma_HeI, float *Ngamma_HeII);


void do_previous_SN(int p, int centralgal, int halonr, double dt);
void do_contemporaneous_SN(int p, int centralgal, double dt, double *stars, double *reheated_mass, double *mass_metals_new, double *mass_stars_recycled, double *ejected_mass);
void do_IRA_SN(int p, int centralgal, double *stars, double *reheated_mass, double *mass_metals_new, double *mass_stars_recycled, double *ejected_mass);
void do_previous_recycling(int p, int centralgal, int step, double dt); 
void calculate_Delta_Eta(double m_low, double m_high, double *Delta_Eta, double *Delta_m);
double calculate_reheated_mass(double Delta_Eta, double stars, double Vmax);
double calculate_reheated_energy(double Delta_Eta, double stars, double Vmax);
double calculate_ejected_mass(double *reheated_mass, double reheated_energy, double Vvir);
double calculate_coreburning(double t);
void update_from_SN_feedback(int p, int centralgal, double reheated_mass, double ejected_mass, double mass_stars_recycled, double mass_metals_new, double dt);
void  update_stars_array(int p, double stars, double dt, int tree, int step, int ngal);
