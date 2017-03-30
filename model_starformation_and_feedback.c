#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "core_allvars.h"
#include "core_proto.h"


void do_previous_SN(int p, int centralgal, int halonr)
{
  double t, N_SFH_local, metallicity;
  double reheated_mass_previous, reheated_energy_previous, ejected_mass_previous, ejected_metals_previous, mass_stars_nova_previous;

  if (N_SFH == -1) // If we have set to do dynamic calculation for N_SFH.
  {
    while (t < 40)
    {
      t += (Age[Gal[p].SnapNum - N_SFH - 1] - Age[Gal[p].SnapNum - N_SFH]) * UnitTime_in_Megayears / Hubble_h;
      ++N_SFH_local;

    }
    --N_SFH_local; // If the time between snapshots is > 40Myr then N_SFH_local = 1 from the above loop but should be 0.  Hence need to subtract 1 off.
  }  else
  {
    N_SFH_local = N_SFH; // Otherwise set the value to the user defined one.
  }
 
  if(SupernovaRecipeOn == 1) // Calculating the ejected mass due to previous star formation episodes.
  {  
    calculate_previous_SN(p, N_SFH_local, &reheated_mass_previous, &reheated_energy_previous, centralgal, &ejected_metals_previous, &mass_stars_nova_previous); // Calculate the amount of energy and mass reheated from old SN events.
    calculate_ejected_mass(p, centralgal, &reheated_mass_previous, reheated_energy_previous, &ejected_mass_previous); // Calculate the amount of mass ejected from old SN events.
//    calculate_previous_ejected_metals(p, N_SFH_local, &ejected_metals_previous, &mass_stars_nova_previous); // Calculate the mass of metals added to the ISM from old SN (metal enrichment).  Also calcualte the mass from old SN that is added to the cold ISM (stellar mass recyling). 

    if(N_SFH_local == 0 && (reheated_mass_previous != 0 || reheated_energy_previous != 0 || ejected_mass_previous != 0 || ejected_metals_previous != 0 || mass_stars_nova_previous != 0))
    {
      fprintf(stderr, "Wanted to only use instantaneous SN feedback but there was non-zero values from previous episodes of SN.\n");
      fprintf(stderr, "Previous_ejected = %.4e \t Previous_Reheated = %.4e \t previous_metals = %.4e \t mass_stars_nova_previous = %.4e\n", ejected_mass_previous, reheated_mass_previous, ejected_metals_previous, mass_stars_nova_previous);
      exit(EXIT_FAILURE);
    }

 
  /*  
  printf("\n");
  printf("Reheated_Mass = %.4e \t Reheated_Energy = %.4e \t Ejected_Mass = %.4e \t Ejected_Metals = %.4e \t Mass_Stars_Nova = %.4e\n", reheated_mass_previous, reheated_energy_previous, ejected_mass_previous, ejected_metals_previous, mass_stars_nova_previous); 
  printf("Gal[p].ColdGas = %.4e \t Gal[p].HotGas = %.4e \t Gal[p].MetalsColdGas = %.4e \t StellarMass = %.4e\n",  Gal[p].ColdGas, Gal[p].HotGas, Gal[p].MetalsColdGas, Gal[p].StellarMass);
  */   
    metallicity = get_metallicity(Gal[p].ColdGas, Gal[p].MetalsColdGas);
    Gal[p].ColdGas += mass_stars_nova_previous; // Add back the mass of previous SN episodes directly to the cold ISM.
    Gal[p].MetalsColdGas += metallicity * mass_stars_nova_previous;
    Gal[p].StellarMass -= mass_stars_nova_previous;
    Gal[p].MetalsStellarMass -= metallicity * mass_stars_nova_previous;

    if(reheated_mass_previous > Gal[p].ColdGas)
    	reheated_mass_previous = Gal[p].ColdGas;
	 
    metallicity = get_metallicity(Gal[p].ColdGas, Gal[p].MetalsColdGas); // Need to recalculate metallicity as we have added more metals from previous SN episodes.
 
    // Update the amount of cold gas from the amount reheated and update the hot gas from the amount ejected. 
    update_from_feedback(p, centralgal, reheated_mass_previous, ejected_mass_previous, metallicity);
    Gal[p].PreviousReheatedMass[Gal[p].SnapNum] = reheated_mass_previous;

    // Add in the metals from previous SN. 
    if(Gal[p].ColdGas > 1.0e-8)
    {
      Gal[p].MetalsColdGas += ejected_metals_previous; 
      Gal[centralgal].MetalsHotGas += ejected_metals_previous; 
    }
    else
      Gal[centralgal].MetalsHotGas += ejected_metals_previous;
  /*	
  printf("AFTER UPDATING\n");
  printf("Reheated_Mass = %.4e \t Reheated_Energy = %.4e \t Ejected_Mass = %.4e \t Ejected_Metals = %.4e \t Mass_Stars_Nova = %.4e\n", reheated_mass_previous, reheated_energy_previous, ejected_mass_previous, ejected_metals_previous, mass_stars_nova_previous); 
  printf("Gal[p].ColdGas = %.4e \t Gal[p].HotGas = %.4e \t Gal[p].MetalsColdGas = %.4e \t StellarMass = %.4e\n",  Gal[p].ColdGas, Gal[p].HotGas, Gal[p].MetalsColdGas, Gal[p].StellarMass);
  */
  } 
}

void starformation_and_feedback(int p, int centralgal, double time, double dt, int halonr, int step)
{
  double reff, tdyn, strdot, stars, fac, metallicity;
  double cold_crit;
  double FracZleaveDiskVal;
  double reheated_mass_current = 0.0, reheated_energy_current = 0.0, ejected_mass_current = 0.0, ejected_metals_current = 0.0, mass_stars_nova_current = 0.0; 

  // Initialise variables
  strdot = 0.0;

  // star formation recipes 
  if(SFprescription == 0)
  {
    // we take the typical star forming region as 3.0*r_s using the Milky Way as a guide
    reff = 3.0 * Gal[p].DiskScaleRadius;
    tdyn = reff / Gal[p].Vvir;

    // from Kauffmann (1996) eq7 x piR^2, (Vvir in km/s, reff in Mpc/h) in units of 10^10Msun/h 
    cold_crit = 0.19 * Gal[p].Vvir * reff;
    if(Gal[p].ColdGas > cold_crit && tdyn > 0.0)
    	strdot = SfrEfficiency * (Gal[p].ColdGas - cold_crit) / tdyn;
    else
    	strdot = 0.0;
  }
  else
  {
    printf("No star formation prescription selected!\n");
    ABORT(0);
  }

  stars = strdot * dt;
  if(stars < 0.0)
    stars = 0.0;
   
  if(SupernovaRecipeOn == 1)
  {
    calculate_current_SN(p, stars, &reheated_mass_current, &reheated_energy_current, &ejected_metals_current, &mass_stars_nova_current);
    calculate_ejected_mass(p, centralgal, &reheated_mass_current, reheated_energy_current, &ejected_mass_current);	
  }
  else
  {
    reheated_mass_current = 0.0;
    mass_stars_nova_current = 0.0;
    ejected_metals_current = 0.0;
  }

  /*
  if(Gal[p].StellarMass < 0.1) 
  {
  printf("\n");
  printf("Reheated_Mass = %.4e \t Reheated_Energy = %.4e \t Ejected_Mass = %.4e \t Ejected_Metals = %.4e \t Mass_Stars_Nova = %.4e \t reheated_mass_\n", reheated_mass_current, reheated_energy_current, ejected_mass_current, ejected_metals_current, mass_stars_nova_current); 
  printf("Gal[p].ColdGas = %.4e \t Gal[p].HotGas = %.4e \t Gal[p].MetalsColdGas = %.4e \t StellarMass = %.4e \t stars = %.4e\n",  Gal[p].ColdGas, Gal[p].HotGas, Gal[p].MetalsColdGas, Gal[p].StellarMass, stars);
  printf("\n");
  }
  */
 
  assert(reheated_mass_current >= 0.0);

  if(Gal[p].StellarMass < 0.0)
    Gal[p].StellarMass = 0.0;
  if (Gal[p].MetalsStellarMass < 0.0)
    Gal[p].MetalsStellarMass = 0.0;

  // cant use more cold gas than is available! so balance SF and feedback 
  if((stars + reheated_mass_current) > Gal[p].ColdGas && (stars + reheated_mass_current) > 0.0)
  {
    fac = Gal[p].ColdGas / (stars + reheated_mass_current);
    stars *= fac;
    reheated_mass_current *= fac;
  }

  metallicity = get_metallicity(Gal[p].ColdGas, Gal[p].MetalsColdGas);
  Gal[p].ColdGas += mass_stars_nova_current; // Add back the mass of this SN episode directly to the cold ISM.
  Gal[p].MetalsColdGas += metallicity * mass_stars_nova_current;
  Gal[p].StellarMass -= mass_stars_nova_current;
  Gal[p].MetalsStellarMass -= metallicity * mass_stars_nova_current;

  Gal[p].SNStars[Gal[p].SnapNum] += stars;

//  printf("reheated_mass_previous = %.4e \t reheated_mass_current = %.4e \t ejected_mass_previous = %.4e \t ejected_mass_current = %.4e \t ejected_energy_previous = %.4e \t ejected_energy_current = %.4e\n", reheated_mass_previous, reheated_mass_current, ejected_mass_previous, ejected_mass_current, reheated_energy_previous, reheated_energy_current); 
  
  // update the star formation rate 
  Gal[p].SfrDisk[step] += stars / dt;
  Gal[p].SfrDiskColdGas[step] = Gal[p].ColdGas;
  Gal[p].SfrDiskColdGasMetals[step] = Gal[p].MetalsColdGas;

  // update for star formation 
  metallicity = get_metallicity(Gal[p].ColdGas, Gal[p].MetalsColdGas);
  update_from_star_formation(p, stars, metallicity, mass_stars_nova_current);

  // recompute the metallicity of the cold phase
  metallicity = get_metallicity(Gal[p].ColdGas, Gal[p].MetalsColdGas);

  // update from SN feedback 
  update_from_feedback(p, centralgal, reheated_mass_current, ejected_mass_current, metallicity);

  // check for disk instability
  if(DiskInstabilityOn)
    check_disk_instability(p, centralgal, halonr, time, dt, step);

  // New metals injected form this snapshot's SN episode.  
  if(Gal[p].ColdGas > 1.0e-8)
  {
    Gal[p].MetalsColdGas += ejected_metals_current; 
    Gal[centralgal].MetalsHotGas += ejected_metals_current; 
  }
  else
    Gal[centralgal].MetalsHotGas += Yield * stars;
    // Gal[centralgal].MetalsEjectedMass += Yield * stars;
   
}

void update_from_star_formation(int p, double stars, double metallicity, double mass_stars_nova)
{
  // update gas and metals from star formation 
  Gal[p].ColdGas -= (stars - mass_stars_nova);
  Gal[p].MetalsColdGas -= metallicity * (stars - mass_stars_nova);
  Gal[p].StellarMass += (stars - mass_stars_nova);
  Gal[p].MetalsStellarMass += metallicity * (stars - mass_stars_nova);
  if (Gal[p].ColdGas < 0.0)
    Gal[p].ColdGas = 0.0;
  if (Gal[p].MetalsColdGas < 0.0)
    Gal[p].MetalsColdGas = 0.0;
}

void update_from_feedback(int p, int centralgal, double reheated_mass, double ejected_mass, double metallicity)
{
  double metallicityHot;

	assert(!(reheated_mass > Gal[p].ColdGas && reheated_mass > 0.0));

  if(SupernovaRecipeOn == 1)
  {
    Gal[p].ColdGas -= reheated_mass;
    Gal[p].MetalsColdGas -= metallicity * reheated_mass;

    Gal[centralgal].HotGas += reheated_mass;
    Gal[centralgal].MetalsHotGas += metallicity * reheated_mass;

    if(ejected_mass > Gal[centralgal].HotGas)
      ejected_mass = Gal[centralgal].HotGas;
    metallicityHot = get_metallicity(Gal[centralgal].HotGas, Gal[centralgal].MetalsHotGas);

    Gal[centralgal].HotGas -= ejected_mass;
    Gal[centralgal].MetalsHotGas -= metallicityHot * ejected_mass;
    Gal[centralgal].EjectedMass += ejected_mass;
    Gal[centralgal].MetalsEjectedMass += metallicityHot * ejected_mass;

    Gal[p].OutflowRate += reheated_mass;   

    //Gal[centralgal].deltaEjectedMass[Gal[centralgal].SnapNum] += ejected_mass;
    //Gal[centralgal].deltaMetalsEjectedMass[Gal[centralgal].SnapNum] += ejected_mass * metallicityHot;

  }

}

void calculate_previous_SN(int p, int SnapHistory, double *reheated_mass, double *reheated_energy, int centralgal, double *ejected_metals, double *mass_stars_nova)
{
  int i;
  double m_low, m_high, t;
  double a = 0.7473, b = -2.6979, c = -4.7659, d = 0.5934;
  double Delta_Eta, Delta_m;
  *reheated_mass = 0.0; 
  *reheated_energy = 0.0;
  *ejected_metals = 0.0;
  *mass_stars_nova = 0.0;


  for(i = Gal[p].SnapNum - SnapHistory; i < Gal[p].SnapNum - 1; ++i)
  {
    t = (Age[i] - Age[Gal[p].SnapNum]) * UnitTime_in_Megayears / Hubble_h; // Time between snapshot i and the current snapshot (in Myr). 
    m_low = pow(10, a/log10(t) + b * exp(c/log10(t)) + d); // Low mass limit for stars that can go nova between snapshots i and current snapshot. In units of Msun.
    if(m_low < 8)
	m_low = 8; // Force a minimum m_low of 8 Msun.
    if (m_low > 120)
	Delta_Eta = 0;
    else
    {
      m_high = 120;
      if (i - Gal[p].SnapNum == 2)
      {
      	Delta_Eta = 0.00960467; // Chabrier.
	//Delta_Eta = 0.00430701; // Salpeter.

      	Delta_m = 0.254446; // Chabrier.
	//Delta_m = 0.114107; // Salpeter. 
	
      }
      else if(i - Gal[p].SnapNum == 3)
      {
	Delta_Eta = 0.0134431; // Chabrier.
	//Delta_Eta = 0.00602864; // Salpeter.

	Delta_m = 0.294533; // Chabrier. 
	//Delta_m = 0.132085; // Salpeter. 
      }
      else if(i - Gal[p].SnapNum == 4)
      {
	Delta_Eta = 0.0165724; // Chabrier.
	//Delta_Eta = 0.00743198; // Salpeter.

	Delta_m = 0.321102; // Chabrier.
	//Delta_m = 0.14417; // Salpeter.

      }
      else
      {
      	Delta_Eta = calculate_Delta_Eta(m_low, m_high,0); 
      	Delta_m = calculate_Delta_Eta(m_low, m_high,1);
      } 
    }

    //printf("History = %.4e \t Current = %.4e\n", Gal[centralgal].VmaxHistory[Gal[centralgal].SnapNum], Gal[centralgal].Vmax);

    //double epsilon_energy = alpha_energy * (0.5 + pow(Gal[centralgal].VmaxHistory[Gal[centralgal].SnapNum]/V_energy, -beta_energy));
    //double epsilon_mass = alpha_mass * (0.5 + pow(Gal[centralgal].VmaxHistory[Gal[centralgal].SnapNum]/V_mass, -beta_mass));

    double epsilon_energy = alpha_energy * (0.5 + pow(Gal[centralgal].Vmax/V_energy, -beta_energy));
    double epsilon_mass = alpha_mass * (0.5 + pow(Gal[centralgal].Vmax/V_mass, -beta_mass));


    if (epsilon_mass > epsilon_mass_max)
	epsilon_mass = epsilon_mass_max;

    epsilon_mass = FeedbackReheatingEpsilon;
    epsilon_energy = FeedbackEjectionEfficiency;

    *reheated_mass += Delta_Eta/Eta_SNII * Gal[p].SNStars[i] * epsilon_mass; 
    *reheated_energy += Delta_Eta * (Gal[p].SNStars[i]*1.0e10/Hubble_h) * EtaSN * EnergySN * epsilon_energy; 

    *ejected_metals += Delta_m / m_SNII * Gal[p].SNStars[i] * Yield; 
    *mass_stars_nova += Delta_m * Gal[p].SNStars[i]; // Mass of stars that have gone nova.  Note we convert to code units of 1.0e10 Msun/h.
    Gal[p].SNStars[i] -= *mass_stars_nova; // The stars in snapshot i have gone nova, hence they cannot be used for further SN. 
  } 
}

/*
void calculate_previous_ejected_metals(int p, int SnapHistory, double *ejected_metals, double *mass_stars_nova)
{

  int i;
  double m_low, m_high, t;
  double a = 0.7473, b = -2.6979, c = -4.7659, d = 0.5934;
  double Delta_m;
  *ejected_metals = 0.0;
  *mass_stars_nova = 0.0;

  for(i = Gal[p].SnapNum - SnapHistory; i < Gal[p].SnapNum - 1; ++i)
  {
    t = (Age[i] - Age[Gal[p].SnapNum]) * UnitTime_in_Megayears / Hubble_h; // Time between snapshot i and the current snapshot (in Myr). 
    m_low = pow(10, a/log10(t) + b * exp(c/log10(t)) + d); // Low mass limit for stars that can go nova between snapshots i and current snapshot. In units of Msun.
    if(m_low < 8)
	m_low = 8; // Force a minimum m_low of 8 Msun.
    if (m_low > 120)
	Delta_m = 0;
    else
    {
      m_high = 120;
      if (i - Gal[p].SnapNum == 2)
      {
      	//Delta_m = 0.254446; // Chabrier.
	Delta_m = 0.114107; // Salpeter. 
      } 
      else if(i - Gal[p].SnapNum == 3)
      {
	//Delta_m = 0.294533; // Chabrier. 
	Delta_m = 0.132085; // Salpeter. 
      }
      else if(i - Gal[p].SnapNum == 4)
      {
	//Delta_m = 0.321102; // Chabrier.
	Delta_m = 0.14417; // Salpeter.
      }
      else
      	Delta_m = calculate_Delta_Eta(m_low, m_high,1); 
    }
	
    *ejected_metals += Delta_m / m_SNII * Gal[p].SNStars[i] * Yield; 
    *mass_stars_nova += Delta_m * Gal[p].SNStars[i]; // Mass of stars that have gone nova.  Note we convert to code units of 1.0e10 Msun/h. 
  }

}
*/


void calculate_current_SN(int p, double stars, double *reheated_mass, double *reheated_energy, double *ejected_metals, double *mass_stars_nova)
{ 
  double m_low, m_high, t;
  double a = 0.7473, b = -2.6979, c = -4.7659, d = 0.5934;
  double Delta_Eta, Delta_m;
  *reheated_mass = 0.0;
  *reheated_energy = 0.0; 
  *ejected_metals = 0.0;
  *mass_stars_nova = 0.0;
 
  t = (Age[Gal[p].SnapNum] - Age[Gal[p].SnapNum + 1]) * UnitTime_in_Megayears / Hubble_h; // Time between current snapshot and the next one (in Myr). 
  m_low = pow(10, a/log10(t) + b * exp(c/log10(t)) + d); // Low mass limit for stars that can go nova between snapshots i and current snapshot. In units of Msun.
  m_low = 1;
  if(m_low < 8)
	m_low = 8; // Force a minimum m_low of 8 Msun.
  if (m_low > 120)
	Delta_Eta = 0;
  else
  {
    m_high = 120;

    // Fraction of stars that go nova in the current snapshot.
    Delta_Eta = 0.004784; // Chabrier.
    //Delta_Eta = 0.00214473; // Salpeter.

    // Mass fraction of stars that go nova in the current snapshot.
    Delta_m = 0.18301; // Chabrier.
    //Delta_m = 0.0820719; // Salpeter.

  }

  double epsilon_energy = alpha_energy * (0.5 + pow(Gal[p].Vmax/V_energy, -beta_energy));
  double epsilon_mass = alpha_mass * (0.5 + pow(Gal[p].Vmax/V_mass, -beta_mass));

  if (epsilon_mass > epsilon_mass_max)
      epsilon_mass = epsilon_mass_max;

  epsilon_mass = FeedbackReheatingEpsilon;
  epsilon_energy = FeedbackEjectionEfficiency;

  *reheated_mass += Delta_Eta/Eta_SNII * stars * epsilon_mass; // Only doing these proportional to the stars formed in the current snapshot. 
  *reheated_energy += Delta_Eta * (stars * 1.0e10 / Hubble_h) * EnergySN * EtaSN * epsilon_energy; 
  
  *ejected_metals = Delta_m / m_SNII * stars * Yield;
  *mass_stars_nova = Delta_m * stars;

}

/*
void calculate_current_ejected_metals(int p, double stars, double *ejected_metals, double *mass_stars_nova) 
{ 
  double m_low, m_high, t;
  double a = 0.7473, b = -2.6979, c = -4.7659, d = 0.5934;
  double Delta_m;
  *ejected_metals = 0.0;
  *mass_stars_nova = 0.0;
 
  t = (Age[Gal[p].SnapNum] - Age[Gal[p].SnapNum + 1]) * UnitTime_in_Megayears / Hubble_h; // Time between current snapshot and the next one (in Myr). 
  m_low = pow(10, a/log10(t) + b * exp(c/log10(t)) + d); // Low mass limit for stars that can go nova between snapshots i and current snapshot. In units of Msun.
  if(m_low < 8)
	m_low = 8; // Force a minimum m_low of 8 Msun.
  if (m_low > 120)
	Delta_m = 0;
  else
  {
    m_high = 120;
  }

  *ejected_metals = Delta_m / m_SNII * stars * Yield;
  *mass_stars_nova = Delta_m * stars;
 // printf("Delta_m = %.4e \t stars = %.4e \t mass_stars_nova = %.4e\n", Delta_m, stars, *mass_stars_nova);
}
*/

double calculate_Delta_Eta(double m_low, double m_high, int function)
{
// Function is a flag that controls which function to integrate.
// 0: Integrate phi(m).
// 1: Integrate m*phi(m).
#define WORKSIZE 1000
  gsl_function F;
  gsl_integration_workspace *workspace;
  double result, abserr;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  if (function == 0)
  	F.function = &integrand_Delta_Eta;
  else if (function == 1)
  	F.function = &integrand_Delta_m;

  gsl_integration_qag(&F, m_low, m_high, 1e-5, 0.15, 
    WORKSIZE, GSL_INTEG_GAUSS21, workspace, &result, &abserr);

  gsl_integration_workspace_free(workspace);

  return result;
}

double integrand_Delta_Eta(double m, void *param)
{
  return IMF_norm*pow(m, -2.35); 
}

double integrand_Delta_m(double m, void *param)
{
  return m*IMF_norm*pow(m, -2.35); 
}

void calculate_ejected_mass(int p, int centralgal, double *reheated_mass, double reheated_energy, double *ejected_mass)
{
  if(Gal[centralgal].Vvir > 0.0)
  {
    double Delta_Ehot = (0.5 * (*reheated_mass * 1.0e10 / Hubble_h * SOLAR_MASS / 1.0e3) * (Gal[centralgal].Vvir * 1.0e3) * (Gal[centralgal].Vvir * 1.0e3)) * 1.0e7; // Change in the thermal energy of the hot resevoir IN ERG.
//    printf("reheated_energy = %.4e \t Delta_Ehot = %.4e \t Difference in energy = %.4e \t Gal[centralgal].Vvir = %.4e\n", reheated_energy, Delta_Ehot, reheated_energy - Delta_Ehot, Gal[centralgal].Vvir);
    if(reheated_energy > Delta_Ehot)
    {
    	*ejected_mass = (reheated_energy - Delta_Ehot) * 1.0e-7/(0.5 * Gal[centralgal].Vvir * 1.0e3 * Gal[centralgal].Vvir * 1.0e3); // Balance between the excess thermal energy and the thermal energy of the hot gas.
        *ejected_mass = *ejected_mass * 1.0e3 / SOLAR_MASS / 1.0e10 * Hubble_h; // Convert back from kg to 1.0e10*Msun/h. 
    }
    else
    {
	*ejected_mass = 0.0;
//	printf("Reheated_mass before = %.4e ", *reheated_mass);
        *reheated_mass = reheated_energy * 1.0e-7 / (0.5 * Gal[centralgal].Vvir * 1.0e3 * Gal[centralgal].Vvir * 1.0e3);
	*reheated_mass = *reheated_mass * 1.0e3 / SOLAR_MASS / 1.0e10 * Hubble_h; // Convert from kg to 1.0e10*Msun/h.
//	printf("Reheated_mass after = %.4e\n", *reheated_mass);
    }
  } 
					
  else
    *ejected_mass = 0.0;
		
  if(*ejected_mass < 0.0)
    *ejected_mass = 0.0;

}

