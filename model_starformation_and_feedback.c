#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <stdbool.h>

#include "core_allvars.h"
#include "core_proto.h"


void update_from_SN_feedback(int p, int centralgal, double reheated_mass, double ejected_mass, double mass_stars_recycled, double mass_metals_new)
{

  double metallicity, metallicityHot;

  Gal[p].StellarMass -= mass_stars_recycled; // The supernova remnants are instantly moved back into the cold ISM. 
  Gal[p].ColdGas += mass_stars_recycled; 

  if(Gal[p].ColdGas > 1.0e-10)
      Gal[p].MetalsColdGas += mass_metals_new; // ISM is enriched by new metals.
  else
      Gal[centralgal].MetalsHotGas += mass_metals_new;

  if(reheated_mass > Gal[p].ColdGas) // Just perform this check again to be sure.
      reheated_mass = Gal[p].ColdGas; // Because SF takes gas out of cold resevoir.

  if(ejected_mass > Gal[centralgal].HotGas)
      ejected_mass = Gal[centralgal].HotGas; 

  metallicity = get_metallicity(Gal[p].ColdGas, Gal[p].MetalsColdGas);

  Gal[p].ColdGas -= reheated_mass;
  Gal[p].MetalsColdGas -= reheated_mass * metallicity;

  Gal[centralgal].HotGas += reheated_mass;
  Gal[centralgal].MetalsHotGas += reheated_mass * metallicity;

  metallicityHot = get_metallicity(Gal[centralgal].HotGas, Gal[centralgal].MetalsHotGas);

  Gal[centralgal].HotGas -= ejected_mass;
  Gal[centralgal].MetalsHotGas -= ejected_mass * metallicityHot;
  Gal[centralgal].EjectedMass += ejected_mass;
  Gal[centralgal].MetalsEjectedMass += ejected_mass * metallicity;

  Gal[p].OutflowRate += reheated_mass;   

  if(Gal[p].ColdGas < 0.0) // Some final checks. 
    Gal[p].ColdGas = 0.0;
  if(Gal[p].MetalsColdGas < 0.0)
    Gal[p].MetalsColdGas = 0.0;
  if(Gal[p].StellarMass < 0.0)
    Gal[p].StellarMass = 0.0;
  if(Gal[centralgal].HotGas < 0.0)
    Gal[centralgal].HotGas = 0.0;
  if(Gal[p].MetalsHotGas < 0.0)
    Gal[p].MetalsHotGas = 0.0; 
 
}


void update_from_star_formation(int p, double stars, double dt, int step, bool ismerger) 
{

  double metallicity = get_metallicity(Gal[p].ColdGas, Gal[p].MetalsColdGas);
  if(!ismerger)
  {
    Gal[p].SfrDisk[step] += stars / dt;
    Gal[p].SfrDiskColdGas[step] = Gal[p].ColdGas;
    Gal[p].SfrDiskColdGasMetals[step] = Gal[p].MetalsColdGas;
  } else
  {
    
    Gal[p].SfrBulge[step] += stars / dt;
    Gal[p].SfrBulgeColdGas[step] += Gal[p].ColdGas;
    Gal[p].SfrBulgeColdGasMetals[step] += Gal[p].MetalsColdGas;
  
    Gal[p].BulgeMass += stars;
    Gal[p].MetalsBulgeMass += metallicity * stars;

  }
 
  Gal[p].SNStars[Gal[p].SnapNum] += stars;

  // update gas and metals from star formation 
  Gal[p].ColdGas -= stars;
  Gal[p].MetalsColdGas -= metallicity * stars;
  Gal[p].StellarMass += stars;
  Gal[p].MetalsStellarMass += metallicity * stars;
  if (Gal[p].ColdGas < 0.0)
    Gal[p].ColdGas = 0.0;
  if (Gal[p].MetalsColdGas < 0.0)
    Gal[p].MetalsColdGas = 0.0;
}


void starformation_and_feedback(int p, int centralgal, double time, double dt, int halonr, int step)
{
  double reff, tdyn, strdot, stars;
  double cold_crit;
  double FracZleaveDiskVal;
  double reheated_mass = 0.0, mass_metals_new = 0.0, mass_stars_recycled = 0.0, ejected_mass = 0.0;

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
    {
    	strdot = 0.0;
//	fprintf(stderr, "Gal[p].ColdGas = %.4e \t cold_crit = %.4e\n", Gal[p].ColdGas, cold_crit);
    }
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
    do_current_SN(p, centralgal, halonr, &stars, &reheated_mass, &mass_metals_new, &mass_stars_recycled, &ejected_mass);    
  } else
  { 
    if(stars > Gal[p].ColdGas) // We do this check in 'do_current_SN()' but if SupernovaRecipeOn == 0 then we don't do the check.
    {
      double factor = Gal[p].ColdGas / stars; 
      stars *= factor; 
    }
  }

  update_from_star_formation(p, stars, dt, step, false); 
  update_from_SN_feedback(p, centralgal, reheated_mass, ejected_mass, mass_stars_recycled, mass_metals_new);

  // check for disk instability
  if(DiskInstabilityOn)
    check_disk_instability(p, centralgal, halonr, time, dt, step);

}


// This function will calculate the fraction of stars formed in snapshot i that go nova in snapshot j.
// Will also calculate the mass fraction of stars formed in snapshot i that go nova in snapshot j.
// Based off equation (17) and (24) of Mutch et al. (2016).
//
// INPUT: The bounds for the integral (m_low and m_high). ## UNITS: Msun.
// 	: Pointers to modify the values for number and mass fraction of stars that go nova (*Delta_Eta and *Delta_m). ## UNITS: Msun^-1 (Delta_Eta) and unitless (Delta_m).
//
// OUTPUT: None but delta_Eta, the number fraction of stars, and delta_m, the mass fraction of stars, will be modified.

void calculate_Delta_Eta(double m_low, double m_high, double *Delta_Eta, double *Delta_m)
{
  *Delta_Eta = IMF_norm / (IMF_slope + 1.0) * (pow(m_high, IMF_slope + 1.0) - pow(m_low, IMF_slope + 1.0));
  *Delta_m = IMF_norm / (IMF_slope + 1.0) * (pow(m_high, IMF_slope + 2.0) - pow(m_low, IMF_slope + 2.0));
}

// This function determines the amount of mass reheated for a supernova event.
//
// INPUT: The number fraction of stars that have gone nova (Delta_Eta) ## UNITS: Msun^-1.
// 	: The mass of stars formed (stars). ## UNITS: 1.0e10Msun/h (code units).
// 	: The maximum velocity of the galaxy disc (Vmax). ## UNITS: km/s.
//
// OUTPUT: The amount of mass reheated by the supernova event. ## UNITS: 1.0e10Msun/h (code units). 

double calculate_reheated_mass(double Delta_Eta, double stars, double Vmax)
{

 // if(stars < 1e-10)
 //   return 0.0;

  double epsilon_mass = alpha_mass * (0.5 + pow(Vmax/V_mass, -beta_mass));
 
  if (epsilon_mass > epsilon_mass_max)
      epsilon_mass = epsilon_mass_max; // We enforce a maximum value for the mass loading factor. 

  //epsilon_mass = FeedbackReheatingEpsilon;


  double reheated_mass = Delta_Eta/Eta_SNII * stars * epsilon_mass;
  // fprintf(stderr, "Delta_Eta = %.4e \t Eta_SNII = %.4e \t stars = %.4e \t epsilon_mass = %.4e \t reheated_mass = %.4e\n", Delta_Eta, Eta_SNII, stars, epsilon_mass, reheated_mass);

  return reheated_mass; 

}

// This function determines the amount of energy injected from a supernova event. 
//
// INPUT: The number fraction of stars that have gone nova (Delta_Eta) ## UNITS: Msun^-1.
// 	: The mass of stars formed (stars). ## UNITS: 1.0e10Msun/h (code units).
// 	: The maximum velocity of the galaxy disc (Vmax). ## UNITS: km/s.
//
// OUTPUT: The amount of energy from a supernova event. ## UNITS: erg.  

double calculate_reheated_energy(double Delta_Eta, double stars, double Vmax)
{

//  if(stars < 1.0e-10)
//    return 0.0;

  double epsilon_energy = alpha_energy * (0.5 + pow(Vmax/V_energy, -beta_energy));
  //epsilon_energy = FeedbackEjectionEfficiency;
  //double reheated_energy = (Delta_Eta / 1.0e10 * Hubble_h) * stars * epsilon_energy * EnergySN * EtaSN; // Delta_Eta was in Msun so need to convert to 1e10/h Msun. 

  //if(N_SFH == 0)
  //    Delta_Eta = 1.0; 

  double reheated_energy = Delta_Eta * (stars * 1.0e10 / Hubble_h) * epsilon_energy * EnergySN * EtaSN; // Delta_Eta was in Msun so need to convert to stars to Msun. 

//  fprintf(stderr, "reheated_energy = %.4e \t Delta_Eta = %.4e \t stars = %.4e epsilon_energy = %.4e \t EnergySN = %.4e \t EtaSN = %.4e\n", reheated_energy, Delta_Eta /1.0e10 * Hubble_h, stars, epsilon_energy, EnergySN, EtaSN);
 
  return reheated_energy; 

}

// This function calculates the mass of stars that would have gone supernova in time t.
// Function and fit taken from Portinari et al. (1998).
//
// INPUT: Time (in Myr).
//
// OUTPUT: Mass of stars (in Msun) that would have gone Nova in time t.

double calculate_coreburning(double t)
{

  double a = 0.7473; // Fits from Portinari et al. (1998). 
  double b = -2.6979;
  double c = -4.7659;
  double d = 0.5934;

  double m = pow(10, a/log10(t) + b * exp(c/log10(t)) + d); 

  if (m < 8.0) // We enforce that SN-II only occur in stars M > 8Msun. 
    m = 8.0;

  return m; 
}

// If the cold gas that has been reheated has enough energy, it is possible to unbind some of the gas in the hot halo.
// This function determines this threshhold and how much mass is unbound and ejected.
//
// INPUT: Amount of cold gas reheated by supernova events (reheated_mass). ## UNITS: 1.0e10Msun/h (code units).
// 	: Amount of energy injected from the supernova events (reheated_energy). ## UNITS: erg.
// 	: Virial velocity of the background FoF halo (Vvir). ## UNITS: km/s.
//
// OUTPUT: Mass of ejected material. ## UNITS: 1.0e10Msun/h (code units).
//
// Note: Be aware here that a number of unit changes have occurred.
// This is necessary because of how Joules are defined. 
// Possible TODO: Derive a factor so the jumping around of units is unecessary. 

double calculate_ejected_mass(double *reheated_mass, double reheated_energy, double Vvir) 
{
  double ejected_mass;
  double Delta_Ehot;
  if(Vvir > 0.0)
  {
    Delta_Ehot = (0.5 * (*reheated_mass * 1.0e10 / Hubble_h * SOLAR_MASS / 1.0e3) * (Vvir * 1.0e3) * (Vvir * 1.0e3)) * 1.0e7; // Change in the thermal energy of the hot resevoir in erg..
    // Note the units here. Want final answer in Ergs (which is 1e-7 Joules) which has SI units of kg m^2 s^-2.
    // We change the mass from 1.0e10Msun/h to kg.  Change virial velocity from km/s to m/s. Finally change Joules to erg. 
  
//    printf("reheated_energy = %.4e \t Delta_Ehot = %.4e \t Difference in energy = %.4e \t Gal[centralgal].Vvir = %.4e\n", reheated_energy, Delta_Ehot, reheated_energy - Delta_Ehot, Vvir); 
    if(reheated_energy > Delta_Ehot) // Have enough energy to have some ejected mass.
    {
    	ejected_mass = (reheated_energy - Delta_Ehot) * 1.0e-7/(0.5 * Vvir * 1.0e3 * Vvir * 1.0e3); // Balance between the excess thermal energy and the thermal energy of the hot gas.
	// Energy changed from erg to Joules (kg m^2 s^-2) then divided by virial velocity (converted to m/s) giving final units of kg.
        ejected_mass = ejected_mass * 1.0e3 / SOLAR_MASS / 1.0e10 * Hubble_h; // Convert back from kg to 1.0e10*Msun/h. 
    }
    else // Not enough energy to eject mass. 
    {
	ejected_mass = 0.0;
//	printf("Reheated_mass before = %.4e ", *reheated_mass);
        *reheated_mass = reheated_energy * 1.0e-7 / (0.5 * Vvir * 1.0e3 * Vvir * 1.0e3); // Amount of mass that can be reheated to virial temperature of the hot halo.
	*reheated_mass = *reheated_mass * 1.0e3 / SOLAR_MASS / 1.0e10 * Hubble_h; 
//	printf("Reheated_mass after = %.4e\n", *reheated_mass);
    }
  } 					
  else
    ejected_mass = 0.0;
		
  if(ejected_mass < 0.0)
    ejected_mass = 0.0;

  //printf("reheated_energy = %.4e \t Delta_Ehot = %.4e \t Difference in energy = %.4e \t Gal[centralgal].Vvir = %.4e \t Ejected_mass = %.4e\n", reheated_energy, Delta_Ehot, reheated_energy - Delta_Ehot, Vvir, ejected_mass); 
  return ejected_mass;

}

void do_previous_recycling(int p, int centralgal, int SnapHistory)
{

  double MassWeightedStellarAge_Num = 0.0, MassWeightedStellarAge_Denom = 0.0, MassWeightedStellarAge; 
  double t_low, t_high;
  double m_low, m_high;
  double Delta_Eta, Delta_m;
  double mass_stars_recycled;
  
  int i;

  for(i = 0; i < Gal[p].SnapNum - SnapHistory; ++i)
  { 
    MassWeightedStellarAge_Denom += Gal[p].SNStars[i];
    MassWeightedStellarAge_Num += Gal[p].SNStars[i] * Age[i]; 

  }
 
  if(MassWeightedStellarAge_Denom != 0)
  {
    MassWeightedStellarAge = MassWeightedStellarAge_Num / MassWeightedStellarAge_Denom;

    t_high = (MassWeightedStellarAge - Age[Gal[p].SnapNum - 1]) * UnitTime_in_Megayears / Hubble_h;
    m_high = calculate_coreburning(t_high);

    t_low = (MassWeightedStellarAge - Age[Gal[p].SnapNum]) * UnitTime_in_Megayears / Hubble_h;
    m_low = calculate_coreburning(t_low);

    calculate_Delta_Eta(m_low, m_high, &Delta_Eta, &Delta_m); // Calculate the number and mass fraction of stars from snapshot i that go nova in the snapshot we are evolving FROM. 
    mass_stars_recycled = Delta_m * MassWeightedStellarAge_Denom; // Update the amount of stellar mass recycled from previous stars that have gone nova.

    //fprintf(stderr, "m_low = %.4e \t m_high = %.4e \t t_low = %.4e \t t_high = %.4e \t MassWeightedStellarAge = %.4e \t N_Stars = %.4e \t Delta_m = %.4e \t mass_stars_recycled = %.4e\n", m_low, m_high, t_low, t_high, MassWeightedStellarAge * UnitTime_in_Megayears / Hubble_h, MassWeightedStellarAge_Denom, Delta_m, mass_stars_recycled);

    update_from_SN_feedback(p, centralgal, 0.0, 0.0, mass_stars_recycled, 0.0);
  }
}

void do_previous_SN(int p, int centralgal, int halonr)
{
  double t;

  int i, SnapHistory = 0;
  double m_low, m_high, t_low, t_high;
  double Delta_Eta, Delta_m;
  
  double reheated_mass = 0.0; 
  double reheated_energy = 0.0;
  double mass_stars_recycled = 0.0;
  double mass_metals_new = 0.0;
  double ejected_mass;

  if (N_SFH == -1) // If we have set to do dynamic calculation for N_SFH.
  {
    while (t < 40)
    {
      t += (Age[Gal[p].SnapNum - N_SFH - 1] - Age[Gal[p].SnapNum - N_SFH]) * UnitTime_in_Megayears / Hubble_h;
      ++SnapHistory;

    }
    //--N_SFH_local; // If the time between snapshots is > 40Myr then N_SFH_local = 1 from the above loop but should be 0.  Hence need to subtract 1 off.
  }  else
  {
    SnapHistory = N_SFH; // Otherwise set the value to the user defined one.
  }

  if(SupernovaRecipeOn == 1) // Calculating the ejected mass due to previous star formation episodes.
  {  

    
    do_previous_recycling(p, centralgal, SnapHistory);

    for(i = Gal[p].SnapNum - SnapHistory; i < Gal[p].SnapNum - 1; ++i)
    {

      // Here we calculate the properties (energy/mass etc) for supernova using stars that formed in snapshot i and have finally gone nova in snapshot we are evolving FROM (Gal[p].SnapNum - 1). 
      // NOTE: Close attention should be given to the snapshot numbering here.
      // Gal[p].SnapNum is the snapshot we are evolving TO.  
      // Hence 'Gal[p].SnapNum - 1' can actually be thought of the snapshot where the current SF is occurring. I call this snapshot the Snapshot we are evolving FROM. 

      // First calculate the smallest star (which formed in snapshot i) which would have expended its H and He and gone nova. 
      // This defines the mass boundary below which a star will not go nova in the snapshot we are evolving from.
      t_low = ((Age[i] + Age[i + 1]) / 2.0 - Age[Gal[p].SnapNum]) * UnitTime_in_Megayears / Hubble_h; // This is the time between the SF episode of Snapshot i and the BEGINNING of the Snapshot we are evolving TO (Gal[p].SnapNum).
      m_low = calculate_coreburning(t_low);    

      // Next we calcualte the largest stellar mass (which formed in snapshot i) which would have expended its H and He and gone nova in the Snapshot we are evolving FROM.
      // This defines the mass boundary beyond which a star would have gone nova BEFORE the snapshot we are evolving from.
      t_high = ((Age[i] + Age[i + 1]) / 2.0 - Age[Gal[p].SnapNum - 1]) * UnitTime_in_Megayears / Hubble_h; // This is the time between between the SF episode of Snapshot i and the BEGINNING of the Snapshot we are evolving FROM (Gal[p].SnapNum - 1).
      m_high = calculate_coreburning(t_high);

      if (m_high < 8) // In this instance every star that has formed in snapshot i has already gone nova and accounted for.
          continue;

      calculate_Delta_Eta(m_low, m_high, &Delta_Eta, &Delta_m); // Calculate the number and mass fraction of stars from snapshot i that go nova in the snapshot we are evolving FROM. 
      reheated_mass += calculate_reheated_mass(Delta_Eta, Gal[p].SNStars[i], Gal[centralgal].Vmax); // Update the amount of mass reheated from previous stars that have gone nova.
//fprintf(stderr, "alpha_mass = %.4e \t Gal[centralgal].Vmax = %.4e \t V_mass = %.4e \t -beta_mass = %.4e\n", alpha_mass, Gal[centralgal].Vmax, V_mass, -beta_mass); 
//fprintf(stderr, "epsilon_mass = %.4e \t Gal[p].ColdGas = %.4e \t Vmax = %.4e \t reheated_mass_calculate = %.4e \t Delta_Eta = %.4e \t reheated_mass = %.4e \t ejected_mass = %.4e \t Gal[p].SNStars[i] = %.4e\n", epsilon_mass, Gal[p].ColdGas, Gal[centralgal].Vmax, Delta_Eta/Eta_SNII * Gal[p].SNStars[i] * epsilon_mass, Delta_Eta, reheated_mass, ejected_mass, Gal[p].SNStars[i]);

      reheated_energy += calculate_reheated_energy(Delta_Eta, Gal[p].SNStars[i], Gal[centralgal].Vmax); // Update the energy injected from previous stars that have gone nova. 
      mass_stars_recycled += Delta_m * Gal[p].SNStars[i]; // Update the amount of stellar mass recycled from previous stars that have gone nova.
      mass_metals_new += Delta_m / m_SNII * Yield * Gal[p].SNStars[i]; // Update the amount of new metals that the supernova has enriched the ISM with.
//      fprintf(stderr, "Delayed SN: Delta_ETa = %.4e \t stars = %.4e \t Gal[p].ColdGas = %.4e \t reheated_mass = %.4e \t reheated_energy = %.4e\n", Delta_Eta, Gal[p].SNStars[i], Gal[p].ColdGas, reheated_mass, reheated_energy);
    } 
  } 

  if(reheated_mass > Gal[p].ColdGas) // Can't reheated more cold gas than we currently have.
      reheated_mass = Gal[p].ColdGas;
    
  assert(reheated_mass >= 0.0); // Just make sure we're doing this right.
  assert(reheated_energy >= 0.0);
  assert(mass_stars_recycled >= 0.0);
  assert(mass_metals_new >= 0.0);

  //fprintf(stderr, "Calling ejected mass in delayed\n");
  ejected_mass = calculate_ejected_mass(&reheated_mass, reheated_energy, Gal[centralgal].Vvir); // Calculate the amount of mass ejected from supernova events. 
//  fprintf(stderr, "Delayed: Tot_Reheated_Energy = %.4e \t tot_reheated_mass = %.4e \t Ejected_mass = %.4e\n", reheated_energy, reheated_mass, ejected_mass);
  update_from_SN_feedback(p, centralgal, reheated_mass, ejected_mass, mass_stars_recycled, mass_metals_new); 


 
}

void do_current_SN(int p, int centralgal, int halonr, double *stars, double *reheated_mass, double *mass_metals_new, double *mass_stars_recycled, double *ejected_mass)
{

  *reheated_mass = 0.0; 
  double reheated_energy = 0.0;
  *mass_stars_recycled = 0.0;
  *mass_metals_new = 0.0;
  *ejected_mass = 0.0;

  double m_low = 8.0;
  double m_high = 120.0;
  double t_low;

  double Delta_Eta = Eta_SNII;
  double Delta_m = RecycleFraction;
  double mass_fraction = 1.0; 

  if (N_SFH != 0) // If we have set N_SFH to 0 then we wish to use the IRA then need to calculate the correct recycle fraction. 
  // If we are using the IRA then the declarations we did above will be correct. 
  {
    
    // First calculate the smallest star (which formed in snapshot i) which would have expended its H and He and gone nova. 
    // This defines the mass boundary below which a star will not go nova in the snapshot we are evolving from.
    t_low = ((Age[Gal[p].SnapNum - 1] - Age[Gal[p].SnapNum]) / 2.0) * UnitTime_in_Megayears / Hubble_h; // This is the time between the SF episode of this snapshot and the snapshot we are evolving TO. 
    m_low = calculate_coreburning(t_low);    

    calculate_Delta_Eta(m_low, m_high, &Delta_Eta, &Delta_m); // Calculate the number and mass fraction of stars from snapshot i that go nova in the snapshot we are evolving FROM.      
    mass_fraction = Delta_m / m_SNII; // Mass fraction of stars that went nova. 
  } 
  
  *mass_stars_recycled = Delta_m * (*stars); 
  *reheated_mass = calculate_reheated_mass(Delta_Eta, *stars, Gal[centralgal].Vmax); // Update the amount of mass reheated from stars that have gone nova. 

  if(*reheated_mass > Gal[p].ColdGas) // Can only reheat as much cold gas we have available.
      *reheated_mass = Gal[p].ColdGas;
  if((*reheated_mass + *stars > Gal[p].ColdGas) && (*reheated_mass + *stars > 0.0))
  {
    double factor = Gal[p].ColdGas / (*reheated_mass + *stars);
    *reheated_mass *= factor;
    *stars *= factor; 
  }


  *mass_metals_new = mass_fraction * Yield * (*stars); // Update the amount of new metals that the supernova has enriched the ISM with.
//  fprintf(stderr, "Calling ejected in current\n");
  reheated_energy = calculate_reheated_energy(Delta_Eta, *stars, Gal[centralgal].Vmax); // Update the energy injected from previous stars that have gone nova. 
  
  assert(*reheated_mass >= 0.0); // Just make sure we're doing this right.
  assert(reheated_energy >= 0.0);
  assert(*mass_stars_recycled >= 0.0);
  assert(*mass_metals_new >= 0.0);

  *ejected_mass = calculate_ejected_mass(&(*reheated_mass), reheated_energy, Gal[centralgal].Vvir); // Calculate the amount of mass ejected from supernova events. 
//    fprintf(stderr, "Current SN: stars = %.4e \t Gal[p].ColdGas = %.4e \t reheated_mass = %.4e \t ejected_mass = %.4e \t reheated_energy = %.4e\n", *stars, Gal[p].ColdGas, *reheated_mass, *ejected_mass, reheated_energy);
}

