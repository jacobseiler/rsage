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

#define IMF_CONSTANT IMF_norm/(IMF_slope + 1.0) 

// Local Proto-Types //

void  update_SN_stars_array(int p, double stars, double dt, int tree, int ngal);
void update_stellar_tracking(int p, double stars, double dt, int tree, int ngal);

// External Functions //

void update_from_SN_feedback(int p, int centralgal, double reheated_mass, double ejected_mass, double mass_stars_recycled, double mass_metals_new, double NSN, double dt) 
{

  double metallicity, metallicityHot, dust_fraction_cold, dust_fraction_hot;
  double dust_fraction_coldhot, dust_fraction_hotcold; // These are the dust fractions defined by M_colddust / (Mcold + Mhot) and M_hotdust / (Mcold + Mhot). 
  // These second dust fractions are required because when we annihilate dust from SN we want the TOTAL amount to be NSN * eta * Ms.

  // Dust constants, taken mainly from Dayal, Ferrara and Saro 2011.
  double yd = 0.5 / 1.0e10 * Hubble_h; // Amount of dust mass created from each supernova event. Units of Msun.
  double eta = 0.4; // Coupling between shock energy and annihilation of dust.
  double Ms = 6.8e3 / 1.0e10 * Hubble_h; // Amount of mass that is accelerated to 100km/s by SN blastwave. Units of Msun. 
  double dust_added = 0.0, dust_removed = 0.0;

  XASSERT((reheated_mass >= 0.0) && (ejected_mass >= 0.0) && (mass_stars_recycled >= 0.0) && (mass_metals_new >= 0.0), "When trying to update from the supernova feedback we got negative masses.\nReheated Mass = %.4e (1.0e10 Msun/h) \t Ejected Mass = %.4e (1.0e10 Msun/h) \t Mass of Stars Recycled = %.4e (1.0e10 Msun/h) \t Mass of Metals Added = %.4e (1.0e10 Msun/h)\n", reheated_mass, ejected_mass, mass_stars_recycled, mass_metals_new);
 
  Gal[p].StellarMass -= mass_stars_recycled; // The supernova remnants are instantly moved back into the cold ISM. 
  Gal[p].ColdGas += mass_stars_recycled; 

  Gal[p].DustColdGas += yd * NSN; // Dust is created by supernova (something) on cold gas. 
  Gal[p].ColdGas -= yd * NSN; // Hence the created dust is subtracted from the cold gas reservoir.

  dust_added += yd * NSN;

  if(Gal[p].ColdGas > 1.0e-10)
      Gal[p].MetalsColdGas += mass_metals_new; // ISM is enriched by new metals.
  else
      Gal[centralgal].MetalsHotGas += mass_metals_new;

  // Here we do checks to see if we need to rescale the amount of mass that has been ejected/reheated. //

  if(reheated_mass > Gal[p].ColdGas) // Just perform this check again to be sure.
      reheated_mass = Gal[p].ColdGas; // Because SF takes gas out of cold resevoir.
 
  if(ejected_mass > Gal[centralgal].HotGas)
      ejected_mass = Gal[centralgal].HotGas; 

  metallicity = get_metallicity(Gal[p].ColdGas, Gal[p].MetalsColdGas);
  dust_fraction_cold = get_dust_fraction(Gal[p].ColdGas, Gal[p].DustColdGas);
  dust_fraction_coldhot = get_dust_fraction(Gal[p].ColdGas + Gal[p].HotGas + Gal[p].EjectedMass, Gal[p].DustColdGas);

  Gal[p].ColdGas -= reheated_mass;
  Gal[p].MetalsColdGas -= reheated_mass * metallicity;
  Gal[p].DustColdGas -= reheated_mass * dust_fraction_cold; // Dust that has been reheated to the hot reservoir.
  Gal[p].DustColdGas -= NSN * eta * Ms * dust_fraction_coldhot; // Dust that has been annihilated by shocks.

  dust_removed += NSN * eta * Ms * dust_fraction_coldhot;

  Gal[centralgal].HotGas += reheated_mass;
  Gal[centralgal].MetalsHotGas += reheated_mass * metallicity;
  Gal[centralgal].DustHotGas += reheated_mass * dust_fraction_cold; 

  metallicityHot = get_metallicity(Gal[centralgal].HotGas, Gal[centralgal].MetalsHotGas);
  dust_fraction_hot = get_dust_fraction(Gal[p].HotGas, Gal[p].DustHotGas);
  dust_fraction_hotcold = get_dust_fraction(Gal[p].ColdGas + Gal[p].HotGas + Gal[p].EjectedMass, Gal[p].DustHotGas);
  // When dust_fraction_hotcold is calculated is irrelevant (i.e., before or after mass is added/subtracted from ColdGas).
  // This is because mass is conserved between ColdGas and HotGas during this step.

  Gal[centralgal].HotGas -= ejected_mass;
  Gal[centralgal].MetalsHotGas -= ejected_mass * metallicityHot;
  Gal[centralgal].DustHotGas -= ejected_mass * dust_fraction_hot; // Dust that has been ejected from the galaxy.
  Gal[centralgal].DustHotGas -= NSN * eta * Ms * dust_fraction_hotcold; // Dust that has been annihilated by shocks.
  dust_removed += NSN * eta * Ms * dust_fraction_hotcold;
 
  Gal[centralgal].EjectedMass += ejected_mass;
  Gal[centralgal].MetalsEjectedMass += ejected_mass * metallicityHot;

  Gal[centralgal].DustEjectedMass += ejected_mass * dust_fraction_hot; 

  Gal[p].GridOutflowRate[Halo[Gal[centralgal].HaloNr].SnapNum] += ejected_mass / dt; 

  Gal[p].GrandSum -= mass_stars_recycled;

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
  if(Gal[p].GrandSum < 0.0)
    Gal[p].GrandSum = 0.0;

  if (Gal[p].DustColdGas < 0.0)
    Gal[p].DustColdGas = 0.0;
  if (Gal[p].DustHotGas < 0.0)
    Gal[p].DustHotGas = 0.0;
 
}


void update_from_star_formation(int p, double stars, double dt, int step, bool ismerger, int tree, int ngal)
{
 
  double metallicity = get_metallicity(Gal[p].ColdGas, Gal[p].MetalsColdGas);
  double dust_fraction_cold = get_dust_fraction(Gal[p].ColdGas, Gal[p].DustColdGas);
  double time_spanned;

  // If the SF episode was from a merger the stars are placed in different reservoirs. 
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
    Gal[p].MetalsBulgeMass += stars;

  }

  // If we're doing delayed SN we need to update the tracking 
  if (ismerger == true)
  { 
    time_spanned = 0.0;
  }
  else
  {
    time_spanned = dt;
  }
  
  if (IRA == 0)
  { 
    update_SN_stars_array(p, stars, time_spanned, tree, ngal);
  }

  if (PhotonPrescription == 1)
  {
    //update_stellar_tracking(p, stars, time_spanned, tree, ngal);
  }

  // update gas and metals from star formation 
  Gal[p].ColdGas -= stars;
  Gal[p].MetalsColdGas -= metallicity * stars;
  Gal[p].StellarMass += stars;
  Gal[p].MetalsStellarMass += metallicity * stars;
  Gal[p].DustColdGas -= dust_fraction_cold * stars; // When stars form they eat up dust. 
  if (Gal[p].ColdGas < 0.0)
    Gal[p].ColdGas = 0.0;
  if (Gal[p].MetalsColdGas < 0.0)
    Gal[p].MetalsColdGas = 0.0;

}

void starformation_and_feedback(int p, int centralgal, double time, double dt, int halonr, int step, int tree, int ngal)
{
  double reff, tdyn, strdot, stars;
  double cold_crit; 
  double reheated_mass = 0.0, mass_metals_new = 0.0, mass_stars_recycled = 0.0, ejected_mass = 0.0, NSN = 0.0;

  // Initialise variables
  strdot = 0.0;

  // First we check to see if we need to do delayed supernova feedback.

  if(IRA == 0)
  {
    do_previous_recycling(p, centralgal, step, dt); 

    // The total amount that the reservoirs need to be adjusted by for this evolution step (i.e., ALL steps) has already been determined.
    // Within each substep (i.e., each `step` iteration) we ask 'how many iterations are left?' and then adjusts the reservoirs by this amount.
    // After the reservoirs have been updated, we subtract the amount of mass that we added from the running totals.
    //
    // These running totals are required because contemporaneous SN could cause extra mass to be added/removed within the next step iteration. 

    update_from_SN_feedback(p, centralgal, Gal[p].reheated_mass / (STEPS-step), Gal[p].ejected_mass / (STEPS-step), Gal[p].mass_stars_recycled / (STEPS-step), Gal[p].mass_metals_new / (STEPS-step), Gal[p].NSN / (STEPS-step), dt);

    Gal[p].reheated_mass -= Gal[p].reheated_mass / (STEPS-step);
    Gal[p].ejected_mass -= Gal[p].ejected_mass / (STEPS-step);
    Gal[p].mass_stars_recycled -= Gal[p].mass_stars_recycled / (STEPS-step);
    Gal[p].mass_metals_new -= Gal[p].mass_metals_new / (STEPS-step);
    Gal[p].NSN -= Gal[p].NSN / (STEPS-step);

  }

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
    if(IRA == 0)
      do_contemporaneous_SN(p, centralgal, dt, &stars, &reheated_mass, &mass_metals_new, &mass_stars_recycled, &ejected_mass, &NSN);  
    else if(IRA == 1)
    {
      do_IRA_SN(p, centralgal, &stars, &reheated_mass, &mass_metals_new, &mass_stars_recycled, &ejected_mass, &NSN); 
    } 

  } 
  
  if(stars > Gal[p].ColdGas) // we do this check in 'do_current_sn()' but if supernovarecipeon == 0 then we don't do the check.
  {
    double factor = Gal[p].ColdGas / stars; 
    stars *= factor; 
  }
  
  update_from_star_formation(p, stars, dt, step, false, tree, ngal);
 
  update_from_SN_feedback(p, centralgal, reheated_mass, ejected_mass, mass_stars_recycled, mass_metals_new, NSN, dt);

  // check for disk instability
  if(DiskInstabilityOn)
    check_disk_instability(p, centralgal, halonr, time, dt, step, tree, ngal);

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

  int32_t bin_idx_low = round((m_low - m_IMFbins_low) / (m_IMFbins_delta)); 
  int32_t bin_idx_high = round((m_high - m_IMFbins_low) / (m_IMFbins_delta)); 

  if (bin_idx_high == N_massbins)
    --bin_idx_high;

  if (bin_idx_low == N_massbins)
    --bin_idx_low;

  *Delta_Eta = IMF_CONSTANT * (IMF_massgrid_eta[bin_idx_high] - IMF_massgrid_eta[bin_idx_low]);
  *Delta_m = IMF_CONSTANT * (IMF_massgrid_m[bin_idx_high] - IMF_massgrid_m[bin_idx_low]);

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

  double epsilon_mass;

  if(RescaleSN == 0)
  {  
    epsilon_mass = FeedbackReheatingEpsilon;
  } else
  {
    epsilon_mass = alpha_mass * (0.5 + pow(Vmax/V_mass, -beta_mass));
  }



  if (epsilon_mass > epsilon_mass_max)
  {
    epsilon_mass = epsilon_mass_max; // We enforce a maximum value for the mass loading factor. 
  }

  double reheated_mass = Delta_Eta/Eta_SNII * stars * epsilon_mass;

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

  double epsilon_energy;

  if(RescaleSN == 0)
  {
    epsilon_energy = FeedbackEjectionEfficiency;
  } else  
  {
    epsilon_energy = alpha_energy * (0.5 + pow(Vmax/V_energy, -beta_energy));
  }

  double reheated_energy = 0.5 * Delta_Eta * (stars * 1.0e10 / Hubble_h) * epsilon_energy * EnergySN; // Delta_Eta was in Msun so need to convert to stars to Msun. 
 
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

  /*
  double a = 0.7473; // Fits from Portinari et al. (1998). 
  double b = -2.6979;
  double c = -4.7659;
  double d = 0.5934;

  double m = pow(10, a/log10(t) + b * exp(c/log10(t)) + d); 

  return m; 
  */
  int32_t bin_idx = (t - coreburning_tbins_low) / (coreburning_tbins_delta); 
  if (bin_idx < 0) // Time is so short that only stars with mass greater than 120Msun can go nova.
    return 120.0;

  if (bin_idx > N_tbins) // Time is so long that all stars within the IMF range can go nova.
    return 8.0; 

  return coreburning_times[bin_idx];

 
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

    if(reheated_energy > Delta_Ehot) // Have enough energy to have some ejected mass.
    {
    	ejected_mass = (reheated_energy - Delta_Ehot) * 1.0e-7/(0.5 * Vvir * 1.0e3 * Vvir * 1.0e3); // Balance between the excess thermal energy and the thermal energy of the hot gas.
	// Energy changed from erg to Joules (kg m^2 s^-2) then divided by virial velocity (converted to m/s) giving final units of kg.

        ejected_mass = ejected_mass * 1.0e3 / SOLAR_MASS / 1.0e10 * Hubble_h; // Convert back from kg to 1.0e10*Msun/h. 
    }
    else // Not enough energy to eject mass. 
    {
	ejected_mass = 0.0;
        *reheated_mass = reheated_energy * 1.0e-7 / (0.5 * Vvir * 1.0e3 * Vvir * 1.0e3); // Amount of mass that can be reheated to virial temperature of the hot halo.
	*reheated_mass = *reheated_mass * 1.0e3 / SOLAR_MASS / 1.0e10 * Hubble_h; 
    }

  } 					
  else
    ejected_mass = 0.0;
		
  if(ejected_mass < 0.0)
    ejected_mass = 0.0;

  return ejected_mass;

}

// This function answers the question, "How many stars formed in the previous time steps will explode during the current supernova timestep."
// Note: Since the time scale on which we calculate SN feedback is variable (defined by the user in 'TimeResolutionSN') we use the terminology 'supernova timestep'.
// Double note: This function is different to "do_contemporaneous_SN" because that looks at the stars formed during the CURRENT STAR FORMATION time step, not previous ones.
void do_previous_SN(int p, int centralgal, double dt)
{
  double reheated_mass = 0.0; 
  double reheated_energy = 0.0;
  double mass_stars_recycled = 0.0;
  double mass_metals_new = 0.0;
  double ejected_mass = 0.0;
  double NSN = 0.0; // This is the number of supernova events during this time step.

  int i; 
  double m_low = -1.0, m_high = -1.0, t_low = -1.0, t_high = -1.0;
  double Delta_Eta, Delta_m;
  double time_until_next_SN;

  if(dt * UnitTime_in_Megayears / Hubble_h < TimeResolutionSN) // If the star formation time scale is smaller than the time scale on which we do SN feedback
    time_until_next_SN = TimeResolutionSN; // Then the time that we next calculate delayed SN feedback would be given dictated by the time resolution of SN. 
  else // Otherwise the star formation time scale is larger than the SN feedback time scale.
  {
    time_until_next_SN = dt * UnitTime_in_Megayears / Hubble_h;  // Then the time that we next calculate delayed SN feedback is in the next SF timestep.
  } 

  if(SupernovaRecipeOn == 1) // Calculating the ejected mass due to previous star formation episodes.
  {   
    for(i = 1; i < SN_Array_Len; ++i)
    {
      if(Gal[p].SN_Stars[i] < 1e-10)
	    continue;

      // First calculate the smallest star which would have expended its H and He and gone supernova. 
      // This defines the mass boundary below which a star will not go nova in the current SN step.

      t_low = (i * TimeResolutionSN) / 2.0 + time_until_next_SN; // (i * TimeResolutionSN) is the number of stars that were formed that many megayears ago.  time_until_next_SN allows us to ask how many of stars will explode in this SN step.
      if(t_low < 2)
  	    t_low = 2;
 
      m_low = calculate_coreburning(t_low);    

      if (m_low < 8.0) // We enforce that SN-II only occur in stars M > 8Msun. 
        m_low = 8.0;

      if (m_low > 120.0)
        continue;

      // Next we calculate the largest stellar mass which would have expended its H and He and gone supernova.
      // This defines the mass boundary beyond which a star would have gone nova BEFORE the current SN step. 

      t_high = (i * TimeResolutionSN) / 2.0; // (i * TimeResolutionSN) is the number of stars that were formed that many Megayears ago.  The / 2.0 factor arises from the fact that we assume stars were formed in a single co-eval burst in the middle of this period.
      if(t_high < 2)
        t_high = 2;

      m_high = calculate_coreburning(t_high);


      if (m_high < 8) // In this instance every star that has formed in i Myr ago has already gone supernova and accounted for.
          continue;

      if (m_high > 120.0)
          continue;

      calculate_Delta_Eta(m_low, m_high, &Delta_Eta, &Delta_m); // Calculate the number and mass fraction of stars i Myr ago that go supernova in the current SF step. 
      reheated_mass += calculate_reheated_mass(Delta_Eta, Gal[p].SN_Stars[i], Gal[centralgal].Vmax); // Update the amount of mass reheated from previous stars that have gone nova.

      reheated_energy += calculate_reheated_energy(Delta_Eta, Gal[p].SN_Stars[i], Gal[centralgal].Vmax); // Update the energy injected from previous stars that have gone nova. 
      mass_stars_recycled += Delta_m * Gal[p].SN_Stars[i]; // Update the amount of stellar mass recycled from previous stars that have gone nova.
      mass_metals_new += Delta_m / m_SNII * Yield * Gal[p].SN_Stars[i]; // Update the amount of new metals that the supernova has enriched the ISM with.

      NSN += Delta_Eta * Gal[p].SN_Stars[i] * 1.0e10 / Hubble_h; // The number of supernova events will be simply given by the number fraction of stars that exploded times the mass of stars. 
 
      XASSERT(reheated_mass >= 0.0, "i = %d \t Reheated mass = %.4e \t t_low = %.4e Myr \t m_low = %.4e \t t_high = %.4e Myr \t m_high = %.4e\n", i, reheated_mass, t_low, m_low, t_high, m_high); // Just make sure we're doing this right.

    } 
    
    if(reheated_mass > Gal[p].ColdGas) // Can't reheated more cold gas than we currently have.
        reheated_mass = Gal[p].ColdGas;
    
    XASSERT(reheated_mass >= 0.0, "Reheated mass = %.4e \t t_low = %.4e Myr \t m_low = %.4e \t t_high = %.4e Myr \t m_high = %.4e\n", reheated_mass, t_low, m_low, t_high, m_high); // Just make sure we're doing this right.
    assert(reheated_energy >= 0.0);
    assert(mass_stars_recycled >= 0.0);
    assert(mass_metals_new >= 0.0);

    ejected_mass = calculate_ejected_mass(&reheated_mass, reheated_energy, Gal[centralgal].Vmax); // Calculate the amount of mass ejected from supernova events.

    XASSERT(ejected_mass >= 0.0, "For galaxy %d the ejected mass was %.4e \t reheated_mass = %.4e\n", p, ejected_mass, reheated_mass); 
  }

  Gal[p].reheated_mass = reheated_mass;
  Gal[p].ejected_mass = ejected_mass;
  Gal[p].mass_stars_recycled = mass_stars_recycled;
  Gal[p].mass_metals_new = mass_metals_new;
  Gal[p].NSN = NSN;
 
}

/*
void do_current_SN(int p, int centralgal, int halonr, double *stars, double *reheated_mass, double *mass_metals_new, double *mass_stars_recycled, double *ejected_mass)
{

  *reheated_mass = 0.0; 
  double reheated_energy = 0.0;
  *mass_stars_recycled = 0.0;
  *mass_metals_new = 0.0;
  *ejected_mass = 0.0;

  double Delta_Eta = Eta_SNII;
  double Delta_m = RecycleFraction;
  double mass_fraction = 1.0; 

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
  reheated_energy = calculate_reheated_energy(Delta_Eta, *stars, Gal[centralgal].Vmax); // Update the energy injected from previous stars that have gone nova. 
  
  assert(*reheated_mass >= 0.0); // Just make sure we're doing this right.
  assert(reheated_energy >= 0.0);
  assert(*mass_stars_recycled >= 0.0);
  assert(*mass_metals_new >= 0.0);

  *ejected_mass = calculate_ejected_mass(&(*reheated_mass), reheated_energy, Gal[centralgal].Vvir); // Calculate the amount of mass ejected from supernova events. 

}
*/

void do_previous_recycling(int p, int centralgal, int step, double dt) 
{
 
  if(Gal[p].StellarAge_Numerator > 0.0)
  {  
    double t_low, t_high;
    double m_low, m_high;
    double Delta_Eta, Delta_m;
    double mass_stars_recycled;
    double NSN; 
 
    double mwmsa = Gal[p].StellarAge_Numerator / Gal[p].StellarAge_Denominator; // The numerator is weighted by the age of the stars with the denominator simply being the number of stars.
    double time_into_snap = dt * step; 


    t_high = (mwmsa - ((Age[Gal[p].SnapNum - 1]) - time_into_snap)) * UnitTime_in_Megayears / Hubble_h;
    if (t_high < 2.0)
      t_high = 2.0;
    m_high = calculate_coreburning(t_high);

    t_low = (mwmsa - (Age[Gal[p].SnapNum] - time_into_snap)) * UnitTime_in_Megayears / Hubble_h;
    if (t_low < 2.0)
      t_low = 2.0;
    m_low = calculate_coreburning(t_low);

    calculate_Delta_Eta(m_low, m_high, &Delta_Eta, &Delta_m); // Calculate the number and mass fraction of stars from snapshot i that go nova in the snapshot we are evolving FROM. 
    mass_stars_recycled = Delta_m * Gal[p].StellarAge_Denominator; // Update the amount of stellar mass recycled from previous stars that have gone nova.
    NSN = Delta_Eta * Gal[p].StellarAge_Denominator * 1.0e10 / Hubble_h;

    update_from_SN_feedback(p, centralgal, 0.0, 0.0, mass_stars_recycled, 0.0, NSN, dt);
  } 
}

// In this function we answer the question, "How many stars formed in the current star formation time step will explode by the time we next calculate our supernova feedback?"
// Note: This function only concerns itself with stars formed in the CURRENT STAR FORMATION time step, unlike "do_previous_SN" which focuses on calculated SN feedback from PREVIOUS STAR FORMATION time steps.

void do_contemporaneous_SN(int p, int centralgal, double dt, double *stars, double *reheated_mass, double *mass_metals_new, double *mass_stars_recycled, double *ejected_mass, double *NSN)
{

  *reheated_mass = 0.0; 
  double reheated_energy = 0.0;
  *mass_stars_recycled = 0.0;
  *mass_metals_new = 0.0;
  *ejected_mass = 0.0;
 
  double m_low, m_high, t_low;
  double Delta_Eta, Delta_m;
  double fac;
 
  if(dt * UnitTime_in_Megayears / Hubble_h / 2.0 < TimeResolutionSN) // If the star formation time scale is smaller than the time scale on which we do SN feedback
    t_low = (TimeResolutionSN - Gal[p].Total_SN_SF_Time); // Then our 'sub-grid' SN feedback time will be the time from this SF episode until we next calculate SN feedback (divided by 2 as we assume the stars are formed in the middle of the interval).
  else // Otherwise the star formation time scale is larger than the SN feedback time scale.
    t_low  = (dt * UnitTime_in_Megayears / Hubble_h) / 2.0; // Then the feedback time will be the time from this SF event to the next star formation event. This is because SN feedback is only calculated when star formation occurs regardless of the actual value of 'TimeResolutionSN'. 
   
  if (t_low < 2.0)
    t_low = 2.0;
 
  m_low = calculate_coreburning(t_low);    
  if(m_low < 8.0)
    m_low = 8.0;
 
  if(m_low > 120.0)
    return;      
 
  m_high = 120.0; 

  calculate_Delta_Eta(m_low, m_high, &Delta_Eta, &Delta_m); // Calculate the number and mass fraction of stars i Myr ago that go supernova in the current SF step. 
  *reheated_mass = calculate_reheated_mass(Delta_Eta, *stars, Gal[centralgal].Vmax); // Update the amount of mass reheated from previous stars that have gone nova.

  if((*stars + *reheated_mass) > Gal[p].ColdGas && (*stars + *reheated_mass) > 0.0)
  {
    fac = Gal[p].ColdGas / (*stars + *reheated_mass);
    *stars *= fac;
    *reheated_mass *= fac;
  }   

  reheated_energy += calculate_reheated_energy(Delta_Eta, *stars, Gal[centralgal].Vmax); // Update the energy injected from previous stars that have gone nova. 
  *mass_stars_recycled += Delta_m * (*stars); // Update the amount of stellar mass recycled from previous stars that have gone nova.
  *mass_metals_new += Delta_m / m_SNII * Yield * (*stars); // Update the amount of new metals that the supernova has enriched the ISM with.

  if(*reheated_mass > Gal[p].ColdGas) // Can't reheated more cold gas than we currently have.
    *reheated_mass = Gal[p].ColdGas;

  *ejected_mass = calculate_ejected_mass(&(*reheated_mass), reheated_energy, Gal[centralgal].Vmax); // Calculate the amount of mass ejected from supernova events.    
  *NSN = Delta_Eta * (*stars * 1.0e10 / Hubble_h);

  assert(*reheated_mass >= 0.0); // Just make sure we're doing this right.
  assert(reheated_energy >= 0.0);
  assert(*mass_stars_recycled >= 0.0);
  assert(*mass_metals_new >= 0.0);
  assert(*ejected_mass >= 0.0);  

}

void do_IRA_SN(int p, int centralgal, double *stars, double *reheated_mass, double *mass_metals_new, double *mass_stars_recycled, double *ejected_mass, double *NSN)
{

  double fac;
  *reheated_mass = 0.0;
  *mass_metals_new = 0.0;
  *mass_stars_recycled = 0.0;
  *ejected_mass = 0.0;

  *reheated_mass = calculate_reheated_mass(Eta_SNII, *stars, Gal[centralgal].Vmax);
  *NSN = Eta_SNII * (*stars);

  if((*stars + *reheated_mass) > Gal[p].ColdGas && (*stars + *reheated_mass) > 0.0)
  {
    fac = Gal[p].ColdGas / (*stars + *reheated_mass);
    *stars *= fac;
    *reheated_mass *= fac;
  }   

  double reheated_energy = calculate_reheated_energy(Eta_SNII, *stars, Gal[centralgal].Vmax); 

  *ejected_mass = calculate_ejected_mass(&(*reheated_mass), reheated_energy, Gal[centralgal].Vmax);
	
  if(*ejected_mass < 0.0)
    *ejected_mass = 0.0;
       
  *mass_metals_new = Yield * (*stars);    
  *mass_stars_recycled = RecycleFraction * (*stars);

}

// Local Functions //

// This function updates the array which holds the amount of stars formed in the past 50 Myr.
// As we only need this array for doing delayed SN, if we are using the IRA then this function will never be called.
// If the SF timestep is less than the resolution on which we do supernova feedback then we will keep track of the total stars formed over these small timesteps and then update them all in one bin.
// If the SF timestep is larger than the resolution on which we do supernova feedback we will spread the stars formed evenly over a number of elements (>= 1). 
//
// INPUT: The index of the galaxy (p).
// 	: The number of stars formed in the current SF episode (stars). ## UNITS: 1.0e10 Msun/h (Code Units). 
// 	: The timestep over which the SF is occuring (dt). ## UNITS: Code Units, multiply by 'UnitTime_in_Megayears / Hubble_h' for Myr.
// 	: The tree currently being used (tree); currently used for debugging purposes.
//
// OUTPUT: None.

void update_SN_stars_array(int p, double stars, double dt, int tree, int ngal)
{

  double time_spanned = dt * UnitTime_in_Megayears / Hubble_h; // The time spanned by this star formation event.

  Gal[p].Total_SN_SF_Time += time_spanned; // How long it has been since we've updated the array?
  Gal[p].Total_SN_Stars += stars; // How many stars we will need to bin once we do update the array?

  if(Gal[p].Total_SN_SF_Time > TimeResolutionSN * SN_Array_Len) // This handles cases in which the time spanned is greater than 50Myr.  In this case we wipe the array clean and push the star formation into a 50Myr bin. 
    Gal[p].Total_SN_SF_Time = TimeResolutionSN * SN_Array_Len;

  if(Gal[p].Total_SN_SF_Time < TimeResolutionSN) // If it hasn't been long enough yet, don't update the array. 
    return;

  int num_shuffled = round(Gal[p].Total_SN_SF_Time/TimeResolutionSN); // How many cells will this SF event span.

  double stars_spread = Gal[p].Total_SN_Stars/num_shuffled; // We spread the stars evenly over time.

  int i;

  XASSERT(Gal[p].IsMalloced == 1, "We are attempting to update the stars array but this galaxy has already had its arrays freed.\nGalaxy %d \t Halo %d \t Tree %d \t Time spanned %.4eMyr\n", p, Gal[p].HaloNr, tree, time_spanned);  

  for(i = SN_Array_Len - 1; i > num_shuffled - 1; --i)
  {
      Gal[p].SN_Stars[i] = Gal[p].SN_Stars[i-num_shuffled]; // Shuffle the current elements of the array far enough along so we can store the new stars.
  }
    XPRINT(p < ngal, "We have the case where the galaxy p is greater than the number of galaxies.  p = %d \t ngal = %d\n", p, ngal);
  for(i = SN_Array_Len - 1; i > (SN_Array_Len - num_shuffled - 1); --i)
  {
    Gal[p].StellarAge_Numerator += Gal[p].SN_Stars[i] * (Age[Gal[p].SnapNum - 4]);
    Gal[p].StellarAge_Denominator += Gal[p].SN_Stars[i]; 
  }

  for(i = 0; i < num_shuffled; ++i) 
  {
    Gal[p].StellarAge_Numerator += Gal[p].SN_Stars[i] * (Age[Gal[p].SnapNum - 4]);
    Gal[p].SN_Stars[i] = stars_spread; // Update the vacated elements with the new stars.
    Gal[p].GrandSum += stars_spread; 
  } 

  Gal[p].Total_SN_SF_Time = 0.0; // We've updated so reset our variables.
  Gal[p].Total_SN_Stars = 0.0;
}


void update_stellar_tracking(int p, double stars, double dt, int tree, int ngal)
{

  double time_spanned = dt * UnitTime_in_Megayears / Hubble_h; // The time spanned by this star formation event.

  Gal[p].Total_Stellar_SF_Time += time_spanned; // How long it has been since we've updated the array?
  Gal[p].Total_Stellar_Stars += stars; // How many stars we will need to bin once we do update the array?

  if(Gal[p].Total_Stellar_SF_Time > TimeResolutionStellar * StellarTracking_Len) 
    Gal[p].Total_Stellar_SF_Time = TimeResolutionStellar * StellarTracking_Len;

  if(Gal[p].Total_Stellar_SF_Time < TimeResolutionStellar) // If it hasn't been long enough yet, don't update the array. 
    return;

  int num_shuffled = round(Gal[p].Total_Stellar_SF_Time/TimeResolutionStellar); // How many cells will this SF event span.

  double stars_spread = Gal[p].Total_Stellar_Stars/num_shuffled; // We spread the stars evenly over time.

  int i;

  XASSERT(Gal[p].IsMalloced == 1, "We are attempting to update the Stellar Tracking array but this galaxy has already had its arrays freed.\nGalaxy %d \t Halo %d \t Tree %d \t Time spanned %.4eMyr\n", p, Gal[p].HaloNr, tree, time_spanned);  

  for(i = StellarTracking_Len - 1; i > num_shuffled - 1; --i)
  {
      Gal[p].Stellar_Stars[i] = Gal[p].Stellar_Stars[i-num_shuffled]; // Shuffle the current elements of the array far enough along so we can store the new stars.
  }
  XPRINT(p < ngal, "We have the case where the galaxy p is greater than the number of galaxies.  p = %d \t ngal = %d\n", p, ngal);

  for(i = 0; i < num_shuffled; ++i) 
  {
    Gal[p].Stellar_Stars[i] = stars_spread; // Update the vacated elements with the new stars.  
  } 

  Gal[p].Total_Stellar_SF_Time = 0.0; // We've updated so reset our variables.
  Gal[p].Total_Stellar_Stars = 0.0;
}

