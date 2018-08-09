#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#include "core_allvars.h"
#include "core_proto.h"



void reincorporate_gas(int centralgal, double dt)
{
  double reincorporated, metallicity, dust_fraction_ejected;

  // SN velocity is 630km/s, and the condition for reincorporation is that the 
  // halo has an escape velocity greater than this, i.e. V_SN/sqrt(2) = 445.48km/s
  double Vcrit = 445.48 * ReIncorporationFactor;  
  double channel_fractionSN, channel_fractionQSO;
  channel_fractionSN = channel_fractionQSO = 0.0;
	
  if(Gal[centralgal].Vvir > Vcrit)
  {
    reincorporated = 
      ( Gal[centralgal].Vvir / Vcrit - 1.0 ) *
				Gal[centralgal].EjectedMass / (Gal[centralgal].Rvir / Gal[centralgal].Vvir) * dt; 

    if(reincorporated > Gal[centralgal].EjectedMass)
      reincorporated = Gal[centralgal].EjectedMass;

    metallicity = get_metallicity(Gal[centralgal].EjectedMass, Gal[centralgal].MetalsEjectedMass);
    dust_fraction_ejected = get_dust_fraction(Gal[centralgal].EjectedMass, Gal[centralgal].DustEjectedMass);
   
    if (Gal[centralgal].EjectedMassSN > 0.0)
    {
      channel_fractionSN = Gal[centralgal].EjectedMassSN / Gal[centralgal].EjectedMass; 
    }

    if (Gal[centralgal].EjectedMassQSO > 0.0)
    { 
      channel_fractionQSO =  Gal[centralgal].EjectedMassQSO /  Gal[centralgal].EjectedMass; 
    }

 
    Gal[centralgal].EjectedMass -= reincorporated;
    Gal[centralgal].MetalsEjectedMass -= metallicity * reincorporated;
    Gal[centralgal].DustEjectedMass -= dust_fraction_ejected * reincorporated;
    Gal[centralgal].HotGas += reincorporated;
    Gal[centralgal].MetalsHotGas += metallicity * reincorporated;
    Gal[centralgal].DustHotGas += dust_fraction_ejected * reincorporated;

    Gal[centralgal].EjectedMassSN = Gal[centralgal].EjectedMass * channel_fractionSN;
    Gal[centralgal].EjectedMassQSO = Gal[centralgal].EjectedMass * channel_fractionQSO;
    
    if (Gal[centralgal].DustEjectedMass < 0.0)
      Gal[centralgal].DustEjectedMass = 0.0;

    //printf("Reincorporated Gal %d\tHalo %d\tEjected %.4e\tSN %.4e\tQSO %.4e\tSNFrac %.4e\tQSOFrac %.4e\n", centralgal, Gal[centralgal].HaloNr, Gal[centralgal].EjectedMass, Gal[centralgal].EjectedMassSN, Gal[centralgal].EjectedMassQSO, channel_fractionSN, channel_fractionQSO);

    if (Gal[centralgal].EjectedMass > 1e-10)
      XASSERT(Gal[centralgal].EjectedMass / (Gal[centralgal].EjectedMassSN + Gal[centralgal].EjectedMassQSO) > 0.95, "EjectedMass %.4e\tSN %.4e\tQSO %.4e\tRatio %.4e\n", Gal[centralgal].EjectedMass, Gal[centralgal].EjectedMassSN, Gal[centralgal].EjectedMassQSO, Gal[centralgal].EjectedMass / (Gal[centralgal].EjectedMassSN + Gal[centralgal].EjectedMassQSO));
  }

}
