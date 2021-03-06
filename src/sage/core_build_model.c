#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>
#include <sys/stat.h>

#include "temporal_array.h"
#include "core_allvars.h"
#include "core_proto.h"

int32_t update_quasar_tracking(int32_t gal, int32_t step, float substep_dt);

void construct_galaxies(int halonr, int tree, int filenr)
{
  static int halosdone = 0;
  int prog, fofhalo, ngal;
 
  HaloAux[halonr].DoneFlag = 1;
  halosdone++;

  prog = Halo[halonr].FirstProgenitor;

  while(prog >= 0)
  {
    if(HaloAux[prog].DoneFlag == 0)
      construct_galaxies(prog, tree, filenr);
    prog = Halo[prog].NextProgenitor;
  }

  fofhalo = Halo[halonr].FirstHaloInFOFgroup;
  if(HaloAux[fofhalo].HaloFlag == 0)
  {
    HaloAux[fofhalo].HaloFlag = 1;
    while(fofhalo >= 0)
    {
      prog = Halo[fofhalo].FirstProgenitor;
      while(prog >= 0)
      {
        if(HaloAux[prog].DoneFlag == 0)
          construct_galaxies(prog, tree, filenr);
        prog = Halo[prog].NextProgenitor;
      }

      fofhalo = Halo[fofhalo].NextHaloInFOFgroup;
    }
  }

  // At this point, the galaxies for all progenitors of this halo have been
  // properly constructed. Also, the galaxies of the progenitors of all other 
  // halos in the same FOF group have been constructed as well. We can hence go
  // ahead and construct all galaxies for the subhalos in this FOF halo, and
  // evolve them in time. 

  fofhalo = Halo[halonr].FirstHaloInFOFgroup;
  if(HaloAux[fofhalo].HaloFlag == 1)
  {
    ngal = 0;
    HaloAux[fofhalo].HaloFlag = 2;

    while(fofhalo >= 0)
    {
      ngal = join_galaxies_of_progenitors(tree, fofhalo, ngal, filenr);
      fofhalo = Halo[fofhalo].NextHaloInFOFgroup;
    }

    evolve_galaxies(Halo[halonr].FirstHaloInFOFgroup, ngal, tree);
  }

}

int join_galaxies_of_progenitors(int treenr, int halonr, int ngalstart, int32_t filenr)
{
  int ngal, prog, i, j, first_occupied, lenmax, lenoccmax, centralgal;
  double previousMvir, previousVvir, previousVmax;
  int step;

  lenmax = 0;
  lenoccmax = 0;
  first_occupied = Halo[halonr].FirstProgenitor;
  prog = Halo[halonr].FirstProgenitor;

  if(prog >=0)
    if(HaloAux[prog].NGalaxies > 0)
    lenoccmax = -1;

  // Find most massive progenitor that contains an actual galaxy
  // Maybe FirstProgenitor never was FirstHaloInFOFGroup and thus has no galaxy

  while(prog >= 0)
  {
    if(Halo[prog].Len > lenmax)
    {
      lenmax = Halo[prog].Len;      
    }
    if(lenoccmax != -1 && Halo[prog].Len > lenoccmax && HaloAux[prog].NGalaxies > 0)
    {
      lenoccmax = Halo[prog].Len;
      first_occupied = prog;
    }
    prog = Halo[prog].NextProgenitor;
  }

  ngal = ngalstart;
  prog = Halo[halonr].FirstProgenitor;

  while(prog >= 0)
  {
    for(i = 0; i < HaloAux[prog].NGalaxies; i++)
    {
      if(ngal == (FoF_MaxGals-1)) 
      {
        printf("Current FoF_MaxGals = %d. Reallocing.\n", FoF_MaxGals);
        FoF_MaxGals += 10000;
        Gal = myrealloc(Gal, FoF_MaxGals * sizeof(struct GALAXY), (FoF_MaxGals - 1000) * sizeof(struct GALAXY));
      }
	
      assert(ngal < FoF_MaxGals);

      // This is the cruical line in which the properties of the progenitor galaxies 
      // are copied over (as a whole) to the (temporary) galaxies Gal[xxx] in the current snapshot 
      // After updating their properties and evolving them 
      // they are copied to the end of the list of permanent galaxies HaloGal[xxx] 

      Gal[ngal] = HaloGal[HaloAux[prog].FirstGalaxy + i];
      HaloGal[HaloAux[prog].FirstGalaxy + i].IsMerged = 0; // This is done because there is an if check within 'free_galaxies_and_tree' that checks for IsMerged != -1.
                                                           // This galaxy exists at multiple redshifts so we need to be careful that we only free it once.
      Gal[ngal].HaloNr = halonr;

      Gal[ngal].dT = -1.0;

      // this deals with the central galaxies of (sub)halos 
      if(Gal[ngal].Type == 0 || Gal[ngal].Type == 1)
      {
        // this halo shouldn't hold a galaxy that has already merged; remove it from future processing
        if(Gal[ngal].mergeType != 0)
        {
          Gal[ngal].Type = 3;
          continue;
        }

        // remember properties from the last snapshot
        previousMvir = Gal[ngal].Mvir;
        previousVvir = Gal[ngal].Vvir;
        previousVmax = Gal[ngal].Vmax;

        if(prog == first_occupied)
        {
          // update properties of this galaxy with physical properties of halo 
          Gal[ngal].MostBoundID = Halo[halonr].MostBoundID;

          for(j = 0; j < 3; j++)
          {
            Gal[ngal].Pos[j] = Halo[halonr].Pos[j];
            Gal[ngal].Vel[j] = Halo[halonr].Vel[j];
          }
					
          Gal[ngal].Len = Halo[halonr].Len;
          Gal[ngal].Vmax = Halo[halonr].Vmax;

					Gal[ngal].deltaMvir = get_virial_mass(halonr) - Gal[ngal].Mvir;

          if(get_virial_mass(halonr) > Gal[ngal].Mvir)
          {
            Gal[ngal].Rvir = get_virial_radius(halonr);  // use the maximum Rvir in model
            Gal[ngal].Vvir = get_virial_velocity(halonr);  // use the maximum Vvir in model 
          }
          Gal[ngal].Mvir = get_virial_mass(halonr);

          Gal[ngal].Cooling = 0.0;
          Gal[ngal].Heating = 0.0;
          Gal[ngal].QuasarModeBHaccretionMass = 0.0;
          Gal[ngal].OutflowRate = 0.0;

          for(step = 0; step < STEPS; step++)
          {
            Gal[ngal].SfrDisk[step] = Gal[ngal].SfrBulge[step] = 0.0;
            Gal[ngal].SfrDiskColdGas[step] = Gal[ngal].SfrDiskColdGasMetals[step] = 0.0;
            Gal[ngal].SfrBulgeColdGas[step] = Gal[ngal].SfrBulgeColdGasMetals[step] = 0.0;
          }

          if(halonr == Halo[halonr].FirstHaloInFOFgroup)
          {
            // a central galaxy
            Gal[ngal].mergeType = 0;
            Gal[ngal].mergeIntoID = -1;
            Gal[ngal].MergTime = 999.9;            

            Gal[ngal].DiskScaleRadius = get_disk_radius(halonr, ngal);

            Gal[ngal].Type = 0;
          }
          else
          {
            // a satellite with subhalo
            Gal[ngal].mergeType = 0;
            Gal[ngal].mergeIntoID = -1;

            if(Gal[ngal].Type == 0)  // remember the infall properties before becoming a subhalo
            {
              Gal[ngal].infallMvir = previousMvir;
              Gal[ngal].infallVvir = previousVvir;
              Gal[ngal].infallVmax = previousVmax;
            }

            if(Gal[ngal].Type == 0 || Gal[ngal].MergTime > 999.0)
              // here the galaxy has gone from type 1 to type 2 or otherwise doesn't have a merging time.
              Gal[ngal].MergTime = estimate_merging_time(halonr, Halo[halonr].FirstHaloInFOFgroup, ngal);
            
            Gal[ngal].Type = 1;
          }
        }
        else
        {
          // an orhpan satellite galaxy - these will merge or disrupt within the current timestep
          Gal[ngal].deltaMvir = -1.0*Gal[ngal].Mvir;
          Gal[ngal].Mvir = 0.0;

          if(Gal[ngal].MergTime > 999.0 || Gal[ngal].Type == 0)
          {
            // here the galaxy has gone from type 0 to type 2 - merge it!
            Gal[ngal].MergTime = 0.0;
          
            Gal[ngal].infallMvir = previousMvir;
            Gal[ngal].infallVvir = previousVvir;
            Gal[ngal].infallVmax = previousVmax;
          }

          Gal[ngal].Type = 2;
        }
      }

      ngal++;

    }

    prog = Halo[prog].NextProgenitor;
  }

  if(ngal == 0)
  {
    // We have no progenitors with galaxies. This means we create a new galaxy. 
    init_galaxy(ngal, halonr, treenr, filenr);
    ngal++;
  }

  // Per Halo there can be only one Type 0 or 1 galaxy, all others are Type 2  (orphan)
  // In fact, this galaxy is very likely to be the first galaxy in the halo if 
	// first_occupied==FirstProgenitor and the Type0/1 galaxy in FirstProgenitor was also the first one 
  // This cannot be guaranteed though for the pathological first_occupied!=FirstProgenitor case 

  for(i = ngalstart, centralgal = -1; i < ngal; i++)
  {
    if(Gal[i].Type == 0 || Gal[i].Type == 1)
    {
			assert(centralgal == -1);
      centralgal = i;
    }
  }

  for(i = ngalstart; i < ngal; i++)
  {
    Gal[i].CentralGal = centralgal;
    XASSERT((Gal[centralgal].Type == 0 || Gal[centralgal].Type == 1) && (Gal[centralgal].mergeType == 0), "The central galaxy of this galaxy must have Type 0 (a proper central) or Type 1 (a satellite with a subhalo) and have mergeType 0 (will not merge in the current timestep).\nFor galaxy %d in Halo %d (there are %d galaxies inside this halo), the central galaxy instead has Type %d and mergeType %d\n", i, halonr, ngal, Gal[centralgal].Type, Gal[centralgal].mergeType); 

  }


  return ngal;

}


void evolve_galaxies(int halonr, int ngal, int tree)	// Note: halonr is here the FOF-background subhalo (i.e. main halo) 
{
  int p, i, step, centralgal, merger_centralgal, currenthalo, offset;
  double infallingGas, coolingGas, deltaT, substep_dt, time, galaxyBaryons, currentMvir;
  
  centralgal = Gal[0].CentralGal;

  XASSERT(Gal[centralgal].Type == 0 && (Gal[centralgal].HaloNr == halonr), "We require that Gal[centralgal].Type = 0 and Gal[centralgal].HaloNr = halonr.\nWe have Gal[centralgal].Type = %d, Gal[centralgal].HaloNr = %d, halonr = %d, centralgal = %d\n", Gal[centralgal].Type, Gal[centralgal].HaloNr, halonr, centralgal);
  XASSERT(Gal[centralgal].mergeType == 0 || Gal[centralgal].Type == 1, "The central galaxy should not merge within this timestep (mergeType = 0) - if it does it is a satellite galaxy (Type = 1).\n The central galaxy is %d and has mergeType %d with Type %d\n", centralgal, Gal[centralgal].mergeType, Gal[centralgal].Type); 

  infallingGas = infall_recipe(centralgal, ngal, ZZ[Halo[halonr].SnapNum], halonr);

  // We first want to determine how much delayed SN feedback needs to be applied for each galaxy over this evolution step.

  for(p = 0; p < ngal; p++)
  {
    deltaT = Age[Gal[p].SnapNum] - Age[Halo[halonr].SnapNum];
    if(IRA == 0 && (Gal[p].Total_SN_SF_Time + (deltaT * UnitTime_in_Megayears / Hubble_h) > TimeResolutionSN))
    {
      do_previous_SN(p, centralgal, deltaT);
    }

    // If we're using the Quasar fescPrescription, then quasars can go off in the middle of a snapshot.
    // In this case the fesc is boosted for the rest of the snapshot.  However after that snapshot, we want to put the boost to the full value.
    // For galaxies that had a quasar in a previous snapshot, their boost value will be fractional. 
    if (fescPrescription == 2)
    {
      if (Gal[p].QuasarActivityToggle > 0 && Gal[p].QuasarFractionalPhotons < 1.0)
      {
        Gal[p].QuasarFractionalPhotons = 1.0;
      }
      if (Gal[p].QuasarActivityToggle == 0)
      {
        Gal[p].QuasarFractionalPhotons = 0.0;
      }
    }
  }

  // We integrate things forward by using a number of intervals equal to STEPS 
  for(step = 0; step < STEPS; step++)
  {

    // Loop over all galaxies in the halo 
    for(p = 0; p < ngal; p++)
    {
      // Don't treat galaxies that have already merged 
      if(Gal[p].mergeType > 0)
        continue;
      deltaT = Age[Gal[p].SnapNum] - Age[Halo[halonr].SnapNum];
      substep_dt = deltaT / STEPS;

      time = Age[Gal[p].SnapNum] - (step + 0.5) * substep_dt; 

      Gal[p].Age = Age[Gal[p].SnapNum] - Age[Halo[halonr].SnapNum];

      if(Gal[p].dT < 0.0)
        Gal[p].dT = deltaT;

      // If we're doing the quasar boosted fescPrescription, update the tracking.
      if (fescPrescription == 2)
      {
        update_quasar_tracking(p, step, substep_dt);
      } 

      // For the central galaxy only 
      if(p == centralgal)
      {
        add_infall_to_hot(centralgal, infallingGas / STEPS);
        Gal[p].GridInfallRate[Gal[p].SnapNum] += infallingGas / STEPS / substep_dt; 
        if(ReIncorporationFactor > 0.0)
          reincorporate_gas(centralgal, substep_dt); 
      }
			else 
				if(Gal[p].Type == 1 && Gal[p].HotGas > 0.0)
					strip_from_satellite(halonr, centralgal, p);

      // Determine the cooling gas given the halo properties 
      coolingGas = cooling_recipe(p, substep_dt); 
      cool_gas_onto_galaxy(p, coolingGas);

      starformation_and_feedback(p, centralgal, time, substep_dt, halonr, step, tree, ngal);

      for(int32_t snap = 0; snap < MAXSNAPS; ++snap) {
        if(Gal[p].MUV[snap] > -0.001 && Gal[p].MUV[snap] < 0.001) {
            fprintf(stderr, "In Evolve: Snap %d has MUV %.4e\n", snap, Gal[p].MUV[snap]);
        }
      }

    }

    // check for satellite disruption and merger events 
    for(p = 0; p < ngal; p++)
    {
      
      if((Gal[p].Type == 1 || Gal[p].Type == 2) && Gal[p].mergeType == 0)  // satellite galaxy!
      {
				assert(Gal[p].MergTime < 999.0);

        deltaT = Age[Gal[p].SnapNum] - Age[Halo[halonr].SnapNum];
        substep_dt = deltaT / STEPS;
        Gal[p].MergTime -= substep_dt; 

        // only consider mergers or disruption for halo-to-baryonic mass ratios below the threshold
        // or for satellites with no baryonic mass (they don't grow and will otherwise hang around forever)
        currentMvir = Gal[p].Mvir - Gal[p].deltaMvir * (1.0 - ((double)step + 1.0) / (double)STEPS);
        galaxyBaryons = Gal[p].StellarMass + Gal[p].ColdGas;

        if((galaxyBaryons == 0.0) || (galaxyBaryons > 0.0 && (currentMvir / galaxyBaryons <= ThresholdSatDisruption)))        
        {
          if(Gal[p].Type==1) 
            merger_centralgal = centralgal;
          else
            merger_centralgal = Gal[p].CentralGal;
          if(Gal[merger_centralgal].mergeType > 0)
          {         
            merger_centralgal = Gal[merger_centralgal].CentralGal;
          }

          Gal[p].mergeIntoID = NumGals + merger_centralgal;  // position in output 


          XASSERT(p != merger_centralgal, "A galaxy is trying to merge into itself!\nGalaxy = %d \t Halonr = %d \t Snapshot = %d\n", p, halonr, Gal[p].SnapNum);
          if(Gal[p].MergTime > 0.0)  // disruption has occured!
          {

            disrupt_satellite_to_ICS(merger_centralgal, p, tree);
            update_temporal_array(p, halonr, step); // Updates the temporal arrays before it's added to the merger list.
            add_galaxy_to_merger_list(p);
	
          }
          else
          {
            if(Gal[p].MergTime <= 0.0)  // a merger has occured! 
            {
              time = Age[Gal[p].SnapNum] - (step + 0.5) * substep_dt;

              deal_with_galaxy_merger(p, merger_centralgal, centralgal, time, substep_dt, halonr, step, tree, ngal);
              update_temporal_array(p, halonr, step); // Updates the temporal arrays before it's added to the merger list.
              add_galaxy_to_merger_list(p);

            }
          } 
        }
        
      }
    }

  } // Go on to the next STEPS substep

  for(p = 0; p < ngal; ++p)
  { 
    if (Gal[p].mergeType > 0) 
    {
      continue; // Merged galaxies have already had their grid properties updated.
    }
    update_temporal_array(p, halonr, STEPS); // Otherwise update the temporal properties of the galaxy.
  }

  // Extra miscellaneous stuff before finishing this halo
	Gal[centralgal].TotalSatelliteBaryons = 0.0;
  deltaT = Age[Gal[0].SnapNum] - Age[Halo[halonr].SnapNum];
	
  for(p = 0; p < ngal; p++)
  {

    // Don't bother with galaxies that have already merged 
    if(Gal[p].mergeType > 0)
      continue;
		
    Gal[p].Cooling /= deltaT;
    Gal[p].Heating /= deltaT;
		
    if(p != centralgal)
			Gal[centralgal].TotalSatelliteBaryons += 
				(Gal[p].StellarMass + Gal[p].BlackHoleMass + Gal[p].ColdGas + Gal[p].HotGas);

  }

  // Attach final galaxy list to halo 
  offset = 0;
  for(p = 0, currenthalo = -1; p < ngal; p++)
  {
    if(Gal[p].HaloNr != currenthalo)
    {      
      currenthalo = Gal[p].HaloNr;
      HaloAux[currenthalo].FirstGalaxy = NumGals;
      HaloAux[currenthalo].NGalaxies = 0;
    }

    // Merged galaxies won't be output. So go back through its history and find it
    // in the previous timestep. Then copy the current merger info there.
    offset = 0;
    i = p-1;
    while(i >= 0)
    {
     if(Gal[i].mergeType > 0) 
       if(Gal[p].mergeIntoID > Gal[i].mergeIntoID)
         offset++;  // these galaxies won't be kept so offset mergeIntoID below
     i--;
    }
    
    i = -1;
    if(Gal[p].mergeType > 0)
    {
      i = HaloAux[currenthalo].FirstGalaxy - 1;
      while(i >= 0)
      {
        if(HaloGal[i].GalaxyNr == Gal[p].GalaxyNr)
          break;
        else
          i--;
      }
      
			assert(i >= 0);
      
      HaloGal[i].mergeType = Gal[p].mergeType;
      HaloGal[i].mergeIntoID = Gal[p].mergeIntoID - offset;
      HaloGal[i].mergeIntoSnapNum = Halo[currenthalo].SnapNum;
    }
    
    if(Gal[p].mergeType == 0)
    {
      if(NumGals > MaxGals)
      {
        printf("ngal = %d \t MaxGals = %d\n", NumGals, MaxGals);
        assert(NumGals < MaxGals);
      }

      Gal[p].SnapNum = Halo[currenthalo].SnapNum;
      HaloGal[NumGals++] = Gal[p];
      HaloAux[currenthalo].NGalaxies++;      
    }
  }
}


int32_t update_quasar_tracking(int32_t gal, int32_t step, float substep_dt)
{
  // substep_dt is in code units, as is TargetQuasarTime.


  if (Gal[gal].QuasarActivityToggle > 0)
  {
    Gal[gal].QuasarBoostActiveTime += substep_dt; 
    if (Gal[gal].QuasarBoostActiveTime > Gal[gal].TargetQuasarTime)
    {
      Gal[gal].QuasarActivityToggle -= 1;
      if (Gal[gal].QuasarActivityToggle == 0) // There can be the case where a quasar goes off while the galaxy has it's fesc still boosted. Only turn off boost when ALL quasars have been switched off. 
      {
        Gal[gal].QuasarFractionalPhotons = (step + 1.0) / STEPS; // Turn the quasar off partway through the snapshot.
                                                               // We do the +1.0 because we measure from the END of a substep.
                                                               // i.e., if the quasar boost turns off when step == 4, it was active for 50% of the snapshot, not 40%.
        Gal[gal].QuasarBoostActiveTime = 0.0;
        Gal[gal].TargetQuasarTime = 0.0;  
      }
    }
  }

  return EXIT_SUCCESS;
 
}
