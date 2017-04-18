#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')

import sys
import os
import numpy as np
import pylab as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from numpy import *
from random import sample, seed
import math
import random
import csv
from io import StringIO
from collections import Counter
from matplotlib.colors import LogNorm
import time
from matplotlib.ticker import MultipleLocator
from os.path import getsize as getFileSize
from numpy import inf 

def Read_SAGE_Objects(Model_Name, Object_Desc, Contain_TreeInfo, Dot, First_File, Last_File):
    # Initialize variables.
    TotNTrees = 0
    TotNHalos = 0
    FileIndexRanges = []
    
    print "Determining array storage requirements."
    
    # Read each file and determine the total number of galaxies to be read in
    goodfiles = 0
    for fnr in xrange(First_File, Last_File+1):
        if (Dot == 1):
            fname = Model_Name+'.'+str(fnr)  # Complete filename
        else:
            fname = Model_Name+'_'+str(fnr)  # Complete filename

        if fnr == 33:
	    continue 
        if not os.path.isfile(fname):
            print "File\t%s  \tdoes not exist!  Skipping..." % (fname)
            continue
        
        if getFileSize(fname) == 0:
            print "File\t%s  \tis empty!  Skipping..." % (fname)
            continue
        
        fin = open(fname, 'rb')  # Open the file
        if (Contain_TreeInfo == 1):
            Ntrees = np.fromfile(fin,np.dtype(np.int32),1)  # Read number of trees in file
            TotNTrees = TotNTrees + Ntrees  # Update total sim trees number
        NtotHalos = np.fromfile(fin,np.dtype(np.int32),1)[0]  # Read number of gals in file.
        TotNHalos = TotNHalos + NtotHalos  # Update total sim gals number
        goodfiles = goodfiles + 1  # Update number of files read for volume calculation
        fin.close()
    
    print
    print "Input files contain:\t%d trees ;\t%d Halos/Galaxies." % (TotNTrees, TotNHalos)
    print

    # Initialize the storage array
    G = np.empty(TotNHalos, dtype=Object_Desc)
    
    offset = 0  # Offset index for storage array
    
    # Open each file in turn and read in the preamble variables and structure.
    print "Reading in files."
    for fnr in xrange(First_File,Last_File+1):
	
        if (Dot == 1):
            fname = Model_Name+'.'+str(fnr)  # Complete filename
        else:
            fname = Model_Name+'_'+str(fnr)  # Complete filename
        
        if not os.path.isfile(fname):
            print "Couldn't find file %s" %(fname)
            continue
        
        if getFileSize(fname) == 0:
            print "File %s is empty!" %(fname)
            continue
    
        fin = open(fname, 'rb')  # Open the file
        if (Contain_TreeInfo == 1):
            Ntrees = np.fromfile(fin, np.dtype(np.int32), 1)  # Read number of trees in file
        NtotHalos = np.fromfile(fin, np.dtype(np.int32), 1)[0]  # Read number of gals in file.
        if (Contain_TreeInfo == 1):
            GalsPerTree = np.fromfile(fin, np.dtype((np.int32, Ntrees)),1) # Read the number of gals in each tree
	
        print ":   Reading N=", NtotHalos, "   \t Objects from file: ", fname
 
        GG = np.fromfile(fin, Object_Desc, NtotHalos)  # Read in the galaxy structures

        FileIndexRanges.append((offset,offset+NtotHalos))
        
        # Slice the file array into the global array
        # N.B. the copy() part is required otherwise we simply point to
        # the GG data which changes from file to file
        # NOTE THE WAY PYTHON WORKS WITH THESE INDICES!
        G[offset:offset+NtotHalos]=GG[0:NtotHalos].copy()
        
        del(GG)
        offset = offset + NtotHalos  # Update the offset position for the global array

        fin.close()  # Close the file

    print
    print "Total Objects considered:", TotNHalos
    
    # Convert the Galaxy array into a recarray
    G = G.view(np.recarray)
    
    return G

def ReadHalos(DirName, First_File, Last_File):

    Halo_Desc_full = [
    ('Descendant',          np.int32),
    ('FirstProgenitor',     np.int32),
    ('NextProgenitor',      np.int32),
    ('FirstHaloInFOFgroup', np.int32),
    ('NextHaloInFOFgroup',  np.int32),
    ('Len',                 np.int32),
    ('M_mean200',           np.float32),
    ('Mvir',                np.float32),
    ('M_TopHat',            np.float32),
    ('Pos',                 (np.float32, 3)),
    ('Vel',                 (np.float32, 3)),
    ('VelDisp',             np.float32),
    ('Vmax',                np.float32),
    ('Spin',                (np.float32, 3)),
    ('MostBoundID',         np.int64),
    ('SnapNum',             np.int32),
    ('Filenr',              np.int32),
    ('SubHaloIndex',        np.int32),
    ('SubHalfMass',         np.float32)
                     ]

    print "First_File = %d" %(First_File)
    print "Last_File = %d" %(Last_File)

    names = [Halo_Desc_full[i][0] for i in xrange(len(Halo_Desc_full))]
    formats = [Halo_Desc_full[i][1] for i in xrange(len(Halo_Desc_full))]
    Halo_Desc = np.dtype({'names':names, 'formats':formats}, align=True)

    return Read_SAGE_Objects(DirName, Halo_Desc, 1, 1, First_File, Last_File)


def ReadGals_Post_Processed_SAGE(DirName, First_File, Last_File, MAXSNAPS):

    Galdesc_full = [
                ('SnapNum'                      , np.int32),
                ('Type'                         , np.int32),
                ('GalaxyIndex'                  , np.int64),
                ('CentralGalaxyIndex'           , np.int64),
                ('SAGEHaloIndex'                , np.int32),
                ('SAGETreeIndex'                , np.int32),
                ('SimulationFOFHaloIndex'       , np.int32),
                ('mergeType'                    , np.int32),
                ('mergeIntoID'                  , np.int32),
                ('mergeIntoSnapNum'             , np.int32),
                ('dT'                           , np.float32),
                ('Pos'                          , (np.float32, 3)),
                ('Vel'                          , (np.float32, 3)),
                ('Spin'                         , (np.float32, 3)),
                ('Len'                          , np.int32),
                ('Mvir'                         , np.float32),
                ('CentralMvir'                  , np.float32),
                ('Rvir'                         , np.float32),
                ('Vvir'                         , np.float32),
                ('Vmax'                         , np.float32),
                ('VelDisp'                      , np.float32),
                ('ColdGas'                      , np.float32),
                ('StellarMass'                  , np.float32),
                ('BulgeMass'                    , np.float32),
                ('HotGas'                       , np.float32),
                ('EjectedMass'                  , np.float32),
                ('BlackHoleMass'                , np.float32),
                ('IntraClusterStars'            , np.float32),
                ('MetalsColdGas'                , np.float32),
                ('MetalsStellarMass'            , np.float32),
                ('MetalsBulgeMass'              , np.float32),
                ('MetalsHotGas'                 , np.float32),
                ('MetalsEjectedMass'            , np.float32),
                ('MetalsIntraClusterStars'      , np.float32),
                ('SfrDisk'                      , np.float32),
                ('SfrBulge'                     , np.float32),
                ('SfrDiskZ'                     , np.float32),
                ('SfrBulgeZ'                    , np.float32),
                ('DiskRadius'                   , np.float32),
                ('Cooling'                      , np.float32),
                ('Heating'                      , np.float32),
                ('QuasarModeBHaccretionMass'    , np.float32),
                ('TimeOfLastMajorMerger'         , np.float32),
                ('TimeOfLastMinorMerger'         , np.float32),
                ('OutflowRate'                  , np.float32),
                ('infallMvir'                   , np.float32),
                ('infallVvir'                   , np.float32),
                ('infallVmax'                   , np.float32),
                ('BirthTime', np.float32),
                ('GridHistory', (np.int32, MAXSNAPS)), # Array index 48
                    #('deltaEjectedMass', (np.float32, 64)),
                    #('deltaMetalsEjectedMass', (np.float32, 64)),
                    #('deltaEnergyEjected', (np.float32, 64)),
                    #('deltaSfr', (np.float32, 64)),
                ('GridStellarMass', (np.float32, MAXSNAPS)),
                ('GridSFR', (np.float32, MAXSNAPS)),
                ('GridZ', (np.float32, MAXSNAPS)),
                ('GridCentralGalaxyMass', (np.float32, MAXSNAPS))
                ]
    
    print "Reading in SAGE files (Post processed.)."
    
    names = [Galdesc_full[i][0] for i in xrange(len(Galdesc_full))]
    formats = [Galdesc_full[i][1] for i in xrange(len(Galdesc_full))]
    Gal_Desc = np.dtype({'names':names, 'formats':formats}, align=True)
    

    return (Read_SAGE_Objects(DirName, Gal_Desc, 1, 0, First_File, Last_File), Gal_Desc)

def ReadGals_SAGE_Photons(DirName, First_File, Last_File, MAXSNAPS):

    Galdesc_full = [
         ('SnapNum'                      , np.int32),
         ('Type'                         , np.int32),
         ('GalaxyIndex'                  , np.int64),
         ('CentralGalaxyIndex'           , np.int64),
         ('SAGEHaloIndex'                , np.int32),
         ('SAGETreeIndex'                , np.int32),
         ('SimulationFOFHaloIndex'       , np.int32),
         ('mergeType'                    , np.int32),
         ('mergeIntoID'                  , np.int32),
         ('mergeIntoSnapNum'             , np.int32),
         ('dT'                           , np.float32),
         ('Pos'                          , (np.float32, 3)),
         ('Vel'                          , (np.float32, 3)),
         ('Spin'                         , (np.float32, 3)),
         ('Len'                          , np.int32),
         ('Mvir'                         , np.float32),
         ('CentralMvir'                  , np.float32),
         ('Rvir'                         , np.float32),
         ('Vvir'                         , np.float32),
         ('Vmax'                         , np.float32),
         ('VelDisp'                      , np.float32),
         ('ColdGas'                      , np.float32),
         ('StellarMass'                  , np.float32),
         ('BulgeMass'                    , np.float32),
         ('HotGas'                       , np.float32),
         ('EjectedMass'                  , np.float32),
         ('BlackHoleMass'                , np.float32),
         ('IntraClusterStars'            , np.float32),
         ('MetalsColdGas'                , np.float32),
         ('MetalsStellarMass'            , np.float32),
         ('MetalsBulgeMass'              , np.float32),
         ('MetalsHotGas'                 , np.float32),
         ('MetalsEjectedMass'            , np.float32),
         ('MetalsIntraClusterStars'      , np.float32),
         ('SfrDisk'                      , np.float32),
         ('SfrBulge'                     , np.float32),
         ('SfrDiskZ'                     , np.float32),
         ('SfrBulgeZ'                    , np.float32),                  
         ('DiskRadius'                   , np.float32),                  
         ('Cooling'                      , np.float32),                  
         ('Heating'                      , np.float32),
         ('QuasarModeBHaccretionMass'    , np.float32),
         ('TimeOfLastMajorMerger'         , np.float32),
         ('TimeOfLastMinorMerger'         , np.float32),
         ('OutflowRate'                  , np.float32),
         ('infallMvir'                   , np.float32),
         ('infallVvir'                   , np.float32),
         ('infallVmax'                   , np.float32),
         ('GridHistory', (np.int32, MAXSNAPS)), # Array index 48 
         ('GridStellarMass', (np.float32, MAXSNAPS)),
         ('GridSFR', (np.float32, MAXSNAPS)),
         ('GridZ', (np.float32, MAXSNAPS)),
         ('GridCentralGalaxyMass', (np.float32, MAXSNAPS)),
         ('Photons_HI', (np.float32, MAXSNAPS)),
         ('Photons_HeI', (np.float32, MAXSNAPS)), 
         ('Photons_HeII', (np.float32, MAXSNAPS)), 
         ('MfiltGnedin', (np.dtype('d'), MAXSNAPS)),
         ('MfiltSobacchi', (np.dtype('d'), MAXSNAPS)), 
         ('EjectedFraction', (np.float32, MAXSNAPS)),  
         ('LenHistory', (np.int32, MAXSNAPS)) 
         ]
                             
    print "Reading in SAGE files (Post STARBURST)."
    print len(Galdesc_full)
    names = [Galdesc_full[i][0] for i in xrange(len(Galdesc_full))]
    formats = [Galdesc_full[i][1] for i in xrange(len(Galdesc_full))]
    Gal_Desc = np.dtype({'names':names, 'formats':formats}, align=True)  
 
    return (Read_SAGE_Objects(DirName, Gal_Desc, 1, 0, First_File, Last_File), Gal_Desc)


def ReadGals_SAGE_NoGrid(DirName, First_File, Last_File, MAXSNAPS):

    Galdesc_full = [
         ('SnapNum'                      , np.int32),
         ('Type'                         , np.int32),
         ('GalaxyIndex'                  , np.int64),
         ('CentralGalaxyIndex'           , np.int64),
         ('SAGEHaloIndex'                , np.int32),
         ('SAGETreeIndex'                , np.int32),
         ('SimulationFOFHaloIndex'       , np.int32),
         ('mergeType'                    , np.int32),
         ('mergeIntoID'                  , np.int32),
         ('mergeIntoSnapNum'             , np.int32),
         ('dT'                           , np.float32),
         ('Pos'                          , (np.float32, 3)),
         ('Vel'                          , (np.float32, 3)),
         ('Spin'                         , (np.float32, 3)),
         ('Len'                          , np.int32),
         ('Mvir'                         , np.float32),
         ('CentralMvir'                  , np.float32),
         ('Rvir'                         , np.float32),
         ('Vvir'                         , np.float32),
         ('Vmax'                         , np.float32),
         ('VelDisp'                      , np.float32),
         ('ColdGas'                      , np.float32),
         ('StellarMass'                  , np.float32),
         ('BulgeMass'                    , np.float32),
         ('HotGas'                       , np.float32),
         ('EjectedMass'                  , np.float32),
         ('BlackHoleMass'                , np.float32),
         ('IntraClusterStars'            , np.float32),
         ('MetalsColdGas'                , np.float32),
         ('MetalsStellarMass'            , np.float32),
         ('MetalsBulgeMass'              , np.float32),
         ('MetalsHotGas'                 , np.float32),
         ('MetalsEjectedMass'            , np.float32),
         ('MetalsIntraClusterStars'      , np.float32),
         ('SfrDisk'                      , np.float32),
         ('SfrBulge'                     , np.float32),
         ('SfrDiskZ'                     , np.float32),
         ('SfrBulgeZ'                    , np.float32),                  
         ('DiskRadius'                   , np.float32),                  
         ('Cooling'                      , np.float32),                  
         ('Heating'                      , np.float32),
         ('QuasarModeBHaccretionMass'    , np.float32),
         ('TimeOfLastMajorMerger'         , np.float32),
         ('TimeOfLastMinorMerger'         , np.float32),
         ('OutflowRate'                  , np.float32),
         ('infallMvir'                   , np.float32),
         ('infallVvir'                   , np.float32),
         ('infallVmax'                   , np.float32)
         ]
                             
    print "Reading in SAGE files (Post STARBURST)."
    print len(Galdesc_full)
    names = [Galdesc_full[i][0] for i in xrange(len(Galdesc_full))]
    formats = [Galdesc_full[i][1] for i in xrange(len(Galdesc_full))]
    Gal_Desc = np.dtype({'names':names, 'formats':formats}, align=True)  
 
    return (Read_SAGE_Objects(DirName, Gal_Desc, 1, 0, First_File, Last_File), Gal_Desc)

    
## Using this one ##   
def ReadGals_SAGE_DelayedSN(DirName, First_File, Last_File, MAXSNAPS):

    Galdesc_full = [
         ('SnapNum'                      , np.int32),
         ('Type'                         , np.int32),
         ('GalaxyIndex'                  , np.int64),
         ('CentralGalaxyIndex'           , np.int64),
         ('SAGEHaloIndex'                , np.int32),
         ('SAGETreeIndex'                , np.int32),
         ('SimulationFOFHaloIndex'       , np.int32),
         ('mergeType'                    , np.int32),
         ('mergeIntoID'                  , np.int32),
         ('mergeIntoSnapNum'             , np.int32),
         ('dT'                           , np.float32),
         ('Pos'                          , (np.float32, 3)),
         ('Vel'                          , (np.float32, 3)),
         ('Spin'                         , (np.float32, 3)),
         ('Len'                          , np.int32),
         ('Mvir'                         , np.float32),
         ('CentralMvir'                  , np.float32),
         ('Rvir'                         , np.float32),
         ('Vvir'                         , np.float32),
         ('Vmax'                         , np.float32),
         ('VelDisp'                      , np.float32),
         ('ColdGas'                      , np.float32),
         ('StellarMass'                  , np.float32),
         ('BulgeMass'                    , np.float32),
         ('HotGas'                       , np.float32),
         ('EjectedMass'                  , np.float32),
         ('BlackHoleMass'                , np.float32),
         ('IntraClusterStars'            , np.float32),
         ('MetalsColdGas'                , np.float32),
         ('MetalsStellarMass'            , np.float32),
         ('MetalsBulgeMass'              , np.float32),
         ('MetalsHotGas'                 , np.float32),
         ('MetalsEjectedMass'            , np.float32),
         ('MetalsIntraClusterStars'      , np.float32),
         ('SfrDisk'                      , np.float32),
         ('SfrBulge'                     , np.float32),
         ('SfrDiskZ'                     , np.float32),
         ('SfrBulgeZ'                    , np.float32),                  
         ('DiskRadius'                   , np.float32),                  
         ('Cooling'                      , np.float32),                  
         ('Heating'                      , np.float32),
         ('QuasarModeBHaccretionMass'    , np.float32),
         ('TimeOfLastMajorMerger'         , np.float32),
         ('TimeOfLastMinorMerger'         , np.float32),
         ('OutflowRate'                  , np.float32),
         ('infallMvir'                   , np.float32),
         ('infallVvir'                   , np.float32),
         ('infallVmax'                   , np.float32),
         ('GridHistory', (np.int32, MAXSNAPS)), # Array index 48 
         ('GridStellarMass', (np.float32, MAXSNAPS)),
         ('GridSFR', (np.float32, MAXSNAPS)),
         ('GridZ', (np.float32, MAXSNAPS)),
         ('GridCentralGalaxyMass', (np.float32, MAXSNAPS)),
         ('SNStars', (np.float32, MAXSNAPS)),
         #('Photons_HeI', (np.float32, MAXSNAPS)), 
         #('Photons_HeII', (np.float32, MAXSNAPS)), 
         ('MfiltGnedin', (np.dtype('d'), MAXSNAPS)),
         ('MfiltSobacchi', (np.dtype('d'), MAXSNAPS)), 
         ('EjectedFraction', (np.float32, MAXSNAPS)),  
         ('LenHistory', (np.int32, MAXSNAPS)) 
         ]
                             
    print "Reading in SAGE files (Post STARBURST)."
    print len(Galdesc_full)
    names = [Galdesc_full[i][0] for i in xrange(len(Galdesc_full))]
    formats = [Galdesc_full[i][1] for i in xrange(len(Galdesc_full))]
    Gal_Desc = np.dtype({'names':names, 'formats':formats}, align=True)  
 
    return (Read_SAGE_Objects(DirName, Gal_Desc, 1, 0, First_File, Last_File), Gal_Desc)


def ReadGals_STARBURST_SAGE(DirName, First_File, Last_File, MAXSNAPS):

    Galdesc_full = [
                ('SnapNum'                      , np.int32),
                ('Type'                         , np.int32),
                ('GalaxyIndex'                  , np.int64),
                ('CentralGalaxyIndex'           , np.int64),
                ('SAGEHaloIndex'                , np.int32),
                ('SAGETreeIndex'                , np.int32),
                ('SimulationFOFHaloIndex'       , np.int32),
                ('mergeType'                    , np.int32),
                ('mergeIntoID'                  , np.int32),
                ('mergeIntoSnapNum'             , np.int32),
                ('dT'                           , np.float32),
                ('Pos'                          , (np.float32, 3)),
                ('Vel'                          , (np.float32, 3)),
                ('Spin'                         , (np.float32, 3)),
                ('Len'                          , np.int32),
                ('Mvir'                         , np.float32),
                ('CentralMvir'                  , np.float32),
                ('Rvir'                         , np.float32),
                ('Vvir'                         , np.float32),
                ('Vmax'                         , np.float32),
                ('VelDisp'                      , np.float32),
                ('ColdGas'                      , np.float32),
                ('StellarMass'                  , np.float32),
                ('BulgeMass'                    , np.float32),
                ('HotGas'                       , np.float32),
                ('EjectedMass'                  , np.float32),
                ('BlackHoleMass'                , np.float32),
                ('IntraClusterStars'            , np.float32),
                ('MetalsColdGas'                , np.float32),
                ('MetalsStellarMass'            , np.float32),
                ('MetalsBulgeMass'              , np.float32),
                ('MetalsHotGas'                 , np.float32),
                ('MetalsEjectedMass'            , np.float32),
                ('MetalsIntraClusterStars'      , np.float32),
                ('SfrDisk'                      , np.float32),
                ('SfrBulge'                     , np.float32),
                ('SfrDiskZ'                     , np.float32),
                ('SfrBulgeZ'                    , np.float32),
                ('DiskRadius'                   , np.float32),
                ('Cooling'                      , np.float32),
                ('Heating'                      , np.float32),
                ('QuasarModeBHaccretionMass'    , np.float32),
                ('TimeOfLastMajorMerger'         , np.float32),
                ('TimeOfLastMinorMerger'         , np.float32),
                ('OutflowRate'                  , np.float32),
                ('infallMvir'                   , np.float32),
                ('infallVvir'                   , np.float32),
                ('infallVmax'                   , np.float32),
                ('GridHistory', (np.int32, MAXSNAPS)),
                #('deltaEjectedMass', (np.float32, 64)),
                #('deltaMetalsEjectedMass', (np.float32, 64)),
                #('deltaEnergyEjected', (np.float32, 64)),
                #('deltaSfr', (np.float32, 64)),
                ('GridStellarMass', (np.float32, MAXSNAPS)),
                ('GridSFR', (np.float32, MAXSNAPS)),
                ('GridZ', (np.float32, MAXSNAPS)),
                ('GridCentralGalaxyMass', (np.float32, MAXSNAPS)),
                ('BirthTime', np.float32),
                ('Photons', (np.dtype('d'), MAXSNAPS))
                ]
    
    print "Reading in SAGE files (Post STARBURST)."
    
    names = [Galdesc_full[i][0] for i in xrange(len(Galdesc_full))]
    formats = [Galdesc_full[i][1] for i in xrange(len(Galdesc_full))]
    Gal_Desc = np.dtype({'names':names, 'formats':formats}, align=True)
    
    
    return (Read_SAGE_Objects(DirName, Gal_Desc, 0, 0, First_File, Last_File), Gal_Desc)


def Join_Arrays(Array1, Array2, Desc):
    
    print "Joining arrays."

    G = np.empty(len(Array1) + len(Array2), Desc) # Create an empty array with enough space to hold both arrays.

    G[0:len(Array1)] = Array1[0:len(Array1)].copy() # Slice in the first array.
    G[len(Array1):len(Array1) + len(Array2)] = Array2[0:len(Array2)].copy() # Then append in the second array.

    G = G.view(np.recarray) # Turn into a C-like struct.

    return G

def Create_Snap_Arrays(G, NumSnaps, SnapList):

    print "Creating separate arrays for Snapshots", SnapList

    SnapCount = np.zeros((NumSnaps))

    for i in xrange(0, len(G)):
        for j in xrange(0,NumSnaps):
            if (G.GridHistory[i,j] != -1):
                SnapCount[j] += 1


