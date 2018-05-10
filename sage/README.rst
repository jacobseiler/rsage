**Semi-Analytic Galaxy Evolution (SAGE)**

================
Saved Properties
================

The following tables lists the properties that are saved in the output files.  
Each output file first contains header information, **printed once**. After 
this header information are the properties for each galaxy.  **Pay attention**,
some properties are temporal properties and will contain values for each
snapshot (e.g., ``GridHistory`` will have an entry for each snapshot) whereas
some properties will only have a single number (e.g., ``TreeNr`` will be a
single integer).


Header
------

+-------------------+--------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| **Variable Name** |                 **Definition**                               |  **Data Type** |          **Units/Possible Values**               | **Number Entries** |
+===================+==============================================================+================+==================================================+====================+
| STEPS             | Number of substeps SAGE used between each snapshot.          | 32 bit integer.| Unitless.  Will be greater than 0.               | 1                  |
+-------------------+--------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| MAXSNAPS          | Number of snapshots in the simulation.                       | 32 bit integer.| Unitless.  Will be greater than 0.               | 1                  |
+-------------------+--------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| ZZ                | Redshift of each snapshot.                                   | 64 bit double. | Unitless.  Will be greater than (or equal to) 0. | MAXSNAPS           |
+-------------------+--------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| Hubble_h          | Hubble Parameter of the simulation.                          | 64 bit double. | Unitless. Will be between 0 and 1                | 1                  |
+-------------------+--------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| Omega             | Matter critical density parameter of the simulation.         | 64 bit double  | Unitless.                                        | 1                  |
+-------------------+--------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| OmegaLambda       | Dark energy critical density parameter of the simulation.    | 64 bit double. | Unitless.                                        | 1                  |
+-------------------+--------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| BaryonFrac        | (Cosmic) baryon fraction.                                    | 64 bit double. | Unitless.                                        | 1                  |
+-------------------+--------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| PartMass          | Mass of a single dark matter particle in the simulation.     | 64 bit double. | 1.0e10 Msun/h.                                   | 1                  |
+-------------------+--------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| BoxSize           | Side-length of the simulation box.                           | 64 bit double. | Mpc/h.                                           | 1                  |
+-------------------+--------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| GridSize          | Number of grid cells along one side.                         | 32 bit integer.| Unitless.                                        | 1                  |
+-------------------+--------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| Ntrees            | Number of trees in this file.                                | 32 bit integer.| Unitless. Will be greater than 0.                | 1                  |
+-------------------+--------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| TotGalaxies       | Number of galaxies in this file.                             | 32 bit integer.| Unitless. Will be greater than 0.                | 1                  |
+-------------------+--------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| TreeNgals         | Number of galaxies per tree in this file.                    | 32 bit integer.| Unitless. Will be greater than 0.                | Ntrees             |
+-------------------+--------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+

Galaxy Properties
-----------------

These are the properties **for each galaxy**.  The number of galaxies in this
file is given by ``TotGalaxies``.  For properties that have more than one
entry (e.g., ``GridHistory``), each entry will be the value of the property 
at that snapshot. 

If the galaxy does not exist at a snapshot (e.g., it has not formed yet) then
its value will be given by the **Default:** parameter.

+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| **Variable Name**  |                 **Definition**                                      |  **Data Type** |          **Units/Possible Values**               | **Number Entries** |
+====================+=====================================================================+================+==================================================+====================+
| TreeNr             | Tree galaxy belongs to.                                             | 32 bit integer.| Unitless.  Greater than (or equal to) 0.         | 1                  |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| GridType           | Type of galaxy.  Default: -1                                        | 32 bit integer.| Unitless.  0 = Central.                          |                    |
|                    |                                                                     |                |            1 = Satellite.                        | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| GridFoFHaloNr      | Friends-of-Friends (FoF) halo number.  Default: -1.                 | 32 bit integer.| Unitless.  Greater than 0.                       |                    |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| GridHistory        | 1-Dimensional Grid Index the galaxy is located at. Default: -1.     | 32 bit integer.| Unitless.  Greater than (or equal to) -1.        | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| GridColdGas        | Mass of the cold gas in the galaxy.  Default: 0.0.                  | 32 bit float.  | 1.0e10 Msun/h.                                   | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| GridHotGas         | Mass of the hot gas in the galaxy.  Default: 0.0.                   | 32 bit float.  | 1.0e10 Msun/h.                                   | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| GridEjectedMass    | Mass of the ejected gas in the galaxy.  Default: 0.0.               | 32 bit float.  | 1.0e10 Msun/h.                                   | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| GridDustColdGas    | Mass of the cold dust in the galaxy.  Default: 0.0.                 | 32 bit float.  | 1.0e10 Msun/h.                                   | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| GridDustHotGas     | Mass of the hot dust in the galaxy.  Default: 0.0.                  | 32 bit float.  | 1.0e10 Msun/h.                                   | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| GridDustEjectedMass| Mass of the ejected dust in the galaxy.  Default: 0.0.              | 32 bit float.  | 1.0e10 Msun/h.                                   | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| GridBHMass         | Mass of the black hole in the galaxy.  Default: 0.0.                | 32 bit float.  | 1.0e10 Msun/h.                                   | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| GridStellarMass    | Mass of the stars in the galaxy.  Default: 0.0.                     | 32 bit float.  | 1.0e10 Msun/h.                                   | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| GridSFR            | Total star formation rate (bulge + disc) of the galaxy.             |                |                                                  |                    |
|                    | Default: 0.0                                                        | 32 bit float.  | Msun/yr.                                         | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| GridZ              | Metallicity of the cold gas, i.e., m_cold_metals/m_cold_gas.        | 32 bit float.  | Unitless.  Absolute ratio (**not solar**)        | MAXSNAPS           |
|                    | Default: 0.0.                                                       |                |                                                  |                    |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| GridFoFMass        | Mass of the FoF halo the galaxy belongs to.  Default: 0.0.          | 32 bit float.  | 1.0e10 Msun/h.                                   | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| EjectedFraction    | Fraction of baryons in the ejected reservoir.  Default: 0.0.        | 32 bit float.  | Unitless.  Between 0 and 1 (inclusive).          | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| LenHistory         | Number of dark matter particles of the host halo (**not FOF halo**).|                |                                                  |                    |
|                    | Default: -1.                                                        | 32 bit integer.| Unitless.  Greater than 0.                       | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| QuasarActivity     | Denotes if there was a quasar that ejected all cold + hot gas.      |                |                                                  |                    |
|                    | Default: 0.                                                         | 32 bit integer.| Unitless.  Either 0 or 1.                        | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| QuasarSubstep      | Substep the quasar event occurred.  Default: -1.                    | 32 bit integer.| Unitless.  Between -1 and STEPS.                 | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| DynamicalTime      | Dynamical time of the host halo.  Default: 0.0.                     | 32 bit float.  | Myr.                                             | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| LenMergerGal       | If there was a merger in this galaxy's past, denotes the number of  |                |                                                  |                    |
|                    | dark matter particles of the **merging** halo.  Default: -1.        | 32 bit integer.| Unitless.  Greater than 0.                       | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| GridReionMod       | Value of the reionization modifier.  Default: -1.                   | 32 bit float.  | Unitless.  Greater than (or equal to) 0.         | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| GridNgamma_HI      | Number of HI ionizing photons.  Default: 0.0.                       | 32 bit float.  | log10(Photons/s).                                | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+
| Gridfesc           | Escape fraction.  Default: 0.0.                                     | 32 bit float.  | Unitless.  Greater than (or equal to) 0.         | MAXSNAPS           |
+--------------------+---------------------------------------------------------------------+----------------+--------------------------------------------------+--------------------+


========================
Adding An Extra Property
========================

These steps outline how to add an extra variable for SAGE to track for each snapshot and how to update the Python reader routine.

* Add the variable as a pointer to the `GALAXY` struct in `core_allvars.h`
* Allocate the memory for each galaxy in the `malloc_temporal_arrays` function within `self_consistent/common.c`
* **Important** Ensure the memory is being freed by adding the variable to the `free_temporal_arrays` function in `self_consistent/common.c`
* Initialize the variable in the `init_galaxy` function within `model_misc.c`
* Ensure the variable is being added to the MergedGalaxy list by adding it to the `add_galaxy_to_merger_list` function within `model_mergers.c`
* Save the variable to the output list by adding it to the `write_temporal_arrays` function in `self_consistent/common.c`.  **Note:** There are two macros for saving the property, `WRITE_GRID_PROPERTY` and `WRITE_CONVERTED_GRID_RPOPERTY`.  The first macro writes the property in the internal code units (think hard about what these should be).  The second allows you to write the property in different units; e.g., converting from code units to Msun/yr (as is done for star formation rate).  You will need to specify the conversion factor and the data-type of the written out property. 

At this point the variable is correctly being tracked and output.  Now need to update the ReadScript to properly reflect the change in the galaxy struct

* Add the extra variable to `Galdesc_full` in the `ReadGals_SAGE_DelayedSN` function within `output/ReadScript.py`. Be careful that the order of this is identical to the save order within `core_save.c`. The name used corresponds to what it will be called by the Python scripts. The data type should be identical to that used by SAGE.

The variable should now be tracked and read in properly, accessible via `Galaxy[Variable][SnapNum]` within the python reading scripts. 
