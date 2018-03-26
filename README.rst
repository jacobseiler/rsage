This repository contains the Self-Consistent Semi-Analytic Galaxy Evolution (SC-SAGE) model; an augmented version of the Semi-Analytic Galaxy Evolution (SAGE) model (https://github.com/darrencroton/sage).   


Any issues/questions/ranting can be sent to Jacob Seiler: jseiler@swin.edu.au 

## Adding an extra variable

These steps outline how to add an extra variable for SAGE to track for each snapshot and how to update the Python reader routine.

* Add the variable as a pointer to the `GALAXY` struct in `core_allvars.h`
* Allocate the memory for each galaxy in the `malloc_grid_arrays` function within `core_io_tree.c`
* **Important** Ensure the memory is being freed by adding the variable to the `free_grid_arrays` function in `core_io_tree.c`
* Initialize the variable in the `init_galaxy` function within `model_misc.c`
* Ensure the variable is being added to the MergedGalaxy list by adding it to the `add_galaxy_to_merger_list` function within `model_mergers.c`
* Save the variable to the output list by adding it to the `write_gridarray` function in `core_save.c`

At this point the variable is correctly being tracked and output.  Now need to update the ReadScript to properly reflect the change in the galaxy struct

* Add the extra variable to `Galdesc_full` in the `ReadGals_SAGE_DelayedSN` function within `output/ReadScript.py`. Be careful that the order of this is identical to the save order within `core_save.c`. The name used corresponds to what it will be called by the Python scripts. The data type should be identical to that used by SAGE.

The variable should now be tracked and read in properly, accessible via `Galaxy[Variable][SnapNum]` within the python reading scripts. 

Finally need to update the gridding code (post_processed_grid) to also track this new property.

* Add the variable as a pointer to the `GALAXY_GRID` struct in `core_allvars_grid.h`
* Allocate memory for the variable within the `load_gals` function in `core_io_gals.c`
* **Not Necessary, but useful** Add the variable to the print statement in `core_allvars_grid.h` under the `DEBUG_GALS` ifdef.  Useful to check that you're reading in the galaxy data correctly.
* Ensure the variable is being freed by adding it to the `free_gals` function in `core_io_gals.c`

