import numpy as np
import h5py


def get_halo_fields():

    Halodesc_full = [
                ('FirstProg'                    , np.int64),
                ('NextProg'                     , np.int64),
	        ('FirstHaloInFOFgroup'          , np.int64),
                ('NextHaloInFOFgroup'           , np.int64),
                ('Pos'                          , (np.float64, 3)),
                ('Vel'                          , (np.float64, 3)),
                ('Spin'                         , (np.float64, 3)),
                ('Len'                          , np.int64),
                ('Mvir'                         , np.float64),
                ('Vvir'                         , np.float64),
                ('Vmax'                         , np.float64), 
               	('SnapNum'			, np.int64)
                ]


    print "Getting the halo fields."

    names = [Halodesc_full[i][0] for i in xrange(len(Halodesc_full))]
    formats = [Halodesc_full[i][1] for i in xrange(len(Halodesc_full))]
    Halodesc = np.dtype({'names':names, 'formats':formats})

    return Halodesc

def write_halos(filename, Halos, TotNHalos):

    output_file = open(filename, 'wb')   

    TreeNHalos = np.array([TotNHalos], dtype = np.int32)
    header_array = np.array([1, TotNHalos, TreeNHalos], dtype = np.int32) 
 
    header_array.tofile(output_file)
    Halos.tofile(output_file) 
    output_file.close()

    print "Saved the halos to file", filename

def write_header(filename, TotNHalos):

    output_file = open(filename, 'ab')

    TreeNHalos = np.array([TotNHalos], dtype = np.int32)
    header_array = np.array([1, TotNHalos, TreeNHalos], dtype = np.int32) 

    print "Writing header.  Total trees: 1, TotNHalos:", TotNHalos, "TreeNHalos:", TreeNHalos
    print array
    output_file.seek(0,0)

    header_array.tofile(output_file)
    
    output_file.close() 

def update_progenitors(Halos, Galaxies, prog_to_add, descendant_galaxy):
	
	current_first_progenitor = Halo[descendant_galaxy]['FirstProg']

	if (Galaxies['Mvir'][prog_to_add] > Galaxies['Mvir'][current_first_progenitor]): # If the virial mass of the progenitor we are adding to the chain is greater than the current 'FirstProgenitor',
		Halo[prog_to_add]['NextProg'] = current_first_progenitor # Move the current FirstProgenitor to the NextProgenitor of the progenitor we are adding,
		Halo[descendant_galaxy]['FirstProg'] = prog_to_add # Then add the new galaxy as the first progenitor.
	else:
		greater_mass = 0
		next_progenitor_in_chain = Halo[current_first_progenitor]['NextProg']
		current_progenitor_in_chain = current_first_progenitor 
		while(greater_mass == 0):
			if (next_progenitor_in_chain == -1):
				greater_mass = 1
			elif (Galaxies['Mvir'][prog_to_add] > Galaxies['Mvir'][next_progenitor_in_chain]):
				greater_mass = 1
			else:
				current_progenitor_in_chain = next_progenitor_in_chain 
				next_progenitor_in_chain = Halo[next_progenitor_in_chain]['NextProg']
	
		Halo[current_progenitor_in_chain]['NextProg'] = prog_to_add	
		Halo[prog_to_add]['NextProg'] = next_progenitor_in_chain
			
def update_fof(Halos, Galaxies, halo_to_add, central_galaxy):

	current_first_fof = Halo[central_galaxy]['FirstHaloInFOFgroup']
	assert(current_first_fof == central_galaxy)

	greater_mass = 0
	next_halo_in_chain = Halo[current_first_fof]['NextHaloInFOFgroup']
	current_halo_in_chain = current_first_fof
	while(greater_mass == 0):
		if (next_halo_in_chain == -1):
			greater_mass = 1
		elif (Galaxies['Mvir'][halo_to_add] > Galaxies['Mvir'][next_halo_in_chain]):
			greater_mass = 1
		else:
			current_halo_in_chain = next_halo_in_chain
			next_halo_in_chain = Halo[next_halo_in_chain]['NextProg']

	Halo[current_halo_in_chain]['NextHaloInFOFgroup'] = halo_to_add
	Halo[halo_to_add]['NextHaloInFOFgroup'] = next_halo_in_chain

if __name__ == '__main__':

    Halodesc = get_halo_fields()

    filepath = '/home/msinha/scratch/tao/data_products/output/meraxes/tiamat/Tiamat_meraxes_fixed_merge_intosnapnum.h5'

    with h5py.File(filepath, 'r') as f:

	tree_example = 46624 
 
    	example_gal_idx = np.arange(f['tree_displs'][tree_example], f['tree_displs'][tree_example + 1], dtype = np.int64)
    	print "Creating a tree for galaxies", example_gal_idx

	example_gal_idx_noghost = example_gal_idx[np.where(f['galaxies']['GhostFlag'][example_gal_idx] == 0)[0]]
	print "The nonghost galaxies here are", example_gal_idx_noghost
	#example_gal_idx = example_gal_idx_noghost	

    	Halos = np.full((len(example_gal_idx)), -1, dtype = Halodesc) 
	Galaxies = f['galaxies'][:][example_gal_idx]

	TotNHalos = 0

	for i in xrange(0, len(example_gal_idx)):
	
		print "Doing Halo", i 

		if(Galaxies['GhostFlag'][i] == 1): # If the galaxy we are considering is a ghost,
			continue # Just move on to the next galaxy.

		## Descendant/Progenitors Calculation ##

		ghost_galaxy_descendant = 1
		descendant_galaxy = Galaxies['descendant'][i]

		while(ghost_galaxy_descendant == 1):
			if(Galaxies['GhostFlag'][descendant_galaxy] == 0):
				break
			descendant_galaxy = Galaxies['descendant'][descendant_galaxy] # Galaxy index for the descendant of this galaxy.	

		if (descendant_galaxy != -1): # A galaxy with descendant == -1 is a root galaxy.	
			if (Halos[descendant_galaxy]['FirstProg'] == -1): # If we have yet to assign the descendant galaxy a First Progenitor,
				Halos[descendant_galaxy]['FirstProg'] = i; # Assign it's first progenitor to be the current galaxy.
			else:
				update_progenitors(Halos, Galaxies, i, descendant_galaxy) # Otherwise we need to update the pointers to be mass ordered.

		## FOF Group Claculation ##
		central_galaxy = Galaxies['CentralGal'][i] # Galaxy index for the central galaxy of the FoF group.
	
		assert(Galaxies['GhostFlag'][central_galaxy] == 0) # We assume that a normal non-ghost galaxy cannot have a ghost-galaxy as it's centralgal.  

		Halos[i]['FirstHaloInFOFgroup'] = central_galaxy 

		if (central_galaxy != i): # If the galaxy is its own central galaxy, then we don't need to update NextHaloInFOFgroup.
			if (Halos[central_galaxy]['NextHaloInFOFgroup'] == -1): # If we have yet to assign NextHaloInFOFgroup, assign it to the current galaxy.
				Halos[central_galaxy]['NextHaloInFOFgroup'] = i
			else: 
				update_fof(Halos, Galaxies, i, central_galaxy) # Otherwise update the pointers to be mass ordered.

		Halos[i]['Pos'][0] = Galaxies['posx'][i]
		Halos[i]['Pos'][1] = Galaxies['posy'][i]
		Halos[i]['Pos'][2] = Galaxies['posz'][i]

		Halos[i]['Vel'][0] = Galaxies['velx'][i]
		Halos[i]['Vel'][1] = Galaxies['vely'][i]
		Halos[i]['Vel'][2] = Galaxies['velz'][i]
	
		Halos[i]['Spin'][0] = Galaxies['Spin'][i]
		Halos[i]['Spin'][1] = Galaxies['Spin'][i]
		Halos[i]['Spin'][2] = Galaxies['Spin'][i]

		Halos[i]['Len'] = Galaxies['Len'][i]
		Halos[i]['Mvir'] = Galaxies['Mvir'][i]
		Halos[i]['Vvir'] = Galaxies['Vvir'][i]
		Halos[i]['Vmax'] = Galaxies['Vmax'][i]
		Halos[i]['SnapNum'] = Galaxies['snapnum'][i]

		TotNHalos += 1


    	print Halos
    	print Halos['FirstProg']
    	filename = 'test_halos.0'
	quit()

    #write_halos(filename, Halo, TotNHalos)
    #write_header(filename, TotNHalos)
