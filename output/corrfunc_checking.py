import os
import numpy as np
import h5py
import random

import Corrfunc
from Corrfunc.theory.DD import DD
from Corrfunc.theory.xi import xi 
from Corrfunc.utils import convert_3d_counts_to_cf

import AllVars
import ReadScripts


def get_snap_parts(path, snap, num_files, Npart_to_get, check_file=0):
                   
 
    if check_file:
        file_name = "./snapshot_pos.npz"
        if os.path.isfile(file_name):
            pos = np.load(file_name)

            pos_x = pos["pos_x"]
            pos_y = pos["pos_y"]
            pos_z = pos["pos_z"]

            return pos_x, pos_y, pos_x

    if Npart_to_get:

        Npart = int(Npart_to_get)
        parts_per_file = int(Npart/num_files)
        parts_remaining = Npart
        parts_added = 0

        print("Reading snapshots and grabbing a total of {0} particles. "
              "For {1} files this is approx {2} particles per "
              "file.".format(Npart, num_files, parts_per_file))
    else:

        print("Npart_to_get was not specified, using all particles at this "
              "snapshot.")

        fname = "{0}_{1:03d}/snapshot_{1:03d}.0.hdf5".format(path, snap)
        with h5py.File(fname, "r") as f_in:
            Npart = f_in["Header"].attrs["NumPart_Total"][1]

        parts_added = 0

        print("For snapshot {0} this is {1} particles.".format(snap, Npart))

    pos_x = np.full(Npart, -1, dtype=np.float64)
    pos_y = np.full(Npart, -1, dtype=np.float64)
    pos_z = np.full(Npart, -1, dtype=np.float64)

    for fnr in range(num_files):
        if fnr % 1 == 0:
            print("At file {0}, {1} particles added".format(fnr, parts_added))

        fname = "{0}_{1:03d}/snapshot_{1:03d}.{2}.hdf5".format(path, snap, fnr)

        with h5py.File(fname, "r") as f_in:
            Npart_this_file = int(f_in["Header"].attrs["NumPart_ThisFile"][1])

            if Npart_to_get:
                if fnr == num_files - 1:
                    parts_to_add = parts_remaining
                else:
                    parts_to_add = parts_per_file

                rand_inds = np.random.randint(0, Npart_this_file, parts_to_add) 
                
                for ind in rand_inds: 
                    pos = f_in["PartType1"]["Coordinates"]
  
                    pos_x[parts_added] = pos[ind,0]
                    pos_y[parts_added] = pos[ind,1]
                    pos_z[parts_added] = pos[ind,2]

                    parts_remaining -= 1 
                    parts_added += 1 

            else:
                pos_x[parts_added:parts_added+Npart_this_file] = pos[:,0]
                pos_y[parts_added:parts_added+Npart_this_file] = pos[:,1]
                pos_z[parts_added:parts_added+Npart_this_file] = pos[:,2]

                parts_added += Npart_this_file

    assert(parts_added == Npart)
    '''
    for arr in [pos_x, pos_y, pos_z]:
        w = np.where(arr < 0.0)[0]
        if len(w) > 0:
            print("Indices {0} were not filled.".format(w)) 
            raise RuntimeError
    '''
    if check_file:
        np.savez(file_name, pos_x=pos_x, pos_y=pos_y, pos_z=pos_z)
        print("Pseudo Snapshot positions saved to {0}".format(file_name))

    return pos_x, pos_y, pos_z

def get_pseudo_parts(path, snap, num_files, Npart_to_get, check_file=0):
   
    if check_file:
        file_name = "./pseudo_pos.npz"
        if os.path.isfile(file_name):
            pos = np.load(file_name)

            pos_x = pos["pos_x"]
            pos_y = pos["pos_y"]
            pos_z = pos["pos_z"]

            return pos_x, pos_y, pos_x

    if Npart_to_get:

        Npart = int(Npart_to_get)
        parts_per_file = int(Npart/num_files)
        parts_remaining = Npart
        parts_added = 0

        print("Reading pseudo snapshots and grabbing a total of {0} particles. "
              "For {1} files this is approx {2} particles per "
              "file.".format(Npart, num_files, parts_per_file))
    else:

        print("Npart_to_get was not specified, using all particles at this "
              "snapshot.")

        fname = "{0}_{1:03d}/kali_linker.0.hdf5".format(path, snap)
        with h5py.File(fname, "r") as f_in:
            Npart = f_in["Header"].attrs["NumPart_Total"][1]

        parts_added = 0

        print("For snapshot {0} this is {1} particles.".format(snap, Npart))

    pos_x = np.full(Npart, -1, dtype=np.float64)
    pos_y = np.full(Npart, -1, dtype=np.float64)
    pos_z = np.full(Npart, -1, dtype=np.float64)

    for fnr in range(num_files):
        if fnr % 20 == 0:
            print("At file {0}, {1} particles added".format(fnr, parts_added))
        fname = "{0}_{1:03d}/kali_linker.{2}.hdf5".format(path, snap, fnr)

        with h5py.File(fname, "r") as f_in:
            pos = f_in["PartType1"]["Coordinates"][:]
            Npart_this_file = int(f_in["Header"].attrs["NumPart_ThisFile"][1])

        if Npart_to_get:
            if fnr == num_files - 1:
                parts_to_add = parts_remaining
            else:
                parts_to_add = parts_per_file

            rand_inds = np.random.randint(0, len(pos), parts_to_add) 
               
            pos_x[parts_added:parts_added+parts_to_add] = pos[rand_inds,0]
            pos_y[parts_added:parts_added+parts_to_add] = pos[rand_inds,1]
            pos_z[parts_added:parts_added+parts_to_add] = pos[rand_inds,2]

            max_x = max(pos[rand_inds,0])
            max_y = max(pos[rand_inds,1])
            max_z = max(pos[rand_inds,2])

            parts_remaining -= parts_to_add
            parts_added += parts_to_add

        else:
            pos_x[parts_added:parts_added+Npart_this_file] = pos[:,0]
            pos_y[parts_added:parts_added+Npart_this_file] = pos[:,1]
            pos_z[parts_added:parts_added+Npart_this_file] = pos[:,2]

            parts_added += Npart_this_file

    assert(parts_added == Npart)
    '''
    for arr in [pos_x, pos_y, pos_z]:
        w = np.where(arr < 0.0)[0]
        if len(w) > 0:
            print("Indices {0} were not filled.".format(w)) 
            raise RuntimeError
    '''
    if check_file:
        np.savez(file_name, pos_x=pos_x, pos_y=pos_y, pos_z=pos_z)
        print("Pseudo Snapshot positions saved to {0}".format(file_name))

    return pos_x, pos_y, pos_z


def get_subfind_halo_pos(path, snap, num_files, Nhalos, check_file=0):

    if check_file:
        file_name = "./subfind_pos.npz"
        if os.path.isfile(file_name):
            pos = np.load(file_name)

            pos_x = pos["pos_x"]
            pos_y = pos["pos_y"]
            pos_z = pos["pos_z"]

            return pos_x, pos_y, pos_x

    Nhalos = int(Nhalos)
    halos_per_file = int(Nhalos/num_files)
    halos_remaining = Nhalos
    halos_added = 0

    print("Reading subfind halos and grabbing a total of {0} halos. "
          "For {1} files this is approx {2} halos per file.".format(Nhalos,
                                                                    num_files,
                                                                    halos_per_file))

    pos_x = np.empty(Nhalos, dtype=np.float64)
    pos_y = np.empty(Nhalos, dtype=np.float64)
    pos_z = np.empty(Nhalos, dtype=np.float64)

    halos_read = 0

    for fnr in range(num_files):
        fname = "{0}_{1:03d}.catalog_subgroups_properties/subfind_{1:03d}.catalog_subgroups_properties.{2}" \
                .format(path, snap, fnr)
        Halos = ReadScripts.read_subfind_halos(fname)
        Nhalos_thisfile = len(Halos)

        pos = Halos["position_COM"]

        if fnr == num_files - 1:
            halos_to_add = halos_remaining
        else:
            halos_to_add = halos_per_file

        rand_inds = np.random.randint(0, len(pos), halos_to_add) 
        
        pos_x[halos_added:halos_added+halos_to_add] = pos[rand_inds,0]
        pos_y[halos_added:halos_added+halos_to_add] = pos[rand_inds,1]
        pos_z[halos_added:halos_added+halos_to_add] = pos[rand_inds,2]

        halos_remaining -= halos_to_add
        halos_added += halos_to_add

    assert(halos_added == Nhalos)
    if check_file:        
        np.savez(file_name, pos_x=pos_x, pos_y=pos_y, pos_z=pos_z) 
        print("Subfind Halo positions saved to {0}".format(file_name)) 

    return pos_x, pos_y, pos_z 


if __name__ == "__main__":

    AllVars.Set_Params_Kali()

    #snap = 92 
    snap = 50 
    num_pseudo_files = 256
    num_subfind_files = 64
    Npart = 2e6 
    Nhalos = 5e5
    check_file = 0

    bin_min = 0.1
    bin_max = 10.0
    N_bins = 15    
    bins = np.logspace(np.log10(bin_min), np.log10(bin_max), N_bins) 
    autocorr = 0
    nthreads = 1 
    boxsize = 109.0 

    pseudo_snap_path="/fred/oz004/jseiler/kali/pseudo_snapshots/groups"
    snap_path="/fred/oz004/msinha/simulations/Kali/2400/8364/snapshots/snapshot"
    subfind_halo_path="/fred/oz004/jseiler/kali/subfind_catalogs/subfind"

    print("Grabbing {0} pseudo particles".format(Npart))
    X1, Y1, Z1 = get_pseudo_parts(pseudo_snap_path, snap, num_pseudo_files,
                                  Npart, check_file)

    #print("Grabbing {0} snapshot particles".format(Npart))    
    #X1, Y1, Z1 = get_snap_parts(snap_path, snap, num_pseudo_files,
    #                            Npart, check_file)

    print("Done")
  
    print("Calculating xi") 
    results = xi(boxsize, nthreads, bins, X1, Y1, Z1, output_ravg=True) 
    for r in results: 
        print("{0:10.6f} {1:10.6f} {2:10.6f} {3:10.6f} {4:10d} {5:10.6f}" \
              .format(r['rmin'], r['rmax'],
                      r['ravg'], r['xi'], r['npairs'], r['weightavg']))

    print("Grabbing {0} subfind halos".format(Nhalos))
    X2, Y2, Z2 = get_subfind_halo_pos(subfind_halo_path, snap,
                                      num_subfind_files, Nhalos, check_file)
    print("Done") 

    N_rand = 3*len(X1)

    rand_X = np.random.uniform(0, boxsize, N_rand)
    rand_Y = np.random.uniform(0, boxsize, N_rand)
    rand_Z = np.random.uniform(0, boxsize, N_rand)

    print("Computing Cross-Corr")
    D1D2 = DD(autocorr, nthreads, bins, X1, Y1, Z1, X2=X2, Y2=Y2, Z2=Z2,
              boxsize=boxsize)

    D1R = DD(autocorr, nthreads, bins, X1, Y1, Z1, X2=rand_X, Y2=rand_Y, 
             Z2=rand_Z, boxsize=boxsize)
    
    D2R = DD(autocorr, nthreads, bins, X2, Y2, Z2, X2=rand_X, Y2=rand_Y, 
             Z2=rand_Z, boxsize=boxsize)

    autocorr = 1 
    RR = DD(autocorr, nthreads, bins, X2, Y2, Z2, boxsize=boxsize) 

    cf = convert_3d_counts_to_cf(Npart, Nhalos, N_rand, N_rand,
                                 D1D2, D1R, D2R, RR)

    for xi in cf: 
        print("{0:10.6f}".format(xi))
