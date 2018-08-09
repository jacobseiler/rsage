import os
import numpy as np
import h5py
import random

import Corrfunc
from Corrfunc.theory.DD import DD
from Corrfunc.theory.xi import xi 
from Corrfunc.utils import convert_3d_counts_to_cf

import matplotlib
matplotlib.use('Agg')
import pylab as plt

import AllVars
import ReadScripts
import PlotScripts

output_format = ".png"

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


def plot_corrs(corrfunc_corr_results, corrfunc_crosscorr_results, model_tags, output_tag):

    fig1 = plt.figure(figsize=(8,8))
    ax1 = fig1.add_subplot(111)

    for model_number in range(len(corrfunc_corr_results)):
        results = corrfunc_corr_results[model_number]

        print(model_tags[model_number])
        for i in range(len(results["ravg"])):
            print("Ravg {0}\tXi {1}".format(results["ravg"][i],
                                            results["xi"][i]))

        ax1.plot(results["ravg"], results["xi"], 
                 color = PlotScripts.colors[model_number], 
                 label = model_tags[model_number]) 

        if model_number == 0:
            ax1.plot(results["ravg"], corrfunc_crosscorr_results,
                     color = 'k', ls = '--', label = "Cross-Correlation") 

    ax1.set_yscale("symlog")

    ax1.set_xlabel(r"$\mathbf{log \: r \: [Mpc \: h^{-1}]}$", 
                   fontsize = PlotScripts.global_labelsize)

    ax1.set_ylabel(r"$\mathbf{\xi (r)}$", 
                   fontsize = PlotScripts.global_labelsize)

    leg = ax1.legend(loc='upper right', numpoints=1,
                         labelspacing=0.1)
    leg.draw_frame(False)  # Don't want a box frame
    for t in leg.get_texts():  # Reduce the size of the text
        t.set_fontsize(PlotScripts.global_legendsize-2)

    plt.tight_layout()

    outputFile = "./{0}{1}".format(output_tag,
                                   output_format)

    plt.savefig(outputFile)  # Save the figure

    print('Saved file to {0}'.format(outputFile))


if __name__ == "__main__":

    AllVars.Set_Params_Kali()
    PlotScripts.Set_Params_Plot()

    #snap = 92 
    snap = 50 
    num_pseudo_files = 256
    num_subfind_files = 64
    Npart = 1e6 
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

    print("Grabbing {0} subfind halos".format(Nhalos))
    X2, Y2, Z2 = get_subfind_halo_pos(subfind_halo_path, snap,
                                      num_subfind_files, Nhalos, check_file)
   
    for i in range(len(X2)):
        X2[i] = (X2[i] + 14.28) % boxsize
        Y2[i] = (Y2[i] + 14.28) % boxsize
        Z2[i] = (Z2[i] + 14.28) % boxsize

    #print("Grabbing {0} snapshot particles".format(Npart))    
    #X1, Y1, Z1 = get_snap_parts(snap_path, snap, num_pseudo_files,
    #                            Npart, check_file)
  
    print("Calculating xi") 
    corrfunc_xi_pseudo = xi(boxsize, nthreads, bins, X1, Y1, Z1, output_ravg=True) 
    corrfunc_xi_halos = xi(boxsize, nthreads, bins, X2, Y2, Z2, output_ravg=True) 

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

    plot_corrs([corrfunc_xi_pseudo, corrfunc_xi_halos], cf,
               ["Pseudo Snapshots", "Halos"],
               "pseudo_halos") 

