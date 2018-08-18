import h5py
import numpy as np

import ReadScripts
import AllVars


def get_subfind_halos(path, snap):

    fname = "{0}_{1:03d}.catalog_subgroups_properties/subfind_{1:03d}.catalog_subgroups_properties.0" \
            .format(path, snap)
    Halos = ReadScripts.read_subfind_halos(fname)

    return Halos


def find_MBP(path, snap, num_files, MBP_ID):

    print("Searching for ID {0}".format(MBP_ID))

    for fnr in range(num_files):
        #print("Searching file {0}".format(fnr))

        fname = "{0}_{1:03d}/snapshot_{1:03d}.{2}.hdf5".format(path, snap, fnr)
        #fname = "{0}_{1:03d}/kali_linker.{2}.hdf5".format(path, snap, fnr)

        with h5py.File(fname, "r") as f_in:

            try:
                w = np.where(f_in["PartType1"]["ParticleIDs"][:] == MBP_ID)[0]
            except KeyError:
                continue

            if len(w) == 0:
                continue

            print("Found ID {0} in File {1} at index {2}".format(MBP_ID, fnr,
                                                                 w[0]))

            return f_in["PartType1"]["Coordinates"][w[0]]


if __name__ == "__main__":

    AllVars.Set_Params_Kali()

    subfind_halo_path="/fred/oz004/jseiler/kali/subfind_catalogs/subfind"
    snap_path="/fred/oz004/msinha/simulations/Kali/2400/8364/snapshots/snapshot"
    #snap_path="/fred/oz004/jseiler/kali/pseudo_snapshots/groups"
    snap = 20
    num_files=256
    full_info = 0
 
    Halos = get_subfind_halos(subfind_halo_path, snap)

    large_halo_inds = np.where(Halos["n_particles"] > 10000)[0]
    large_halo_inds = np.arange(0, len(Halos))
    
    for large_halo_idx in large_halo_inds:
        #large_halo_idx = large_halo_inds[2]
        MBP_ID = Halos["id_MBP"][large_halo_idx]

        particle_pos = find_MBP(snap_path, snap, num_files, MBP_ID)
        halo_MBP_pos = Halos["position_MBP"][large_halo_idx]

        if full_info:
            print("One the largest halos is halo {0} consisting of {1} "
                  "particles and has a COM position "
                  "{2}".format(large_halo_idx, 
                               Halos["n_particles"][large_halo_idx], 
                               Halos["position_COM"][large_halo_idx]))
            print("According to the halo file, the MBP position is "
                  "{0}".format(Halos["position_MBP"][large_halo_idx]))

            print("MBP x, y, z Positions") 
            print("Halo\tSnapshot\tDifference\tRatio")
            for i in range(3):
                halo_pos = Halos["position_MBP"][large_halo_idx,i]
                snap_pos = particle_pos[i]
                print("{0:.6f}\t{1:.6f}\t{2:.6f}\t{3:.6f}".format(halo_pos, snap_pos,
                                                                  halo_pos - snap_pos,
                                                                  halo_pos / snap_pos))

            continue


        for i in range(3):
            halo_pos = Halos["position_MBP"][large_halo_idx,i]
            snap_pos = particle_pos[i]
            if (halo_pos - snap_pos > -14.27) or (halo_pos - snap_pos < -14.29):
                # Need to check the 'wrap_around' case.
                if not (halo_pos - snap_pos < 2.92) or (halo_pos - snap_pos > 2.95):
                    print("{0:.6f}\t{1:.6f}\t{2:.6f}\t{3:.6f}".format(halo_pos, snap_pos,
                                                                      halo_pos - snap_pos,
                                                                      halo_pos / snap_pos))
                    raise ValueError
