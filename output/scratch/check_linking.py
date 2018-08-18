import h5py
import numpy as np

def check_linking(snap_path, pseudo_path, snapnum, num_chunks, TypeMax=6):

    print("We are now checking the properties of each particle in the FoF linked list to ensure they match those found in the Snapshot.")
    for chunk_idx in range(0, num_chunks):

        print("Checking chunk {0}.".format(chunk_idx))
        fname_snapshot = "{0}/snapshot_{1:03d}/snapshot_{1:03d}.{2}.hdf5".format(snap_path, snapnum, chunk_idx)
        fname_pseudo = "{0}_{1:03d}/kali_linker.{2}.hdf5".format(pseudo_path, snapnum, chunk_idx)

        with h5py.File(fname_snapshot, "r") as file_snapshot, \
             h5py.File(fname_pseudo, "r") as file_pseudo: 

            snapshot_partid = file_snapshot["PartType1"]["ParticleIDs"][:]
            pseudo_partid = file_pseudo["PartType1"]["ParticleIDs"][:]

            print("There are {0} particles in the snapshot and {1} in the "
                  "pseudo snapshot".format(len(snapshot_partid),
                                           len(pseudo_partid)))
            print("Finding all the IDs from the pseudo in the snapshot.")

            inds_found_snapshot = (np.nonzero(np.in1d(snapshot_partid, pseudo_partid)))[0]
            print("Found {0} IDs in the "
                  "snapshot.".format(len(inds_found_snapshot)))
            snapshot_ids_found = snapshot_partid[inds_found_snapshot]

            print("We found {0} particles from this chunk that were also present in the pseudo_snapshot.".format(len(inds_found_snapshot)))
            assert(len(inds_found_snapshot) == len(pseudo_partid))

            for count, ind in enumerate(inds_found_snapshot):
                if(count % 10000 == 0):
                    print("Finished checking particle {0}".format(count))
                assert(file_snapshot["PartType1"]["Coordinates"][ind, 0] == file_pseudo["PartType1"]["Coordinates"][count, 0])
                assert(file_snapshot["PartType1"]["Coordinates"][ind, 1] == file_pseudo["PartType1"]["Coordinates"][count, 1])
                assert(file_snapshot["PartType1"]["Coordinates"][ind, 2] == file_pseudo["PartType1"]["Coordinates"][count, 2])

                assert(file_snapshot["PartType1"]["Velocities"][ind, 0] == file_pseudo["PartType1"]["Velocities"][count, 0])
                assert(file_snapshot["PartType1"]["Velocities"][ind, 1] == file_pseudo["PartType1"]["Velocities"][count, 1])
                assert(file_snapshot["PartType1"]["Velocities"][ind, 2] == file_pseudo["PartType1"]["Velocities"][count, 2])

    print("All particles have been accounted for.  Good job!!!")
 
if __name__ == "__main__":

    pseudo_path="/fred/oz004/jseiler/kali/pseudo_snapshots/groups"
    snap_path="/fred/oz004/msinha/simulations/Kali/2400/8364/snapshots/"
    snapnum=50  
 
    check_linking(snap_path, pseudo_path, snapnum, 256)    
