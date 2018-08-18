import numpy as np


def get_LHalo_datastruct():
    """
    Generates the LHalo numpy structured array.

    Parameters
    ----------

    None.

    Returns
    ----------

    LHalo_Desc: numpy structured array.  Required.
        Structured array for the LHaloTree data format.
    """

    LHalo_Desc_full = [
        ('Descendant',          np.int32),
        ('FirstProgenitor',     np.int32),
        ('NextProgenitor',      np.int32),
        ('FirstHaloInFOFgroup', np.int32),
        ('NextHaloInFOFgroup',  np.int32),
        ('Len',                 np.int32),
        ('M_Mean200',           np.float32),
        ('Mvir',                np.float32),
        ('M_TopHat',            np.float32),
        ('Posx',                np.float32), 
        ('Posy',                np.float32), 
        ('Posz',                np.float32), 
        ('Velx',                np.float32), 
        ('Vely',                np.float32), 
        ('Velz',                np.float32), 
        ('VelDisp',             np.float32),
        ('Vmax',                np.float32),
        ('Spinx',               np.float32),
        ('Spiny',               np.float32),
        ('Spinz',               np.float32),
        ('MostBoundID',         np.int64),
        ('SnapNum',             np.int32),
        ('Filenr',              np.int32),
        ('SubHaloIndex',        np.int32),
        ('SubHalfMass',         np.float32)
                         ]

    names = [LHalo_Desc_full[i][0] for i in range(len(LHalo_Desc_full))]
    formats = [LHalo_Desc_full[i][1] for i in range(len(LHalo_Desc_full))]
    LHalo_Desc = np.dtype({'names': names, 'formats': formats}, align=True)

    return LHalo_Desc


def shift_halos(path_in, path_out, num_files):

    LHalo_Struct = get_LHalo_datastruct()


    for fnr in range(num_files):

        fname_in = "{0}_{1:03d}.dat".format(path_in, fnr) 
        fname_out = "{0}_{1:03d}.dat".format(path_out, fnr) 

        print("Shifting Tree {0}".format(fnr))

        with open(fname_in, "rb") as f_in, \
             open(fname_out, "wb") as f_out:
            NTrees = np.fromfile(f_in, np.dtype(np.int32), 1)[0]
            NHalos = np.fromfile(f_in, np.dtype(np.int32), 1)[0]
            NHalosPerTree = np.fromfile(f_in,
                                        np.dtype((np.int32, NTrees)), 1)[0]

            # Write the header.
            f_out.write(np.array(NTrees, dtype=np.int32).tobytes()) 
            f_out.write(np.array(NHalos, dtype=np.int32).tobytes()) 
            f_out.write(np.array(NHalosPerTree, dtype=np.int32).tobytes())

            # Go through the input file and write those selected trees.
            for tree in range(NTrees):
                halos = np.fromfile(f_in, LHalo_Struct,
                                    NHalosPerTree[tree])
                for name in ["Posx", "Posy", "Posz"]:
                    halos[name] = (halos[name] + 14.28) % 108.96

                f_out.write(halos.tobytes())

if __name__ == "__main__":
 
    path_in = "/fred/oz004/jseiler/kali/trees/subgroup_trees"
    path_out = "/fred/oz004/jseiler/kali/shifted_trees/subgroup_trees"
    num_files = 64

    shift_halos(path_in, path_out, num_files)
