import numpy as np

import ReadScripts


def calc_gal_photoion(GridPos, PhotoField_fname, GridSize, Precision, debug=0):

    photoion = ReadScripts.read_binary_grid(PhotoField_fname, 
                                            GridSize,
                                            Precision,
                                            reshape=False)

    gal_photoion = photoion[GridPos]

    if debug:
        print("There are {0} Cells with a non-zero photion and there are {1} " 
              "unique Galaxy cells".format(len(photoion[photoion > 1e-16]),
                                           len(np.unique(GridPos))))
    return gal_photoion


def calc_gal_zreion(GridPos, zreion_fname, GridSize, Precision, debug=0):

    zreion = ReadScripts.read_binary_grid(zreion_fname,
                                          GridSize,
                                          Precision,
                                          reshape=False)

    gal_zreion = zreion[GridPos]

    return gal_zreion 

