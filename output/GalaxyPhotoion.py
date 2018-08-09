import numpy as np

import ReadScripts


def calc_gal_photoion(GridPos, PhotoField_fname, GridSize, Precision):

    photoion = ReadScripts.read_binary_grid(PhotoField_fname, 
                                            GridSize,
                                            Precision)    
    photoion = np.reshape(photoion, GridSize*GridSize*GridSize)
    gal_photoion = photoion[GridPos]

    print(min(gal_photoion))

    print("There are {0} Cells with a non-zero photion and there are {1} " 
          "unique Galaxy cells".format(len(photoion[photoion > 1e-16]),
                                       len(np.unique(GridPos))))
    return gal_photoion
