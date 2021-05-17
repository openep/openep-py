'''
Module containing functions to load OpenEP dataset
'''

import numpy as np
import scipy.io as sio
import re
import h5py


def load_mat(filename, exclude_pattern=".*#.*"):
    """
    Attempt to load the given MATLAB file as either an old-style file or as an HDF5 file. The keys in the
    HDF5 file which match `exclude_pattern` will be ignored. This currently does not resolve references
    in the HDF5 file.
    """
    try:
        return sio.loadmat(filename) # try to load old-style
    except NotImplementedError:
        # ignore the exception which should only indicate the file is HDF5
        pass

    # if we got here that means the NotImplementedError was raised, try h5py
    with h5py.File(filename, "r") as f:
        dat={}

        def _visitor(key, value):
            if isinstance(value, h5py.Dataset) and re.search(exclude_pattern, key) is None:
                arr_val = np.array(value)

                # ignore references for now
                if arr_val.dtype!=object:
                    dat[key] = arr_val

        f.visititems(_visitor)  # visit all items in the file and populate dat

    return dat
