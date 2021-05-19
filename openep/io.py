'''
Module containing functions to load OpenEP dataset
'''

import numpy as np
import scipy.io as sio
import re
import h5py

from .case import Case


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


def load_case(filename, name=None, exclude_pattern=".*#.*"):
    """
    Load a Case object from the given file. This assumes a number of names for objects found in the file.
    """
    dat = load_mat(filename, exclude_pattern)

    if name is None:
        name = os.path.basename(filename)

    nodes = dat["userdata/surface/triRep/X"].T
    inds = dat["userdata/surface/triRep/Triangulation"].T - 1
    act, bip = dat["userdata/surface/act_bip"]
    uni, imp, frc = dat["userdata/surface/uni_imp_frc"]

    fields = dict(act=act, bip=bip, uni=uni, imp=imp, frc=frc)
    electric = {}
    rf = {}
    other_data = {}

    for k, v in dat.items():
        _, section, *objname = k.split("/", 2)

        if objname and section in ("electric", "rf"):
            dobj = rf if section == "rf" else electric
            dobj[objname[0]] = v
        else:
            other_data[k] = v

    return Case(name, nodes, inds, fields, electric, rf, other_data)
