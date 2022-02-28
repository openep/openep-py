"""
OpenEP datasets for testing openep
==================================

DATASET_2 was created from openep_dataset_2.mat file in the openep-examples
repository (https://github.com/openep/openep-examples). It was loaded into MATLAB,
the triRep object cast to a regular struct, then saved as a version 7.0 MAT file.

"""
__all__ = [
    "DATASET_2",
]

from pkg_resources import resource_filename

DATASET_2 = resource_filename(__name__,
                                  "OpenEP-MATLAB/openep_dataset_2.mat")
