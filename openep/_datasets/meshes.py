"""
Meshes created for testing openep
=================================

MESH_2 was created from the OpenEP DATASET_2. The largest surface was extracted
and the surface then smoothed.

"""
__all__ = [
    "MESH_2",
]

from pkg_resources import resource_filename

MESH_2 = resource_filename(__name__,
                           "VTK/openep_mesh_2_smooth.vtk")
