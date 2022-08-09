"""
Meshes created for testing openep
=================================

MESH_2 was created from the OpenEP DATASET_2. The largest surface was extracted
and the surface then smoothed.

MESH_2_dense was created from the OpenEP DATASET_2. This mesh has had all nodes
not referenced in the triangulation removed.

"""
__all__ = [
    "MESH_2_SMOOTH",
    "MESH_2_DENSE",
]

from pkg_resources import resource_filename

MESH_2_SMOOTH = resource_filename(__name__,
                           "VTK/openep_mesh_2_smooth.vtk")

MESH_2_DENSE = resource_filename(__name__,
                           "VTK/openep_mesh_2_dense.vtk")
