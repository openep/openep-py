"""
Simple meshes created for testing openep
========================================

CUBE, SPHERE, and BROKEN_SPHERE were originally created using trimesh.creation

TRIANGLES is a set of fives triangles. The first four are connected such that
they form a unit square. The fifth is disconnected from the others.

"""
__all__ = [
    "CUBE",
    "SPHERE",
    "BROKEN_SPHERE",
    "TRIANGLES",
]

from pkg_resources import resource_filename

CUBE = resource_filename(__name__,
                         "PLY/cube.ply")

SPHERE = resource_filename(__name__,
                           "VTK/sphere.vtk")

BROKEN_SPHERE = resource_filename(__name__,
                                  "PLY/broken_sphere.ply")

TRIANGLES = resource_filename(__name__,
                              "PLY/triangles.ply")
