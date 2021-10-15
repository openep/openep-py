"""
Simple meshes created for testing openep
========================================

The meshes were originally created using trimesh.creation

"""
__all__ = [
    "CUBE",
    "SPHERE",
    "BROKEN_SPHERE"
]

from pkg_resources import resource_filename

CUBE = resource_filename(__name__,
                         "PLY/cube.ply")

SPHERE = resource_filename(__name__,
                           "PLY/sphere.ply")

BROKEN_SPHERE = resource_filename(__name__,
                                  "PLY/broken_sphere.ply")
