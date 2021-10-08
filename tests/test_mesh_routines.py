from unittest import TestCase

import numpy as np

import trimesh
from openep.case import Case
from openep.mesh_routines import (
    calculate_field_area, calculate_mesh_volume,
    calculate_vertex_distance, calculate_vertex_path
)


class TestMeshRoutines(TestCase):
    def setUp(self) -> None:
        self.cube = trimesh.creation.box()
        self.sphere = trimesh.creation.icosphere()
        self.broken_sphere = trimesh.Trimesh(self.sphere.vertices, self.sphere.faces[30:])  # sphere with holes

        xfield = self.sphere.vertices[:, 1]

        self.sphere_case = Case("Sphere", self.sphere.vertices, self.sphere.faces, {"xfield": xfield}, {}, {})
        self.broken_case = Case("BrokenSphere", self.broken_sphere.vertices, self.broken_sphere.faces, {}, {}, {})

    def test_calculate_mesh_volume(self):
        vol = calculate_mesh_volume(self.sphere_case, False)

        self.assertAlmostEqual(self.sphere.volume, vol, 2)

        vol = calculate_mesh_volume(self.sphere, False)  # test again with Trimesh object

        self.assertAlmostEqual(self.sphere.volume, vol, 2)

    def test_calculate_mesh_volume_repair(self):
        vol = calculate_mesh_volume(self.broken_case, True)

        self.assertAlmostEqual(self.sphere.volume, vol, 2)

    def test_calculate_field_area(self):
        xbelow0 = self.sphere.triangles[..., 0].mean(axis=1) <= 0  # select every triangle with mean x below YZ plane
        calculated_area = self.sphere.area_faces[xbelow0].sum()

        area = calculate_field_area(self.sphere_case, self.sphere_case.fields["xfield"], 0)

        self.assertAlmostEqual(calculated_area, area, 2)

    def test_calculate_vertex_distance1(self):
        test_dist = calculate_vertex_distance(self.cube, 0, 7)  # far corners are at indices 0 and 7
        self.assertAlmostEqual(test_dist, np.sqrt(3), 5)

    def test_calculate_vertex_distance2(self):
        # arbitrary test values
        start_idx = 10
        end_idx = 33

        start, end = self.sphere.vertices[[start_idx, end_idx]]
        dist = np.linalg.norm(start - end)

        test_dist = calculate_vertex_distance(self.sphere, start_idx, end_idx)

        self.assertAlmostEqual(test_dist, dist, 5)

    def test_calculate_vertex_path(self):
        path = calculate_vertex_path(self.cube, 0, 7)

        self.assertListEqual(path.tolist(), [0, 1, 7])
