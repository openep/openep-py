from unittest import TestCase

import numpy as np
import pyvista

from openep.case import Case
from openep.mesh_routines import (
    calculate_field_area, calculate_mesh_volume,
    calculate_vertex_distance, calculate_vertex_path
)
from openep._simple_meshes.simple_meshes import (
    CUBE, SPHERE, BROKEN_SPHERE
)


class TestMeshRoutines(TestCase):
    def setUp(self) -> None:
        self.cube = pyvista.read(CUBE)
        self.sphere = pyvista.read(SPHERE)
        self.broken_sphere = pyvista.read(BROKEN_SPHERE)  # sphere with holes

        self._sphere_faces = self.sphere.faces.reshape(-1, 4)[:, 1:]
        self._sphere_triangles = self.sphere.points[self._sphere_faces]
        self._sphere_areas = self.sphere.compute_cell_sizes(
            length=False,
            area=True,
            volume=False,
        )['Area']

        test_field = self.sphere.points[:, 1]

        self.sphere_case = Case("Sphere", self.sphere.points, self.sphere.faces, {"test_field": test_field}, {}, {})
        self.broken_case = Case("BrokenSphere", self.broken_sphere.points, self.broken_sphere.faces, {}, {}, {})

    def test_calculate_mesh_volume(self):
        vol = calculate_mesh_volume(self.sphere, False)

        self.assertAlmostEqual(self.sphere.volume, vol, 2)

    def test_calculate_mesh_volume_repair(self):
        vol = calculate_mesh_volume(self.broken_sphere, True)

        self.assertAlmostEqual(self.sphere.volume, vol, 2)

    def test_calculate_field_area(self):

        xbelow0 = self._sphere_triangles[..., 0].mean(axis=1) <= 0  # select every triangle with mean x below YZ plane
        calculated_area = self._sphere_areas[xbelow0].sum()

        area = calculate_field_area(self.sphere, self.sphere_case.fields["test_field"], 0)

        self.assertAlmostEqual(calculated_area, area, 2)

    def test_calculate_vertex_distance1(self):
        test_dist = calculate_vertex_distance(self.cube, 0, 7)  # far corners are at indices 0 and 7
        self.assertAlmostEqual(test_dist, np.sqrt(3), 5)

    def test_calculate_vertex_distance2(self):
        # arbitrary test values
        start_idx = 10
        end_idx = 33

        start, end = self.sphere.points[[start_idx, end_idx]]
        dist = np.linalg.norm(start - end)

        test_dist = calculate_vertex_distance(self.sphere, start_idx, end_idx)

        self.assertAlmostEqual(test_dist, dist, 5)

    def test_calculate_vertex_path(self):
        path = calculate_vertex_path(self.cube, 0, 7)

        self.assertListEqual(path.tolist(), [0, 1, 7])
