from unittest import TestCase
from numpy.testing import assert_allclose

import pyvista

from openep.case.case import Case
from openep._simple_meshes.simple_meshes import CUBE


class CaseTests(TestCase):
    def setUp(self) -> None:
        cube = pyvista.read(CUBE)

        indices = cube.faces.reshape(-1, 4)[:, 1:]  # ignore the number of vertices per face
        self.fakefield = cube.points[:, 0]

        self.case = Case("Test", cube.points, indices, {"fakefield": self.fakefield}, {}, {})

    def test_fields(self):
        assert_allclose(self.fakefield, self.case.fields['fakefield'])

    def test_mesh_creation(self):
        mesh = self.case.create_mesh()
        self.assertIsInstance(mesh, pyvista.PolyData)
