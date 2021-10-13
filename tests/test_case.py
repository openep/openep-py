from unittest import TestCase
from numpy.testing import assert_allclose

import pyvista

from openep.case import Case
from openep._simple_meshes.simple_meshes import (
    CUBE, SPHERE, BROKEN_SPHERE
)

class CaseTests(TestCase):
    def setUp(self) -> None:
        cube = pyvista.read(CUBE)
        
        self.fakefield = cube.points[:, 0]
        self.case = Case("Test", cube.points, cube.faces, {"fakefield": self.fakefield}, {}, {})

    def test_fields(self):
        assert_allclose(self.fakefield, self.case.fields['fakefield'])

    def test_mesh_creation(self):
        mesh = self.case.create_mesh()
        self.assertIsInstance(mesh, pyvista.PolyData)
        
