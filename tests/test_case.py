from unittest import TestCase
import trimesh
from openep.case import Case


class CaseTests(TestCase):
    def setUp(self) -> None:
        box = trimesh.creation.box()
        field = box.vertices[:, 0]

        self.case = Case("Test", box.vertices, box.faces, {"fakefield": field}, {}, {})

    def test_mesh_creation(self):
        mesh = self.case.create_mesh()

        self.assertTrue(mesh._kwargs["parent_obj"] is self.case)
