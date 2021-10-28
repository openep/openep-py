import pyvista
import openep
from openep._datasets.simple_meshes import TRIANGLES

filename = "/Users/paul/github/openep-py/examples/data/new_dataset_2.mat"

case = openep.load_case(filename)
mesh = case.create_mesh()

mesh = pyvista.read(TRIANGLES)
free_boundaries = openep.mesh.get_free_boundaries(mesh)
print(f"Perimeter lengths: {free_boundaries.calculate_lengths()}")
print(f"Cross-sectional areas: {free_boundaries.calculate_areas()}")

plotter = openep.draw.draw_map(mesh, field=case.fields['bip'])
plotter = openep.draw.draw_free_boundaries(free_boundaries, plotter=plotter)
plotter.show()
