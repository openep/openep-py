import pyvista
import openep
from openep._datasets.openep_datasets import DATASET_2

case = openep.load_openep_mat(DATASET_2)
mesh = case.create_mesh()
plotter = openep.draw.draw_map(mesh, field=case.fields.bipolar_voltage)
plotter.show()
