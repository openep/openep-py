from openep import io as openep_io
from openep import case as openep_case
from openep import mesh_routines as openep_mesh
from openep import case_routines as case_routines
from openep import draw_routines as draw

import pyvista as pv
import numpy as np
from matplotlib.cm import jet, rainbow, jet_r, seismic
import trimesh as tm


filename = '/home/jra21/work/source/repos/opep/examples/data/new_dataset_2.mat'
np.set_printoptions(suppress=True)


ep_case = openep_io.load_case(filename)
ep_case_mesh = ep_case.create_mesh()


draw.DrawMap(ep_case,surf_color='magenta',freeboundary_color='black',freeboundary_width=5)


