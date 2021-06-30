from openep import io as openep_io
from openep import case as openep_case
from openep import mesh_routines as openep_mesh
from openep import case_routines as case_routines
from openep import draw_routines as draw

import pyvista as pv
import numpy as np
from matplotlib.cm import jet, rainbow, jet_r, seismic
import trimesh as tm


filename = '/home/jra21/work/source/repos/opep/examples/data/new_dataset_1.mat'
np.set_printoptions(suppress=True)


ep_case = openep_io.load_case(filename)
ep_case_mesh = ep_case.create_mesh()

# hsurf = draw.DrawMap(ep_case,
#                     freeboundary_color='black',
#                     freeboundary_width=5,
#                     cmap='jet_r',
#                     minval=0,
#                     maxval=2,
#                     volt_below_color='brown', 
#                     volt_above_color='magenta', 
#                     nan_color='gray',
#                     plot=True)

hsurf1 = draw.getAnatomicalStructures(ep_case,plot=False)