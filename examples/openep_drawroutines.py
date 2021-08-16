import sys
sys.path.append('/home/jra21/work/source/repos/opep')

from openep import io as openep_io
from openep import case as openep_case
from openep import mesh_routines as openep_mesh
from openep import case_routines as case_routines
from openep import draw_routines as draw
import numpy as np
import matplotlib.pyplot as plt



filename = '/home/jra21/work/source/repos/opep/examples/data/new_dataset_2.mat'

ep_case = openep_io.load_case(filename)

# DrawVoltage Map
hsurf = draw.draw_map(ep_case,
                    volt='bip',
                    freeboundary_color='black',
                    freeboundary_width=5,
                    cmap='jet_r',
                    minval=0,
                    maxval=2,
                    volt_below_color='brown', 
                    volt_above_color='magenta', 
                    nan_color='gray',
                    plot=True)

# GetAnatomicalStructure and DrawFreeBoundaries on mesh
# hsurf1 = draw.get_anatomical_structures(ep_case,plot=True)
