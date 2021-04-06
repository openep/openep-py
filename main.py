# Code written by Jerin Rajan on 18th Mar 2021

import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Declarations
file_path = '../openep/openep-examples-main/openep_dataset_1.mat'
data_tri_X_lst = list()
data_tri_T_lst = list()

# Load the file
# Main MAT File - userdata
main_file = sio.loadmat(file_path)

# Parsing TriRep Anatomy Data
# .X & .Triangulation

# Vertices/Points = .X
data_tri_X = main_file['userdata']['surface_triRep_X'][0][0]
x = data_tri_X[:,0]
y = data_tri_X[:,1]
z = data_tri_X[:,2]

# Triangulation/faces - .T
data_tri_T = main_file['userdata']['surface_triRep_Triangulation'][0][0]
t = data_tri_T-1

# PLot the Mesh
# Create a new fig window
fig = plt.figure()
ax = fig.gca(projection='3d')
p = ax.plot_trisurf(x,y,z, triangles=t, linewidth=0.2, antialiased=False,cmap=cm.RdYlGn)
ax.set_title('OpenEP TriRep Anatomy Data')
plt.axis('off')
# Colour Pallette
fig.colorbar(p,ax=ax)
plt.show()

# Fill Threshold - 0.5mV

# Voltage threshold



# DRAWMAP plots an OpenEP map
# userdata is a Carto data structure
def draw_map(
        userdata, *arg):
    pass
