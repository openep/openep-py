# Code written by Jerin Rajan on 18th Mar 2021
# read mat file

import os
from pathlib import Path
import scipy.io as sio
import numpy as np
import pandas as pd
import trimesh
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib import cm
from mayavi import mlab
from scipy.spatial import Delaunay

# relative path to the sample dataset
# Main MAT File - userdata


# print('cwd:',os.getcwd())
# os.chdir('../repos')

file_path = '../openep/openep-examples-main/openep_dataset_2.mat'
data_tri_X_lst = list()
data_tri_T_lst = list()

# Load the file
main_file = sio.loadmat(file_path)

# Parsing TriRep Data
# Anatomy
# .X & .Triangulation

# Vertices/Points
# data_tri_X = pd.DataFrame(main_file['userdata']['surface_t[riRep_X'])
data_tri_X = main_file['userdata']['surface_triRep_X'][0][0]

# for item in data_tri_X:
#     data_tri_X_lst.append(item.tolist())

print('size_of_X:',len(data_tri_X_lst))
# Triangulation/faces
# data_tri_T = pd.DataFrame(main_file['userdata']['surface_triRep_Triangulation'][0][0])
data_tri_T = main_file['userdata']['surface_triRep_Triangulation'][0][0]

# for item in data_tri_T:
#     data_tri_T_lst.append(item.tolist())
#
print('size_of_T:',len(data_tri_T_lst))

# Convert to a Trimesh object
# mesh = trimesh.Trimesh(vertices=data_tri_X_lst,faces=data_tri_T_lst)
# mesh = trimesh.Trimesh(vertices=[[0, 0, 0], [0, 0, 1], [0, 1, 0]],faces=[[0, 1, 2]])
# print(dir(mesh))

x = data_tri_X[:,0]
y = data_tri_X[:,1]
z = data_tri_X[:,2]
t = data_tri_T-1


# # sample input
# vertices=[[0, 0, 0], [0, 0, 1], [0, 1, 0]]
# faces=[[0, 1, 2]]

print('x:',x)
print('y:',y)
print('z:',z)
print('t',t)
print(type(t))
print('length_of_T:', len(t))
print('max_value_x_in_T:',t.max())
print('min value of x in T:',t.min())
print('lenghth_of_x:', len(x))


# print(data_tri_X[:,0])
# triang = tri.Triangulation(x,y)

# PLot the Mesh
# Create a new fig window

#
# ax = plt.figure().add_subplot(projection='3d')
#
fig = plt.figure()
ax = fig.gca(projection='3d')
p = ax.plot_trisurf(x,y,z, triangles=t, linewidth=0.2, antialiased=True,cmap=cm.RdYlGn)
ax.set_title('OpenEP TriRep Anatomy Data')
fig.colorbar(p,ax=ax)
plt.show()



#
#
#
#
# p = ax.plot_trisurf([0, 0, 0], [0, 0, 1], [0, 1, 0], triangles=faces, linewidth=0.2, antialiased=True,cmap=cm.coolwarm)
# ax.set_title('OpenEP TriRep Anatomy Data')
# fig.colorbar(p,ax=ax)
# plt.show()



# mlab.triangular_mesh(x,y,z,t)
# mlab.show()
# Voltage plots


# Colour Pallette

# Fill Threshold - 0.5mV


# Voltage threshold









# DRAWMAP plots an OpenEP map
# userdata is a Carto data structure
def draw_map(
        userdata, *arg):
    pass





# fig1, ax1 = plt.subplots()
# ax1.set_aspect('equal')
# ax1.triplot(triang, 'bo-', lw=1)
# ax1.set_title('2-D plot')

# fig = plt.figure(figsize=plt.figaspect(0.5))
# ax = fig.add_subplot(1, 2, 1, projection='3d')
# ax.plot_trisurf(x, y, z, triangles=tri, cmap=plt.cm.Spectral)
# ax.set_zlim(-1, 1)

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# X, Y = np.meshgrid(x, y)
# Z = z

# # Plot the surface.
# surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
#
# # Add a color bar which maps values to colors.
# fig.colorbar(surf, shrink=0.5, aspect=5)
#
# plt.show()
