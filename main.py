# Code written by Jerin Rajan on 18th Mar 2021
# Functions to support OpenEp visualisation

import numpy as np
import trimesh as tm
import scipy.io as sio
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import make_interp_spline, BSpline

# Declarations
file_path = '../openep/openep-examples-main/openep_dataset_1.mat'
data_tri_X_lst = []
data_tri_T_lst = []
x_values = []
y_values = []
z_values = []
# voltage threshold to split the colorbar by solid and JET
volt_threshold = 2

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
# Resolving indexing issue - matlab to python
t = data_tri_T-1

data_act_bip = main_file['userdata']['surface'][0][0]['act_bip'][0][0]
voltage_data = data_act_bip[:,1]

# Checking for NaN values and replacing with 0 in voltage data
voltage_new = np.asarray(list(map(lambda x:0 if np.isnan(x) else x, voltage_data)))

# ColorShell
# coloring the shell by the voltage value rather than z-axis height
min_volt = min(voltage_new)
max_volt = max(voltage_new)

# normalisation - normalising the voltage - values
norm = mp.colors.Normalize(vmin=min_volt, vmax=max_volt)

# Plot the Mesh
# Create a new fig window
fig = plt.figure()

# first plot
# get the current 3d axis on the current fig
ax1 = fig.add_subplot(1,1,1, projection='3d')

# defining a custom color scalar field
# if only nodal values are known, then the signal need to be averaged by the triangles
color = np.mean(voltage_new[t], axis=1)

# Voltage Threshold - to split the colorbar axis
# determining the size of the voltage values
diff_volt = max_volt - min_volt

dist_below_threshold = volt_threshold - min_volt
dist_above_threshold = max_volt - volt_threshold

# splitting the colorbar
# color below threshold - size
size_below_threshold = round((dist_below_threshold/diff_volt)*256)
# color above threshold - size
size_above_threshold = 256 - size_below_threshold

# Define Colormap - Reverse JET
cm_jet = plt.cm.get_cmap('jet_r',size_below_threshold)
newcolors = cm_jet(np.linspace(0,1,256))
# Magenta - RGBA array
magenta = np.array([255/255, 0/255, 255/255, 1])
# Gray - RGBA array
gray = np.array([128/255, 128/255, 128/255, 1])

# Creating new colormap columns
# Col 1 - Jet
col1 = cm_jet(np.linspace(0,1,42))
# Col2 - Magenta
col2 = np.zeros((214,4)) + magenta
# Combining Col1 + Col2
final_col = np.concatenate((col1,col2))
# NewColormap
newcmp = ListedColormap(final_col)

# Points - vertices
points = data_tri_X

# Trimesh object
mesh = tm.Trimesh(vertices=points,faces=t)
# print(dir(mesh))

# print(dir(mesh.as_open3d))
# Outline Index nodes
# mesh.outline().show()
freeboundary = mesh.outline().to_dict()
vertex_index = mesh.outline().vertex_nodes
# freeboundary = mesh.outline().show()
# freeboundary

# print('freeboundary\n',freeboundary["entities"][0]['points'])
# print('freebounary\n',freeboundary)

# for item in vertex_index:
#     print('vertex_index\n',item)

# print(dir(mesh.outline()))
# print(dir(mesh.outline().entities))

mesh_outline_entities = mesh.outline().entities
# print('mesh_outline_entities\n',len(mesh_outline_entities))
freeboundary_vertex_nodes = list()
freeboundary_points = list()
for i in range(len(mesh_outline_entities)):
    freeboundary_vertex_nodes.append(freeboundary["entities"][i]['points'])

# print('freeboundary-vertex-nodes',freeboundary_vertex_nodes)



# print(vertex_index)

# freeboundary_vertex = vertex_index[0:len(freeboundary_vertex_nodes)]




# print(dir(freeboundary.show()))
# freeboundary = mesh.outline().vertex_nodes
# freeBoundary_vertices = mesh.outline().vertices

# print('mesh outline \n',mesh.outline())
# print(mesh.outline().vertices)
# print(dir(mesh.outline()))

# mesh Vertices
trimesh_points = mesh.vertices
print('freeboundary-vertex-nodes\n',freeboundary_vertex_nodes)

for i in freeboundary_vertex_nodes:
    # print('list',trimesh_points[i])
    # np.array(x_values.append(trimesh_points[i][:,0]))
    # np.array(y_values.append(trimesh_points[i][:,0]))
    # np.array(z_values.append(trimesh_points[i][:,0]))
    x_values.append(trimesh_points[i][:,0])
    print(x_values)
    y_values.append(trimesh_points[i][:,1])
    print(y_values)
    z_values.append(trimesh_points[i][:,2])
    print(z_values)


# x_values = np.array(x_values)
# y_values = np.array(y_values)
# z_values = np.array(z_values)


    # x_values_array = np.array(x_values).reshape(len(x_values),1)
    # y_values_array = np.array(y_values).reshape(len(y_values),1)
    # z_values_array = np.array(z_values).reshape(len(z_values),1)

    # freeboundary_points.append(np.concatenate((x_values,
    #                                       y_values,
    #                                       z_values),
    #                                       axis=1))

# print('x-values\n',z_values)
# print('x-values\n',x_values[0][:,0])
# print('x-values\n',np.array(x_values[0]).reshape(118,1))

freeboundary_points = []
x_values_array = []
y_values_array = []
z_values_array = []


for i in range(len(mesh_outline_entities)):
    x_values_array.append(x_values[i].reshape(len(x_values[i]),1))
    print('x-values',x_values_array)
    y_values_array.append(y_values[i].reshape(len(y_values[i]),1))
    print('y-values',y_values_array)
    z_values_array.append(z_values[i].reshape(len(z_values[i]),1))
    print('z-values',z_values_array)
    freeboundary_points.append(np.concatenate((x_values_array[i],
                                         y_values_array[i],
                                         z_values_array[i]),
                                         axis=1))
    # freeboundary_points = np.hstack((x_values[i],y_values[i],z_values[i]))
    #
    # freeboundary_points = [x_values[i],y_values[i],z_values[i]]
    # print('free-boundary-points\n',freeboundary_points)


for item in freeboundary_points:
    print('freeboundary-points-set\n',item)
    freeboundary_plot = ax1.plot(item[:,0],
                                 item[:,1],
                                 item[:,2],
                                 c='black',
                                 linewidth=2,
                                 alpha=1,
                                 zorder=10)



    # print('x-values\n',x_values[i])

# print('freeboudary-points\n',freeboundary_points)
# extract the ponts/vertices of the freeBoundary
# for item in freeboundary.vertex_nodes:
# for item in freeboundary_vertex:
#     for i in item:
#         x_values.append(trimesh_points[i][0])
#         y_values.append(trimesh_points[i][1])
#         z_values.append(trimesh_points[i][2])

# Reshape the x,y,z points in size  - (n x 1)
# x_values_array = list()
# print(len(x_values))
# for i in range(len(x_values)):
#     x_values_array[i] = np.array(x_values).reshape(len(x_values),1)
#
# x_values_array = np.array(x_values).reshape(len(x_values),1)
# y_values_array = np.array(y_values).reshape(len(y_values),1)
# z_values_array = np.array(z_values).reshape(len(z_values),1)





# xnew = np.linspace(x_values_array.min(),x_values_array.max(),len(x_values_array))
# ynew = np.linspace(y_values_array.min(),y_values_array.max(),len(y_values_array))
# znew = np.linspace(z_values_array.min(),z_values_array.max(),len(z_values_array))
#
# sp1 =make_interp_spline(xnew,ynew,znew,k=3)
# y_s



# # Join the unique edges together
# freeboundary_points = np.concatenate((x_values_array,
#                                       y_values_array,
#                                       z_values_array),
#                                       axis=1)
# print('freeboundary-points\n',freeboundary_points)


# FreeBoundary Edge Plot

# freeboundary_plot = ax1.plot(freeBoundary_vertices[:,0],
#                                 freeBoundary_vertices[:,1],
#                                 freeBoundary_vertices[:,2],
#                                 c='black',
#                                 linewidth=5,
#                                 alpha=0.9,
#                                 # marker='_'),
#                                 linestyle='-')

# freeboundary_plot = ax1.plot(freeboundary_points[:,0],
#                              freeboundary_points[:,1],
#                              freeboundary_points[:,2],
#                              c='black',
#                              linewidth=5,
                                # alpha=0.5,
                                # marker='_'),
                                # linestyle='-'
                              # )
# print(dir(freeboundary_plot))
# # Trisurface 3-D Mesh Plot
surf1 = ax1.plot_trisurf(x,
                         y,
                         z,
                         triangles=t,
                         linewidth=0.1,
                         antialiased=True,
                         cmap=newcmp,
                         alpha=0.9,
                         )
surf1.set_array(color)
ax1.set_title('OpenEP TriRep Anatomy Data')
plt.axis('off')

#
# # plt.axis('off')
# # print(dir(freeboundary_plot))
# # Colour Pallette, Position Left
# cb = plt.colorbar(mp.cm.ScalarMappable(norm=norm, cmap=newcmp),
#                   ax=[ax1],
#                   location='left',
#                   label='Voltage (mV)')

# Show plot
plt.show()


# DRAWMAP plots an OpenEP map
# userdata is a Carto data structure
def draw_map(
        userdata, *arg):
    pass


# Function to read the file
# return the userdata
def load():
    pass
