# Code written by Jerin Rajan on 18th Mar 2021
# Functions to support OpenEp visualisation

import scipy.io as sio
from scipy.spatial import ConvexHull, Delaunay
from scipy import ndimage as ndi
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from skimage import feature
import trimesh as tm
# from matplotlib.tri import Triangulation


# Declarations
file_path = '../openep/openep-examples-main/openep_dataset_1.mat'
data_tri_X_lst = []
data_tri_T_lst = []
x_values = []
y_values = list()
z_values = list()
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
# get the current 3d axis on the current fig
ax1 = fig.gca(projection='3d')


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
print('x-y-z points\n',points.shape)


# Trimesh object
mesh = tm.Trimesh(vertices=points,faces=t)

# All edges
print('total number of edges:\n',mesh.edges.shape)



# Unique edges
unique_edges = mesh.edges_unique
print('unique edges of the mesh\n',unique_edges.shape)

# Identify the unshared edges
unshared_edges = mesh.face_adjacency_unshared
print('unshared\n',unshared_edges.shape)

unique_edges_in_unshared = []
for item in unshared_edges:
    if item in unique_edges:
        unique_edges_in_unshared.append(item)
    else:
        pass
unique_edges_in_unshared = np.array(unique_edges_in_unshared)
print('unique edges in unshared:\n',unique_edges_in_unshared)

# Identify the unique unshared edges
uniques = {}
count = 0
for a,b in unique_edges_in_unshared:
    if (a,b) not in uniques and (b,a) not in uniques:
        uniques[(a,b)]=(a,b)
    else:
        if (a,b) in uniques:
            del uniques[(a,b)]
        else:
            del uniques[(b,a)]
# print('uniques\n',uniques)

unique_unshared_edges = np.array(list(uniques.values()))
print('unique_unshared_edges:\n',unique_unshared_edges.shape)

# Get the vertex points for unique unshared edges


# extract the ponts/vertices of the unique unshared edges
for item in unique_unshared_edges:
    for i in item:
        x_values.append(points[i][0])
        y_values.append(points[i][1])
        z_values.append(points[i][2])
        # print('unshared_edge_coordinates(x);\n',x_values.append(points[i][0]))
        # print('unshared_edge_coordinates(y);\n',y_values.append(points[i][1]))
        # print('unshared_edge_coordinates(z);\n',z_values.append(points[i][2]))

# Get the unique edge points
x_values_array = np.array(x_values).reshape(len(x_values),1)
y_values_array = np.array(y_values).reshape(len(y_values),1)
z_values_array = np.array(z_values).reshape(len(z_values),1)


xx = np.array(x_values).reshape(1,len(x_values))
yy = np.array(y_values).reshape(1,len(y_values))
zz = np.array(z_values).reshape(1,len(z_values))
print('xx-values\n',xx)
print('yy-values\n',yy.shape)
print('zz-values\n',zz.shape)


print('x_values\n',x_values_array.shape)
print('y_values\n',y_values_array.shape)
print('z_values\n',z_values_array.shape)

print('x-values\n',x)
print('y-values\n',y.shape)
print('z-values\n',z.shape)

# # Join the unique edges together
unique_unshared_edge_points = np.concatenate((x_values_array,
                                              y_values_array,
                                              z_values_array),
                                              axis=1)
print('border_points:\n',unique_unshared_edge_points.shape)

# tri = Triangulation(x_values_array,y_values_array,z_values_array)

surf = ax1.plot_trisurf(xx[0],
                       yy[0],
                       zz[0],
                       linewidth=0.1)

# mesh1 = tm.Trimesh(vertices=unique_unshared_edge_points,faces=t)
# #
# mesh1.show()

# Plotting the border to the axis
# ax1.plot_trisurf(xx[0],
#                  yy[0],
#                  zz[0],
#                  triangles=t,
#                  linewidth=0.5)

surf1 = ax1.plot_trisurf(x,
                         y,
                         z,
                         triangles=t,
                         linewidth=0.5,
                         antialiased=True,
                         cmap=newcmp)


# print('dir:\n',dir(surf1))
# print('typeofsurf',type(surf1))
# surf1.set_edgecolor([0,0,0,1])
surf1.set_array(color)
# surf1.set_markeredgecolor([0,0,0,1])
# surf1.set_fc([0.0,1,0])
ax1.set_title('OpenEP TriRep Anatomy Data')
plt.axis('off')

# Colour Pallette, Position Left
cb = plt.colorbar(mp.cm.ScalarMappable(norm=norm, cmap=newcmp),
                  ax=[ax1],
                  location='left',
                  label='Voltage (mV)')

# Show plot
plt.show()



# DRAWMAP plots an OpenEP map
# userdata is a Carto data structure
def draw_map(
        userdata, *arg):
    pass
