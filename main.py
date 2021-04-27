# Code written by Jerin Rajan on 18th Mar 2021
# Functions to support OpenEp visualisation

import numpy as np
import trimesh as tm
import scipy.io as sio
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap



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

# Outline Index nodes
freeboundary = mesh.outline().vertex_nodes

# mesh Vertices
trimesh_points = mesh.vertices

# extract the ponts/vertices of the freeBoundary
for item in freeboundary:
    for i in item:
        x_values.append(trimesh_points[i][0])
        y_values.append(trimesh_points[i][1])
        z_values.append(trimesh_points[i][2])

# Reshape the x,y,z points in size  - (n x 1)
x_values_array = np.array(x_values).reshape(len(x_values),1)
y_values_array = np.array(y_values).reshape(len(y_values),1)
z_values_array = np.array(z_values).reshape(len(z_values),1)

# # Join the unique edges together
freeboundary_points = np.concatenate((x_values_array,
                              y_values_array,
                              z_values_array),
                              axis=1)
# Trisurface 3-D Mesh Plot
surf1 = ax1.plot_trisurf(x,
                         y,
                         z,
                         triangles=t,
                         linewidth=0.5,
                         antialiased=True,
                         cmap=newcmp)
surf1.set_array(color)
ax1.set_title('OpenEP TriRep Anatomy Data')
plt.axis('off')

# FreeBoundary Edge Plot
ax1.scatter(freeboundary_points[:,0],
            freeboundary_points[:,1],
            freeboundary_points[:,2],
            c='black',
            linewidth=5,
            alpha=0.5,
            marker='.')

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
