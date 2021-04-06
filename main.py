# Code written by Jerin Rajan on 18th Mar 2021

import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt

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
# Resolving indexing issue
t = data_tri_T-1

# ColorShell
# coloring the shell by the voltage value rather than z-axis height
data_act_bip = main_file['userdata']['surface'][0][0]['act_bip'][0][0]
# Get the voltage field data - Second column in act_bip
print('voltage (act_bip):\n', data_act_bip[:,1])
print('length of data:', len(data_act_bip))
print('x values:\n',x)

voltage_data = data_act_bip[:,1]
print(type(voltage_data))
# indexRemove = ~np.isnan(voltage_data)
# print(indexRemove)

# Checking for NaN values in voltage data
# lambda voltage_data: ' ' if np.isnan(voltage_data) else voltage_data
# print('voltage_data:\n', voltage_data)

# voltage_new = lambda x: '' if np.isnan(x) else x
voltage_new = np.asarray(map(lambda x:'' if np.isnan(x) else x, voltage_data))
print(type(voltage_new))
# print('index with Nan values:\n',indexRemove)



# Plot the Mesh
# Create a new fig window
fig = plt.figure()
ax1 = fig.add_subplot(111, projection='3d')

# - Tri-Surface Plot with a rainbow color map
surf1 = ax1.plot_trisurf(x,y,z, triangles=t, linewidth=0.2, edgecolor=None,antialiased=False,cmap=plt.cm.hsv)
ax1.set_title('OpenEP TriRep Anatomy Data')

# Disable axis
plt.axis('off')

# Colour Pallette, Position Left
cb = plt.colorbar(surf1,ax=[ax1],location='left',label='Voltage (mV)')

# Show plot
plt.show()


# Fill Threshold - 0.5mV

# Voltage threshold



# DRAWMAP plots an OpenEP map
# userdata is a Carto data structure
def draw_map(
        userdata, *arg):
    pass
