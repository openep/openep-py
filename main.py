# Code written by Jerin Rajan on 18th Mar 2021

import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import sklearn.preprocessing as sk
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
print(t.shape)


data_act_bip = main_file['userdata']['surface'][0][0]['act_bip'][0][0]
voltage_data = data_act_bip[:,1]
print('voltage data raw',voltage_data)
# indexRemove = ~np.isnan(voltage_data)
# print(indexRemove)

# Checking for NaN values and replacing with 0 in voltage data
voltage_new = np.asarray(list(map(lambda x:0 if np.isnan(x) else x, voltage_data)))



# ColorShell
# coloring the shell by the voltage value rather than z-axis height
min_volt = min(voltage_new)
max_volt = max(voltage_new)
print('minimum Voltage value:',min_volt)
print('maximum voltage value:',max_volt)


# scaling the z - values to the voltage values
# z_scaled = np.asarray(list(map(lambda x: x-)))


# normalisation
norm = mp.colors.Normalize(vmin=min_volt, vmax=max_volt)

# Plot the Mesh
# Create a new fig window
fig = plt.figure()
# ax1 = fig.add_subplot(111, projection='3d')
ax1 = fig.gca(projection='3d')
# my_col = plt.cm.jet(voltage_new)

# defining a custom color scalar field
# if only nodal values are known, then the signal need to be averaged by the triangles
color = np.mean(voltage_new[t], axis=1)


print('colorshape:',color.shape )
print('voltage\n',voltage_new.shape)
print('voltage_new[t]\n',voltage_new[t])
print('t\n',t)
print('color\n',color)


# cmap = plt.get_cmap('jet')
# triang = mp.tri.Triangulation(x,y,t)
# surf = ax1.plot_trisurf(triang,z,antialiased=True,cmap=plt.cm.jet)
# surf.set_array(color)


# - Tri-Surface Plot with a rainbow color map
surf1 = ax1.plot_trisurf(x,y,z,triangles=t,linewidth=0.2, edgecolor=[0,0,0,0],antialiased=True,cmap=plt.cm.jet)
# surf1 = ax1.plot_surface(x,y,z,traingles=t,linewidth=0.2,antialiased=True,cmap=plt.cm.jet)
surf1.set_array(color)


# surf1.set_edgecolor([0,0,0,1])

ax1.set_title('OpenEP TriRep Anatomy Data')

plt.axis('off')


# Colour Pallette, Position Left
# cb = plt.colorbar(surf1,ax=[ax1],location='left',label='Voltage (mV)')
cb = plt.colorbar(mp.cm.ScalarMappable(norm=norm, cmap=plt.cm.jet),ax=[ax1],location='left',label='Voltage (mV)')
# cl.LightSource(azdeg=-60,altdeg=30)

# Show plot
plt.show()


# Fill Threshold - 0.5mV

# Voltage threshold



# DRAWMAP plots an OpenEP map
# userdata is a Carto data structure
def draw_map(
        userdata, *arg):
    pass
