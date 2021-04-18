# Code written by Jerin Rajan on 18th Mar 2021
# Functions to support OpenEp visualisation

import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


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
# print(t.shape)


data_act_bip = main_file['userdata']['surface'][0][0]['act_bip'][0][0]
voltage_data = data_act_bip[:,1]
# print('voltage data raw',voltage_data)
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
# z_scaled = np.asarray(list(map(lambda x: x-)))cm.get_cmap('jet')


# normalisation
norm = mp.colors.Normalize(vmin=min_volt, vmax=max_volt)

# Plot the Mesh= plt.cm.ge
# Create a new fig window
fig = plt.figure()
# ax1 = fig.add_subplot(111, projection='3d')
ax1 = fig.gca(projection='3d')


# defining a custom color scalar field
# if only nodal values are known, then the signal need to be averaged by the triangles
color = np.mean(voltage_new[t], axis=1)

# print('colorshape:',color.shape )
# print('voltage\n',voltage_new)
# print('voltage_new[t]\n',voltage_new[t])
# print('t\n',t)
# print('color\n',color)




# Define Colormap
cm_jet = plt.cm.get_cmap('jet_r',42)
newcolors = cm_jet(np.linspace(0,1,256))
# print('newcolors array\n',newcolors.shape)
magenta = np.array([255/255, 0/255, 255/255, 1])

newcolors[50:100,:] = magenta
newcmp = ListedColormap(newcolors)

threshold = 2

tAboveCAxis = np.linspace(0,1,256)
tbelowCAxis = np.linspace(0,1,256)

diff_volt = max_volt - min_volt
d_below_threshold = threshold - min_volt
d_above_threshold = max_volt - threshold

scale_below_threshold = (d_below_threshold/diff_volt)*256
scale_above_threshold = 256 - scale_below_threshold
print('scale-below_threshold',type(round(scale_below_threshold)))
print('scale-above-threshold',round(scale_above_threshold))

scale = round((d_below_threshold/diff_volt),2)
print(scale)

col1 = cm_jet(np.linspace(0,1,42))

print(col1.shape)

# newcol = cm_jet(col1)
# col22 = np.array(magenta)
# print('col22\n',col22)

# col2 = np.linspace(0.3,1,214)
col2 = np.zeros((214,4)) + magenta
print(col2.shape)
col33 = np.concatenate((col1,col2))
# print('magenta\n',magenta)

# col3 = np.linspace(0,1,256)
print(col1)
print(col1.shape)
print('col2\n',col2)
print(col33)
print(col33.shape)
newcmp1 = ListedColormap(col33)
# newcmp2 = newcmp1.colors[::-1]



# tAboveCAxis = np.linspace(0,1,scale__threshold)



# cm_jet = plt.cm.get_cmap('jet',256)
# print('jet colors\n',cm_jet(range(256)))
# print(voltage_new.size)
# newcolors = cm_jet(voltage_new)
# magenta = np.array([255, 0, 255, 1])
# for item in newcolors:
#     print(item)


# newcolors[:25, :] = pink
# newcmp = ListedColormap(newcolors)

# cmap = plt.get_cmap('jet')
# triang = mp.tri.Triangulation(x,y,t)
# surf = ax1.plot_trisurf(triang,z,antialiased=True,cmap=plt.cm.jet)
# surf.set_array(color)


# Fill Threshold - 0.5mV
# Qualitative color map - JET

# mycmap = mp.cm.get_cmap('viridis',12)
# mycmap = mp.cm.get_cmap('jet')
# print('cmap size:',mycmap)
# print('mycmap.colors\n', mycmap.colors)
# print('mycmap(range(12))\n', mycmap(range(12)))
# print('mycmap(np.linspace(0, 1, 12))\n', mycmap(np.linspace(0, 1, 12)))












# - Tri-Surface Plot with a rainbow color map
surf1 = ax1.plot_trisurf(x,y,z,triangles=t,linewidth=5, edgecolor=[0,0,0,0],antialiased=True,cmap=newcmp1)
# surf1 = ax1.plot_surface(x,y,z,traingles=t,linewidth=0.2,antialiased=True,cmap=plt.cm.jet)
surf1.set_array(color)


# surf1.set_edgecolor([0,0,0,1])


ax1.set_title('OpenEP TriRep Anatomy Data')
plt.axis('off')
# Colour Pallette, Position Left
# cb = plt.colorbar(surf1,ax=[ax1],location='left',label='Voltage (mV)')
cb = plt.colorbar(mp.cm.ScalarMappable(norm=norm, cmap=newcmp1),ax=[ax1],location='left',label='Voltage (mV)')

# # Thresholding
# tAboveCAxis = zeros(size(d));
# tBelowCAxis = zeros(size(d));

# taboveaxis is >2mv
# tbelowaxis is <2mv



# cl.LightSource(azdeg=-60,altdeg=30)


# fix color scale according to the logic in colorshell.m

# interpolated shading of faces

# drawing the free boundary


# Show plot
plt.show()





# Voltage threshold



# DRAWMAP plots an OpenEP map
# userdata is a Carto data structure
def draw_map(
        userdata, *arg):
    pass
