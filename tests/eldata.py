# OpenEP
# Copyright (c) 2021 OpenEP Collaborators
## This file is part of OpenEP.
## OpenEP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
## OpenEP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
## You should have received a copy of the GNU General Public License along
# with this program (LICENSE.txt). If not, see <http://www.gnu.org/licenses/>


import scipy.io as sio
from scipy.interpolate import LinearNDInterpolator as linterp
from scipy.interpolate import NearestNDInterpolator as nearest
from scipy.stats import gaussian_kde
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, LightSource
from mpl_toolkits.mplot3d import Axes3D
import trimesh


# interpolation function
'''
Parameters:
------------
x0 – Cartesian co-ordinates of the data points
d0 – Data values; voltage values
x1 – the locations at which you want data (mesh nodes)

Return:
interpolated values (voltage values at each mesh node)
'''

class openEpDataInterpolator():
    def __init__(self):
        pass
                 # method,
                 # interMethod,
                 # exterMethod,
                 # distanceThreshold,
                 # smoothingLength,
                 # fillWith,
                 # rbfConstant):

        # self.method = method
        # self.interMethod = interMethod
        # self.exterMethod = exterMethod
        # self.distanceThreshold = distanceThreshold
        # self.smoothingLength = smoothingLength
        # self.fillWith = fillWith
        # self.rbfConstant = rbfConstant

    def interpolate(self,x0,d0,x1,*args):
        self.x0 = x0
        self.d0 = d0
        self.x1 = x1
        # self.funcinterp = linterp(self.x0, self.d0)
        # self.funcnearest = nearest(self.x0, self.d0)
        # print(self.x0.shape)
        # print(self.d0.shape)
        # z = self.funcinterp(*args)
        F = linterp(self.x0,self.d0,*args)
        self.interp = F(self.x1)
        print('out-values\n',self.interp)

        # self.F =
        chk = np.isnan(self.interp)
        if chk.any():
            return np.where(chk, nearest(self.x0,self.d0,*args),self.interp)
        else:
            return self.interp

        # return self.out






class LinearNDInterpolatorExt(object):
    def __init__(self, points, values):
        self.funcinterp = linterp(points, values)
        self.funcnearest = nearest(points, values)

    def __call__(self, *args):
        z = self.funcinterp(*args)
        chk = np.isnan(z)
        if chk.any():
            return np.where(chk, self.funcnearest(*args), z)
        else:
            return z




# DEFINITIONS
amplitude = 'max'
fill_threshold = 2
volt_threshold = 2
magenta = np.array([255/255, 0/255, 255/255, 1])
magenta_1 = [255, 0, 255, 255]
# Loading Electroanatomic mapping data
file_path = '../openep_py/tests/data/openep_dataset_1.mat'



# Load the file
# Main MAT File - userdata
main_file = sio.loadmat(file_path)

# Anatomical Data - mesh descriptions(vertices & triangles)
# Vertices/Points = .X
data_tri_X = main_file['userdata']['surface_triRep_X'][0][0]
x = data_tri_X[:,0]
y = data_tri_X[:,1]
z = data_tri_X[:,2]

# Triangulation/faces - .T
data_tri_T = main_file['userdata']['surface_triRep_Triangulation'][0][0]
# Resolving indexing issue - matlab to python
t = data_tri_T-1



# Electric data (Userdata.Electric)
electric =  main_file['userdata']['electric'][0][0][0][0]

# Locations - cartesian coordinates
egmSurfX = np.array(electric[13])
print('egmSurfX\n', egmSurfX)
point_x = egmSurfX[:,0]
point_y = egmSurfX[:,1]
point_z = egmSurfX[:,2]


print('point_x-type\n',type(point_x))
from matplotlib import cm
print('egmsurfx-type\n',type(egmSurfX))

# Voltage - time series of electrogram voltages
egm = np.array(electric[4])
print('egm\n',egm)
# duration - milliseconds
egm_duration = egm.shape[1]

# Window of interest
# Annotations Struct (UserData.Electric.Annotations)
annotations = electric[10][0][0]

# Reference annotation
ref = np.array(annotations[1])
print('window of interest - ref\n',ref)

# woi
woi = np.array(annotations[0])
print('woi\n',woi)

# Extract the region of the electrogram of interest for each mapping point
region_of_interest = ref + woi
print('region of interest\n',region_of_interest)

# For each mapping point, n, find the voltage amplitude
#  - Apply the user defined function (for example max(egm) – min(egm))
max_volt=[]
min_volt=[]

for volt in egm:
    max_volt.append(np.amax(volt))
    min_volt.append(np.amin(volt))

max_volt = np.array(max_volt)
min_volt = np.array(min_volt)

#
amplitude_volt = max_volt - min_volt
print('amplitude_volt\n',amplitude_volt.shape)

# amplitude_volt = np.reshape(amplitude_volt,(1366, 1))

x0 = egmSurfX
d0 = amplitude_volt
x1 = data_tri_X


print('max_amplitude\n',np.amax(d0))
print('min_amplitude\n',np.amin(d0))


print('x0\n',x0.shape)
print('d0\n',d0.shape)
print(x1)



# INTERPOLATION
# int = openEpDataInterpolator()
# d1 = int.interpolate(x0=egmSurfX,d0=amplitude_volt,x1=data_tri_X)
# print('d1\n',np.array(d1))
F = LinearNDInterpolatorExt(points=(point_x,point_y,point_z),
                            values=d0)
d1 = F(x,y,z)

# for item in d1:
#     print(item)




# ColorShell
# coloring the shell by the voltage value rather than z-axis height
min_volt = np.amin(d1)
max_volt = max(d1)
print(min_volt)

# normalisation - normalising the voltage - values
norm = mp.colors.Normalize(vmin=min_volt, vmax=max_volt)

# defining a custom color scalar field
# if only nodal values are known, then the signal need to be averaged by the triangles
color = np.mean(d1[t], axis=1)
print('color-arry\n',color.shape)

# Voltage Threshold - to split the colorbar axis
# determining the size of the voltage values
diff_volt = max_volt - min_volt
print('diff_volt\n',diff_volt)
dist_below_threshold = volt_threshold - min_volt
dist_above_threshold = max_volt - volt_threshold

# splitting the colorbar
# color below threshold - size
size_below_threshold = int(round((dist_below_threshold/diff_volt)*256))
# color above threshold - size
size_above_threshold = int(round(256 - size_below_threshold))
print('size-above-threshold:',size_above_threshold)

# Define Colormap - Reverse JET
cm_jet = plt.cm.get_cmap('jet_r',size_below_threshold)
newcolors = cm_jet(np.linspace(0,1,256))
# Magenta - RGBA array
magenta = np.array([255/255, 0/255, 255/255, 1], dtype=np.float64)
# magenta = np.array([139/255, 0/255, 139/255, 1])
# Gray - RGBA array
gray = np.array([128/255, 128/255, 128/255, 1])

# Creating new colormap columns
# Col 1 - Jet
col1 = cm_jet(np.linspace(0,1,size_below_threshold))
print('col1\n',col1)
# Col2 - Magenta
col2 = np.zeros((size_above_threshold,4))+magenta
# col2 = magenta
print(col2)
# col3 = np.ones(size_above_threshold)
# Combining Col1 + Col2
final_col = np.concatenate((col1,col2))
print('final_col\n',final_col.shape)
# NewColormap
newcmp = ListedColormap(final_col)








mesh_3d = trimesh.Trimesh(vertices=data_tri_X,faces=t)
mesh_3d.visual.face_colors=magenta
freeboundary = mesh_3d.outline()



mesh_scene = trimesh.scene.Scene(geometry=mesh_3d)
mesh_scene1 = trimesh.scene.Scene(geometry=freeboundary)
# mesh_scene = trimesh.scene.Scene(geometry=open3d_mesh)

# trimesh.viewer.SceneViewer(scene=mesh_scene,
#                     smooth=True)
# mesh_visual = trimesh.visual.color.ColorVisuals(mesh=mesh_3d, face_colors=magenta, vertex_colors=None)

# scene = trimesh.scene(mesh_visual)
# trimesh.viewer.SceneViewer(mesh_scene)
# mesh_scene.show()
# mesh_scene1.show()















fig = plt.figure(constrained_layout=True,
                 frameon=True,
                 figsize=(40,20), dpi=80)
# # first plot
# # # get the current 3d axis on the current fig
ax1 = fig.add_subplot(1,1,1, projection='3d')
# plt.scatter(point_x,point_y,point_z,d0)

# Trisurface 3-D Mesh Plot
surf1 = ax1.plot_trisurf(x,
                         y,
                         z,
                         triangles=t,
                         linewidth=0.1,
                         antialiased=True,
                         cmap=cm.viridis
                         )

# for item in d1_new[t]:
#     print('item',item)
# color = np.array(d1[t])
color = np.mean(d1[t], axis=1)

print(color.shape)
surf1.set_array(color)
# ax1 = fig.add_subplot(1,2,2, projection='3d')
# plt.scatter(x,y,z,d1)

# min_volt_d0 = min(d0)
# print('min-volt-before-interp',min_volt_d0)
#
# max_volt_d0 = max(d0)
# print('max-volt-before-interp',max_volt_d0)

min_volt = min(d1)
print('min_volt_interp_volt\n',min_volt)
max_volt = max(d1)
print('max_volt_interp_volt\n',max_volt)
norm = mp.colors.Normalize(vmin=min_volt, vmax=max_volt)

cb = plt.colorbar(mp.cm.ScalarMappable(norm=norm, cmap=cm.viridis),
                  ax=[ax1],
                  location='left',
                  label='Voltage (mV)')

# plt.clim(0,10)

plt.show()
