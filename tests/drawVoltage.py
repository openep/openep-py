import numpy as np
import trimesh as tm
import scipy.io as sio
import matplotlib as mp
import matplotlib.pyplot as plt
from scipy.interpolate import NearestNDInterpolator as nearInt
from scipy.interpolate import LinearNDInterpolator as linearInt


# DEFINITIONS
file_path = '../openep_py/tests/data/openep_dataset_1.mat'
data_tri_X_lst = []
data_tri_T_lst = []
x_values = []
y_values = []
z_values = []
default_color = 'gray'

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



# Voltage Data
data_act_bip = main_file['userdata']['surface'][0][0]['act_bip'][0][0]
voltage_data = data_act_bip[:,1]
print('voltag_bip\n',voltage_data.shape)

voltage_new = np.reshape(a=voltage_data,newshape=(len(voltage_data),1))
print('voltage-new\n',voltage_new.shape)

# Checking for NaN values and replacing with 0 in voltage data
voltage_new = np.asarray(list(map(lambda x:0 if np.isnan(x) else x, voltage_data)))


# DRAW THE SURFACE
# Create a new fig window
fig = plt.figure(num=1,
                 figsize=[10,8],
                 facecolor='#F0F0F0',
                 # constrained_layout=True,
                 frameon=True,
                 tight_layout=False,
                 dpi=500)



# get the current 3d axis on the current fig
ax1 = fig.add_subplot(1,1,1, projection='3d')

# print(dir(ax1))
# Trisurface 3-D Mesh Plot
surf1 = ax1.plot_trisurf(x,
                         y,
                         z,
                         triangles=t,
                         linewidth=0.1,
                         antialiased=False,
                         color=default_color,
                         # cmap=newcmp,
                         shade=True,
                         # Z=voltage_new
                         )
print(dir(surf1))

# Trimesh Plot - CENTRE
# ax1.spines.set_position('centre')
# ax1.extent=(20,80,20,80)
# surf1.set_offset_position('data')

# Trimesh - Set default azim and elev
# ax1.view_init(elev=18,azim=37)


# Controls to Zoom Trimesh Plot
# ax1.xaxis.zoom(2.2)
# ax1.yaxis.zoom(2.2)
# ax1.zaxis.zoom(2.2)

ax1.set_title('OpenEP TriRep Anatomy Data')
plt.axis('on')
plt.show()


# DRAW THE FREE-BOUNDARY























# x_new = np.reshape(a=x,newshape=(len(x),1))
# y_new = np.reshape(a=y,newshape=(len(y),1))
# z_new = np.reshape(a=z,newshape=(len(z),1))
# #
# # print('x-newshape\n',x_new.shape)
# # print('x-newshape\n',y_new.shape)
# # print('x-newshape\n',z_new.shape)
#
#
#
# xx = np.arange(x_new.min(),x_new.max(),0.01)
# yy = np.arange(y_new.min(),y_new.max(),0.01)
#
# X,Y = np.meshgrid(xx,yy)
# interp = linearInt(list(zip(x, y)), z)
# Z = interp(X, Y)
#
# print('Z-shape\n',Z.shape)
# ax.scatter(X, Y, Z)
# # ax.plot_trisurf(x,y,z)
# # plt.colorbar()
# # plt.axis("equal")
# plt.show()
#
# # zz = np.arange(z_new.min(),z_new.max(),0.1)
#
# # xx,yy,zz = np.meshgrid(xx,yy,zz,sparse=False)
#
#
# # surf = ax.plot_surface(xx, yy, zz, cmap=mp.cm.coolwarm,
# #                        linewidth=0, antialiased=False)
#
# # plt.show()


def drawsurface(nodes,indices,surface_color,plot_title):
    # DRAW THE SURFACE
    # Create a new fig window
    fig = plt.figure(num=1,
                     figsize=[10,8],
                     facecolor='#F0F0F0',
                     # constrained_layout=True,
                     frameon=True,
                     tight_layout=False,
                     dpi=500)



    # get the current 3d axis on the current fig
    ax1 = fig.add_subplot(1,1,1, projection='3d')

    # print(dir(ax1))
    # Trisurface 3-D Mesh Plot
    surf1 = ax1.plot_trisurf(nodes[:,0],
                             nodes[:,1],
                             nodes[:,2],
                             triangles=indices,
                             linewidth=0.1,
                             antialiased=False,
                             color=surface_color,
                             # cmap=newcmp,
                             shade=True
                             )

    # Trimesh Plot - CENTRE
    # ax1.spines.set_position('centre')
    # ax1.extent=(20,80,20,80)
    surf1.set_offset_position('data')

    # Trimesh - Set default azim and elev
    ax1.view_init(elev=18,azim=37)


    # Controls to Zoom Trimesh Plot
    ax1.xaxis.zoom(2.2)
    ax1.yaxis.zoom(2.2)
    ax1.zaxis.zoom(2.2)

    ax1.set_title('OpenEP TriRep Anatomy Data')
    plt.axis('off')
    plt.show()

    return surf1


def colorShell(hsurf,pts,data,t,*args):
    # COLORSHELL shades the surface hSurf with data at pts
#     Where:
#    hSurf - see plotVelocityGeometry
#    pts - the data position coordinates. Size mpts x ndim
#    data - the data at each location of pts. Size mpts x ndata
#    t - threshold distance from pts greater than which shell will not be
#    coloured

#   dataField - is the data which has been drawn
#   hColBar - is a handle to the color bar

    pass
