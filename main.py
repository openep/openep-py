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

# Declarations
file_path = '../openep/openep-examples-main/openep_dataset_1.mat'
data_tri_X_lst = list()
data_tri_T_lst = list()
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
# print('voltage data raw',voltage_data)
# indexRemove = ~np.isnan(voltage_data)
# print(indexRemove)

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

# print('type_of_t\n',type(t))

# freeboundary
# print(x.size)
# print(y.size)
# print(z.size)
# print(t.size)
# points = np.array([x,y,z])
# print(points.shape)
points = data_tri_X
print('x-y-z points\n',points.shape)
# tri= Delaunay(points)
print('sizeofT:\n',t.shape)


# Trimesh object
mesh = tm.Trimesh(vertices=points,faces=t)




# Identify the unique edge points
unique_edges = mesh.edges_unique
print('Unieq edges of the mesh:\n',unique_edges.shape)



# Identify the unshared edges
unshared_edges = mesh.face_adjacency_unshared
print('unshared\n',unshared_edges.shape)

# Pick out the unique edges in the un-shared edges
unique_unshared_edges = []

for item in unshared_edges:
    if item in unique_edges:
        unique_unshared_edges.append(item)
    else:
        pass


print('unique_unshared_edges:\n',np.array(unique_unshared_edges).shape)

# unshared_edges_points = mesh.vertices[unshared_edges]
# print('vertices-unshared:\n',unshared_edges_points)


# unshared_edges = tm.graph.face_adjacency_unshared(mesh)
# print('unshared_edges\n',unshared_edges)

# scene = tm.Scene([mesh,unshared_edges])
#
# scene.show()


x_values = list()
y_values = list()
z_values = list()
# bound_edges = list()

#
# fac_bound_edges = mesh.facets_boundary
#
# print('mesh edges',mesh.edges)
# print('Unique edges of the mesh:\n',mesh.edges_unique)
# print('facets boundary - edges\n',type(fac_bound_edges))

# for item in fac_bound_edges:
#     for j in range(len(item)):
#         bound_edges.append(item[j])
#         # print('elements\n',item[j])
#     # for j in item:
#     #     print('edges-boundary:\n',item[j])
#
# print('boundary-edges:\n',np.array(bound_edges).shape)
# edgeverts=mesh.vertices[bound_edges]
# print('edgeverts:\n',edgeverts[0][0])
#
 # [[1138 1311]
 # [1311 1531]
 # [1138 1768]
 # [1531 1768]]

# uniques={}
#
# for a,b in mesh.edges:
#     if (a,b) not in uniques and (b,a) not in uniques:
#         uniques[(a,b)]=(a,b)
#     else:
#         if (a,b) in uniques:
#             del uniques[(a,b)]
#         else:
#             del uniques[(b,a)]
#
# print('unique edges\n',uniques)
# edgeverts=mesh.vertices[np.array(list(uniques))]
# print('edgeverts:\n',edgeverts[0])


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
xx = np.array(x_values).reshape(len(x_values),1)
yy = np.array(y_values).reshape(len(y_values),1)
zz = np.array(z_values).reshape(len(z_values),1)
print('xx-values\n',xx.shape)
print('yy-values\n',yy.shape)
print('zz-values\n',zz.shape)


# # Join the unique edges together
unique_unshared_edge_points = np.concatenate((xx,yy,zz),axis=1)
print('x_values:\n',unique_unshared_edge_points.shape)


# ax1.scatter3D(edgeverts[:0],edgeverts[:1],edgeverts[:2], cmap='Greens')

# xx = list(np.array(x_values).reshape(len(x_values),1))
# yy = list(np.array(y_values).reshape(len(y_values),1))
# zz = list(np.array(z_values).reshape(len(z_values),1))

# points1 = np.array([xx,yy,zz])
# print('points:\n',points)
# print('points1:\n',points1)
# new_points = np.array(x,y,z)
# print('new_points:\n',new_points)
# print('x-values:\n',xx)

# mesh1 = tm.Trimesh(vertices=unique_unshared_edge_points,faces=t)
# mesh1.show()
# mesh.show()
# mesh.as_open3d()



# Return the vertex index of the two vertices not in the shared edge between two adjacent face
# edges = tm.graph.face_adjacency_unshared(mesh)
# print('edges:\n',edges.shape)
# print('points of index 4281:\n',points[4281])

#
# f_boundary = mesh.facets_boundary()
# print('facets_ Boundary:\n',f_boundary)



# x_values = list()
# y_values = list()
# z_values = list()
# loc1 = list()
# loc2 = list()

# for item in edges:
#     loc1.append(item[0])
#     loc2.append(item[1])
    # x(3*item+1) = loc1(0)
    # y(3*item+1) = loc1(1)
    # x(3*item+1) = loc1(2)

# x_values = np.array([x[loc1], x[loc2]])
# y_values = [y[loc1], y[loc2]]
# z_values = [z[loc1], z[loc2]]



    # print('index points:',item[1])
    # x_values.append([x[item[0]],x[item[1]]])
    # y_values.append([y[item[0]],y[item[1]]])
    # z_values.append([z[item[0]],z[item[1]]])

# print('loc1:\n', loc1)
# print('x-values:\n',x_values.shape)
# print('y-values:\n',np.array(y_values))
# print('z-values:\n',np.array(z_values))






# print('Delaunay_Tri\n',tri.simplices.shape)

# load Triangulation matrix
# Get the simplices
# Detect free boundary - where it only has one simplex
# plot the results

# faces = points[t]
# print('faces of triangle\n',faces.shape)



# Convex Hull
# hull = ConvexHull(points)
# print('size_hull',hull.simplices.shape)
#
# for simplex in hull.simplices:
#     print('Simplex points\n',points[simplex,0],points[simplex,1],points[simplex,2])
#     plt.plot(points[simplex,0], points[simplex,1], points[simplex,2],'k-')





# print('simplices shape:\n',tri.simplices.shape)

# ax1.plot(x_values,y_values,z_values)

# for item in t:
#     # print('triangles\n',points[item].size)
#     if pointsset_markeredgecolor[item].size == 2:
#         print('hurray')
#     else:
#         pass
    # for i in len(points[item])



    # if points[item][]==0:
    #     print('found you', points[item])
    # else:
    #     pass


# for i in range(tri.simplices.size):
#     # if points[tri.simplices[i,:]
#     # print('points for index ',i,':\n',tri.simplices[i,:])
#     # print('index', i,':\n',tri.simplices[i,:4])
#     tri_supplices = np.array(tri.simplices[i,:4])
#     if tri_supplices[0] == 0 or tri_supplices[1] == 0 or tri_supplices[2] == 0 or tri_supplices[3] == 0:
#         print('facet is on the free boundary:\n', tri_supplices)
#     else:
#         pass

    # print(tri_supplices[0])
    # if (tri_supplices == 0):
    #     print('facet is on the free boundary:\n', tri_supplices)
    # else:
    #     pass

    # if tri_supplices[i,:4] == 0:
    #     print('facet is on the free boundary:\n', tri_supplices)
    # else:
    #     pass
    # # if ( (tri.simplices[i,:0] == 0) or (tri.simplices[i,:1] == 0) or (tri.simplices[i,:2] == 0) or (tri.simplices[i,:3] == 0) ):
    # if tri.simplices[i,:] == 0:
    #     print('facet is on the free boundary:\n', tri.simplices[i,:])
    # else:
    #     pass


# print(dir(tri))




# hull = ConvexHull(t,incremental=True)
# print(type(hull))

# plt.plot(t[:,0], t[:,1], t[:,2], 'o')
# # hull.simplices gives the indices of the points forming the simplices in the triangulation
# for simplex in tri.simplices:
#     # plt.plot(tri[simplex, 0], tri[simplex, 1], tri[simplex,2],'r-')
#     plt.plot(tri[simplex, 0], tri[simplex, 1], tri[simplex,2],'r-')
# print('data_X values',data_tri_X.shape)




# - Tri-Surface Plot with a rainbow color map
# surf1 = ax1.plot_trisurf(xx,
#                          yy,
#                          zz,
#                          linewidth=0.5,
#                          antialiased=True,
#                          color='black')

surf1 = ax1.plot_trisurf(x,
                         y,
                         z,
                         triangles=t,
                         linewidth=0.5,
                         antialiased=True,
                         cmap=newcmp)


# print('dir:\n',dir(surf1))
# print('typeofsurf',type(surf1))
surf1.set_edgecolor([0,0,0,1])
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
