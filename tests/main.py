# Code written by Jerin Rajan on 18th Mar 2021
# Functions to support OpenEp visualisation

from enum import Enum, auto
import numpy as np
import trimesh as tm
import scipy.io as sio
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, LightSource
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import make_interp_spline, BSpline, griddata, interp1d, NearestNDInterpolator
from matplotlib.cm import jet_r
import vtk

# DEFINITIONS
file_path = '../openep_py/tests/data/openep_dataset_1.mat'
data_tri_X_lst = []
data_tri_T_lst = []
x_values = []
y_values = []
z_values = []
# voltage threshold to split the colorbar by solid and JET
volt_threshold = 2
below_color=(0, 0, 0, 255)
above_color=(255, 0, 255, 255)
nan_color=(50/255, 50/255, 50/255, 255/255)
random_color = (6, 202, 255, 0.9)
color_map=jet_r
minval=0,
maxval=1,

color1 = []


class VisualisationBackend(Enum):
    Matplotlib = auto()
    VTK = auto()


# Switch these for matplotlib/VTK comparison
visualisation_backend = VisualisationBackend.VTK
# visualisation_backend = VisualisationBackend.Matplotlib


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

x_new = np.reshape(a=x,newshape=(len(x),1))
y_new = np.reshape(a=y,newshape=(len(y),1))
z_new = np.reshape(a=z,newshape=(len(z),1))
#
print('x-newshape\n',x_new.shape)
#
# xx = np.zeros(shape=(len(x),1),dtype=np.float64) + x_new
# yy = np.zeros(shape=(len(y),1),dtype=np.float64) + y_new
# zz = np.zeros(shape=(len(z),1),dtype=np.float64) + z_new
#
# print('xx\n',xx.shape)
#
# print(xx.max(),xx.min())
# print(yy.max(),yy.min())
# print(zz.max(),zz.min())
#
# # AXIS vectorspoints, values, (grid_x, grid_y
# # xx = np.linspace(int(x.min()),int(x.max()),len(x),dtype='uint8')
# # yy = np.linspace(int(y.min()),int(y.max()),len(y),dtype='uint8')
# # zz = np.linspace(int(z.min()),int(z.max()),len(z),dtype='uint8')
# #
# #
# # print('xx-shape\n',xx)
# # print('yy-shape\n',yy)
# # print('zz-shape\n',zz)
#
# # Meshgrid
# xxx,yyy = np.meshgrid(xx,yy)
# print('xxx\n',xxx)
# print('yyy\n',yyy)

# x=x_new,y=y_new,z=z_new,kind='near')








# Triangulation/faces - .T
data_tri_T = main_file['userdata']['surface_triRep_Triangulation'][0][0]
# Resolving indexing issue - matlab to python
t = data_tri_T-1
print('t\n',t.shape)

# Voltage Data
data_act_bip = main_file['userdata']['surface'][0][0]['act_bip'][0][0]
voltage_data = data_act_bip[:,1]
print('voltag_bip\n',voltage_data.shape)

voltage_new = np.reshape(a=voltage_data,newshape=(len(voltage_data),1))
print('voltage-new\n',voltage_new.shape)

# Checking for NaN values and replacing with 0 in voltage data
# voltage_new = np.asarray(list(map(lambda x:0 if np.isnan(x) else x, voltage_data)))

# ColorShell
# coloring the shell by the voltage value rather than z-axis height
min_volt = min(voltage_data)
max_volt = max(voltage_data)
print('min_volt',min_volt)

# normalisation - normalising the voltage - values
norm = mp.colors.Normalize(vmin=min_volt, vmax=max_volt)

# Plot the Mesh
# Create a new fig window
fig = plt.figure(constrained_layout=True,
                 frameon=False,
                 tight_layout=False,
                 figsize=(40,20),
                 dpi=150)
# first plot
# get the current 3d axis on the current fig
ax1 = fig.add_subplot(1,1,1, projection='3d')

# defining a custom color scalar field
# if only nodal values are known, then the signal need to be averaged by the triangles

# color2 = voltage_data[t]
# print('color-2',color2)

field = voltage_data
new_field = np.zeros((field.shape[0], 4), np.int32)


# creatig the color array
for idx in range(len(field)):
    val = field[idx]
    # print(val)

    if np.isnan(val):
        color1.append(nan_color)
    elif(val<minval):
        color1.append(below_color)
    elif(val>maxval):
        color1.append(above_color)
    else:
        color1.append(random_color)
        # num = val-minval
        # den = maxval-minval
        # valscaled = num / den



color1 = np.array(color1)
print('color1\n',color1.shape)






















print('voltage_new',voltage_new.shape)
print('t',t)
print('voltage_new[t]',voltage_new[t].shape)

color3d = voltage_new[t]
print('color3d',color3d.shape)


color = np.mean(voltage_new[t], axis=1)
print('color',color.shape)

# color_r = np.linspace(color.min(),color.max(),len(color))
# color_interp = interp1d(color_r)
# print('color-arry\n',color_r)

# Voltage Threshold - to split the colorbar axis
# determining the size of the voltage values
diff_volt = max_volt - min_volt
dist_below_threshold = volt_threshold - min_volt
dist_above_threshold = max_volt - volt_threshold

# splitting the colorbar
# color below threshold - size
size_below_threshold = round((dist_below_threshold/diff_volt)*256)
# color above threshold - size
size_above_threshold = int(256 - size_below_threshold)

# Define Colormap - Reverse JET
cm_jet = plt.cm.get_cmap('jet_r',size_below_threshold)
newcolors = cm_jet(np.linspace(0,1,256))
# Magenta - RGBA array
magenta = np.array([255/255, 0/255, 255/255, 1])
# magenta = np.array([139/255, 0/255, 139/255, 1])
# Gray - RGBA array
gray = np.array([128/255, 128/255, 128/255, 1])

# Creating new colormap columns
# Col 1 - Jet
col1 = cm_jet(np.linspace(0,1,size_below_threshold))
print('col1\n',col1)
# Col2 - Magenta
col2 = np.zeros(shape=(size_above_threshold,4))+magenta
# col3 = np.ones(size_above_threshold)
# Combining Col1 + Col2
final_col = np.concatenate((col1,col2))
print('final_col\n',final_col.shape)
# NewColormap
newcmp = ListedColormap(final_col)

# Points - vertices
points = data_tri_X
print('points\n',points)
# Trimesh object
mesh = tm.Trimesh(vertices=points,faces=t)
# mesh.show()
# mesh.visual.base.Visuals(vertex_colors=final_col)
# mesh.show()



freeboundary = mesh.outline().to_dict()
vertex_index = mesh.outline().vertex_nodes
mesh_outline_entities = mesh.outline().entities

freeboundary_vertex_nodes = []
freeboundary_points = []
for i in range(len(mesh_outline_entities)):
    freeboundary_vertex_nodes.append(freeboundary["entities"][i]['points'])

# mesh Vertices
trimesh_points = mesh.vertices

for i in freeboundary_vertex_nodes:
    x_values.append(trimesh_points[i][:,0])
    y_values.append(trimesh_points[i][:,1])
    z_values.append(trimesh_points[i][:,2])


x_values_array = []
y_values_array = []
z_values_array = []


for i in range(len(mesh_outline_entities)):
    x_values_array.append(x_values[i].reshape(len(x_values[i]),1))
    y_values_array.append(y_values[i].reshape(len(y_values[i]),1))
    z_values_array.append(z_values[i].reshape(len(z_values[i]),1))
    freeboundary_points.append(np.concatenate((x_values_array[i],
                                         y_values_array[i],
                                         z_values_array[i]),
                                         axis=1))

# Plotting the freeboundary edges of the 3-d Mesh
for item in freeboundary_points:
    freeboundary_plot = ax1.plot(item[:,0],
                                 item[:,1],
                                 item[:,2],
                                 c='black',
                                 linewidth=1,
                                 zorder=3
                                 )

# create a light source

# Z = jn(0,z)
ls = LightSource(azdeg=315,altdeg=65)
# Shade data, creating a rgb array
rgb = ls.shade(final_col, cmap=newcmp)
print('rgb\n',rgb.shape)
# plt.imshow(rgb)

if visualisation_backend == VisualisationBackend.VTK:

    def get_point_idx(triangulation_data, triangle_idx, point_idx):
        return triangulation_data[triangle_idx, point_idx]
        # return x[idx], y[idx], z[idx]

    def get_color_from_voltage(v, v_max, lookup_table, nan_color):
        if not np.isnan(v):
            idx = int(255 * v / v_max)
            rgb = lookup_table[idx][:3] # rgb between 0 and 1
        else:
            rgb = nan_color
        return [q * 255 for q in rgb]



    colours = vtk.vtkUnsignedCharArray()
    colours.SetNumberOfComponents(3)
    colours.SetName("Colors")

    points = vtk.vtkPoints()
    polys = vtk.vtkCellArray()
    vertices = vtk.vtkCellArray()
    max_voltage = np.nanmax(voltage_new)
    nan_color = (0, 0, 0)
    for i in range(0, t.shape[0]):
        polygon = vtk.vtkPolygon()
        polygon.GetPointIds().SetNumberOfIds(3)
        for j in range(3):
            idx = 3 * i + j
            point_idx = get_point_idx(t, i, j)
            # point = get_coords_of_point_on_triangle(t, i, j, x, y, z)
            id = points.InsertNextPoint(x[point_idx], y[point_idx], z[point_idx])
            polygon.GetPointIds().SetId(j, idx)
            vertices.InsertNextCell(1)
            vertices.InsertCellPoint(id)
            c = get_color_from_voltage(voltage_new[point_idx], max_voltage, final_col, nan_color)
            colours.InsertNextTuple(c)
        polys.InsertNextCell(polygon)
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(polys)
    # polydata.SetVerts(vertices)
    polydata.GetPointData().SetScalars(colours)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)
    mapper.SetColorModeToDefault()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetPointSize(5)

    # Renderer
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(.2, .3, .4)
    renderer.ResetCamera()

    # Render Window
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)

    # Interactor
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)

    # Begin Interaction
    renderWindow.Render()
    renderWindowInteractor.Start()

else:  # if visualisation_backend == VisualisationBackend.Matplotlib:
    # Trisurface 3-D Mesh Plot
    surf1 = ax1.plot_trisurf(x,
                            y,
                            z,
                            triangles=t,
                            linewidth=0.1,
                            antialiased=True,
                            # color='gray',
                            cmap=newcmp,
                            shade=True
                            )

    # passing the custom color scalar field to point the voltage values not z-values

    surf1.set_array(color)

    # surf1.set_array(color3d[:,0])
    # surf1.set_array(color3d[:,1])
    # surf1.set_array(color3d[:,2])

    # surf1.set_facecolors(rgb)
    ax1.set_title('OpenEP TriRep Anatomy Data')
    plt.axis('off')
    # Colour Pallette
    cb = plt.colorbar(mp.cm.ScalarMappable(norm=norm, cmap=newcmp),
                    ax=[ax1],
                    location='left',
                    label='Voltage (mV)')

    # Show plot

    plt.show()



    # UseCase03
