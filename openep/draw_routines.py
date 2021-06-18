from openep import mesh_routines as openep_mesh

import pyvista as pv
import numpy as np
from matplotlib.cm import jet, rainbow, jet_r, seismic
import trimesh as tm


def DrawMap(ep_case,freeboundary_color,freeboundary_width):
    pts = ep_case.nodes
    tri = ep_case.indices.astype(int)
    volt = ep_case.fields['bip']
    


    np.set_printoptions(suppress=True)
    size_tri = []

    x = tri[:,0].reshape(len(tri[:,0]),1)
    y = tri[:,1].reshape(len(tri[:,1]),1)
    z = tri[:,2].reshape(len(tri[:,2]),1)


    for item in tri:
        size_tri.append(len(item))

    size_tri_list = size_tri
    size_tri_arr = np.array(size_tri,dtype=np.int).reshape(len(size_tri),1)

    face = np.hstack(np.concatenate((size_tri_arr,tri),axis=1)).astype(int)
    mesh = pv.PolyData(pts,face)



    # GetFreeBoundary()
    tm_mesh = tm.Trimesh(vertices=pts,faces=tri)
    freeboundary = tm_mesh.outline().to_dict()
    mesh_outline_entities = tm_mesh.outline().entities

    freeboundary_vertex_nodes = []
    freeboundary_points = []


    for i in range(len(mesh_outline_entities)):
        freeboundary_vertex_nodes.append(freeboundary["entities"][i]['points'])

    x_values = []
    y_values = []
    z_values = []

    # mesh Vertices
    trimesh_points = tm_mesh.vertices

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


    p = pv.Plotter()

    minval=0
    maxval=2
    color_map = jet_r
    below_color=[149, 85, 0]
    # print('below_color',below_color)
    above_color=[255, 0, 255]
    nan_color=[180, 180, 180]

    field = ep_case.fields['bip']
    # new_field = np.zeros((field.shape[0], 3), np.int32)
    new_field = []

    for idx in range(len(field)):
        # print(idx)
        val = field[idx]
        # print(val)

        if np.isnan(val):
            # print('first')
            col = nan_color
        elif val < minval:
            # print('second')
            col = below_color
        elif val > maxval:
            # print('third')
            col = above_color
        else:
            # print('hello')
            valscaled = (val - minval) / (maxval - minval)
            col = color_map(valscaled)
            col = col[0:3]

            if isinstance(col[0], float):
                col = [int(c * 1) for c in col]
        new_field.append(col)
    
    # color_array = openep_mesh.compute_field(mesh=tm_mesh,fieldname='bip',minval=0,maxval=1,color_map=jet_r)
    # print(color_array)

    color_array = new_field
    # print(len(color_array))
    # print(color_array)


    # print(color_array[0])
    # print(color_array)

    # Plot OpenEp mesh
    p.add_mesh(mesh,show_edges=False,smooth_shading=True,scalars=volt,nan_color=nan_color,clim=[minval,maxval],cmap=jet_r,below_color=below_color,above_color=above_color)

    # Plot free Boundary - Lines
    for indx in range(len(freeboundary_points)):
        p.add_lines(freeboundary_points[indx],color=freeboundary_color,width=freeboundary_width)
        # p.show()
        # p.add_lines(freeboundary_points[indx],color='black',width=freeboundary_width)
        # p.add_lines(freeboundary_points[2],color='pink',width=freeboundary_width)
        # p.add_lines(freeboundary_points[3],color='red',width=freeboundary_width)
        # p.add_lines(freeboundary_points[4],color='green',width=freeboundary_width)
        # p.add_lines(freeboundary_points[5],color='orange',width=freeboundary_width)
        # p.add_lines(freeboundary_points[6],color='magenta',width=freeboundary_width)

    # p.show()
    
    return mesh,volt,nan_color,minval,maxval,jet_r,below_color,above_color,freeboundary_points


    def getAnatomicalStructures(ep_case, *args):
        pass
