from openep import mesh_routines as openep_mesh

import pyvista as pv
import numpy as np
from matplotlib.cm import jet, rainbow, jet_r, seismic
import trimesh as tm

def getAnatomicalStructures(ep_case, *args):
        pts = ep_case.nodes
        tri = ep_case.indices.astype(int)
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
        return freeboundary_points



def DrawMap(ep_case,freeboundary_color,freeboundary_width,minval,maxval,volt_below_color, volt_above_color, nan_color):
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
    freebound = getAnatomicalStructures(ep_case)

    p = pv.Plotter()

    # Plot OpenEp mesh
    p.add_mesh(mesh,
               show_edges=False,
               smooth_shading=True,
               scalars=volt,
               nan_color=nan_color,
               clim=[minval,maxval],
               cmap=jet_r,
               below_color=below_color,
               above_color=above_color)

    # Plot free Boundary - Lines
    for indx in range(len(freebound)):
        p.add_lines(freebound[indx],color=freeboundary_color,width=freeboundary_width)

    p.show()
    
    return mesh,volt,nan_color,minval,maxval,jet_r,below_color,above_color,freebound


    