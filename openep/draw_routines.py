from openep import mesh_routines as openep_mesh

import pyvista as pv
import numpy as np
from matplotlib.cm import jet, rainbow, jet_r, seismic
import trimesh as tm


plot = False

def lineLength(h):
    
    '''
    Calculates the Length of a line
    
    Args:
        h: mx3 matrix of cartesian co-ordinates representing the line
        the X,Y, Z data are received directly in a matrix of the form:
         [ x_1  y_1  z_1 ]
         [ x_2  y_2  z_2 ]
         [ ...  ...  ... ]
         [ x_n  y_n  z_n ]
         
    Returns:
        l - length of the line as a float
    '''
#   remove any Nan values
    data = h[~np.isnan(h)].reshape(h.shape[0],h.shape[1])
        
    dx = np.diff(data[:,0])
    dy = np.diff(data[:,1])
    dz = np.diff(data[:,2])
    
    l = 0
    for i in range(len(dx)):
        l = l+ np.sqrt(np.square(dx[i]) + np.square(dy[i]) + np.square(dz[i]))
        
    return l



def getAnatomicalStructures(mesh_case, *args):
    
    '''
    GETANATOMICALSTRUCTURES Returns the free boundaries (anatomical 
    structures) described in ep_case
    
    Args:
        mesh_case : Case
        
    Returns:
        FreeboundaryPoints : m x 3 matrix of the freeboundary coordinate points of each anatomical structure
        l                  : array of lengths [perimeters] of each anatomical structure
        a                  : array of areas of each anatomical structure
        tr                 : array of trimesh objects of each anatomical structure
    
    
    GETANATOMICALSTRUCTURES identifies all the anatomical structures of a 
    given data set. Anatomical structures are boundary regions that have been 
    added to an anatomical model in the clinical mapping system. For example, 
    with respect of left atrial ablation, anatomical structures may represent 
    the pulmonary vein ostia, mitral valve annulus or left atrial appendage 
    ostium.
        
    '''


    x_values = []
    y_values = []
    z_values = []

    tr = []
    a = []
    dist = []
    l = []
    
    pts = mesh_case.nodes
    face = mesh_case.indices.astype(int)
        
    tm_mesh = tm.Trimesh(vertices=pts,faces=face)
    freeboundary = tm_mesh.outline().to_dict()
    mesh_outline_entities = tm_mesh.outline().entities
    freeboundary_vertex_nodes = []
    freeboundary_points = []


    for i in range(len(mesh_outline_entities)):
        freeboundary_vertex_nodes.append(freeboundary["entities"][i]['points'])


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

    freeboundary_points = np.array(freeboundary_points,dtype=object)
    
    

    for i in range(len(freeboundary_points)):
        coords = freeboundary_points[i]
        centre = np.round(np.mean(coords,0),3)
        centre = centre.reshape(1,len(centre))

        # Create a Trirep of the boundary
        X = np.vstack((centre,coords))
        numpts = X.shape[0]
        A = np.zeros(((numpts-1),1),dtype=np.int64)
        B = np.array(range(1,numpts)).T
        B = B.reshape(len(B),1)
        C = np.array(range(2,numpts)).T
        C = C.reshape(len(C),1)
        C = np.vstack((C,1))
        TRI = np.concatenate((A,B,C),axis=1)

        for j in range(len(coords)-1):
            dist.append(np.linalg.norm(coords[j+1]-coords[j]))
            

        tr.append(tm.Trimesh(vertices=X,faces=TRI))
        lineLen = lineLength(coords)
        a.append(round(tr[i].area,4))
        l.append(round(lineLen,4))
        print('Perimeter is :',l[i],'| Area is:',a[i])


    return {'FreeboundaryPoints':freeboundary_points,'Lengths':l, 'Area':a, 'tr':tr}
       



def DrawMap(ep_case,freeboundary_color,freeboundary_width,minval,maxval,volt_below_color, volt_above_color, nan_color, plot,**kwargs):
    '''
    DrawMap - 
    '''

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

    p = pv.Plotter()

    if plot:

        # Plot OpenEp mesh
        sargs = dict(interactive=True, 
                          n_labels=2,
                          label_font_size=18,
                          below_label='  ',
                          above_label='  ')

        freebound = getAnatomicalStructures(ep_case)
        p.add_mesh(mesh,
               scalar_bar_args=sargs,
               show_edges=False,
               smooth_shading=True,
               scalars=volt,
               nan_color=nan_color,
               clim=[minval,maxval],
               cmap=jet_r,
               below_color=volt_below_color,
               above_color=volt_above_color)

    
        # Plot free Boundary - Lines
        for indx in range(len(freebound['FreeboundaryPoints'])):
            p.add_lines((freebound['FreeboundaryPoints'][indx]),
                        color=freeboundary_color,
                        width=freeboundary_width)
        p.show()

    return {'hsurf':p,
            'pyvista-mesh':mesh,
            'volt':volt,
            'nan_color':nan_color,
            'minval':minval,
            'maxval':maxval,
            'cmap':jet_r,
            'volt_below_color':volt_below_color,
            'volt_above_color':volt_above_color}


    