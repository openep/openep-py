from numpy.lib.function_base import _ARGUMENT_LIST
from openep import mesh_routines as openep_mesh


import pyvista as pv
import numpy as np
from matplotlib.cm import jet, rainbow, jet_r, seismic
import trimesh as tm

# GLOBAL VARIABLES
plot = False
volt = 0

def lineLength(h):   
    """
    Calculates the Length of a line
    
    Args:
        h (float): mx3 array of cartesian co-ordinates representing the line
         
    Returns:
        float: mx1 array of length of the line    
    """
    #   remove any Nan values
    data = h[~np.isnan(h)].reshape(h.shape[0],h.shape[1])
        
    dx = np.diff(data[:,0])
    dy = np.diff(data[:,1])
    dz = np.diff(data[:,2])
    
    l = 0
    for i in range(len(dx)):
        l = np.round(l+ np.sqrt(np.square(dx[i]) + np.square(dy[i]) + np.square(dz[i])),4)
        
    return l

def create_pymesh(mesh_case):

    """
    Create a pymesh object
    Args:
        mesh_case (obj): Case object from a given openep file

    Returns:
        Pyvista PolyData (obj) - triangulated surface object from Numpy arrays of the vertices and faces
    """

    pts = mesh_case.nodes
    tri = mesh_case.indices.astype(int)
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

    return mesh

def get_freeboundaries(tri):

    """
    Args:
        tri (obj) - Trimesh Object containing a triangular 3D mesh

    Returns:
           int: freeboundaries - Array of free boundary indices
           float: fb_vertex_nodes - mx3 array of cartesian co-ordinates of the freeboundaries
    """

    freeboundary_vertex_nodes = []
    freeboundaries = []

    tri_outline = tri.outline().to_dict()
    mesh_outline_entities = tri.outline().entities
    no_of_freeboundaries = len(mesh_outline_entities)


     # Find all the freeboundary facets
    for i in range(no_of_freeboundaries):
        freeboundary_vertex_nodes.append(tri_outline["entities"][i]['points'])

        # creating an array of fb values with the index of the points
        for j in range(len(freeboundary_vertex_nodes[i])-1):
            freeboundaries.append([freeboundary_vertex_nodes[i][j], freeboundary_vertex_nodes[i][j+1]]) 

    freeboundaries = np.array(freeboundaries).astype(np.int64)
    
    return {'freeboundary':freeboundaries,
            'fb_vertex_nodes':freeboundary_vertex_nodes}



def get_freeboundary_points(tri,fb):

    trimesh_points = tri.vertices
    fb = fb[:,0]
    coords = np.array(list(trimesh_points[fb])).astype(np.float64)
    
    return coords


def free_boundary(tri):
    
    l = []

    f_fall = get_freeboundaries(tri)
    fb = f_fall['freeboundary']

    test_set=np.array(list(np.zeros((fb.shape[0],fb.shape[1]))))
    test_set[:] = np.nan
    test_set[:,0] = fb[:,0]

    #exctracting the contents of 2nd columm fb array into the testset
    second_col_fb = list(map(lambda x:x, fb[:,1][0:len(fb)-1]))
    second_col_fb.insert(0,fb[:,1][-1])
    second_col_fb = np.array(second_col_fb)
    
    # Assigning the values of the second column to the 2nd column of testset
    test_set[:,1] = second_col_fb

    # isstat
    i_start = np.where(np.diff(a=test_set,n=1,axis=1).astype(np.int64))[0]
    
    ff = []

    if not list(i_start):
        ff.append(np.array(fb))
        coords = get_freeboundary_points(tri,ff[0])
        l.append(lineLength(coords))
    else:
        for i in range(i_start.shape[0]):
            if i<(i_start.shape[0]-1):
                np.array(ff.append(fb[i_start[i]:i_start[i+1]]))
                coords = get_freeboundary_points(tri,ff[i])
                l.append(lineLength(coords))
            else:
                ff.append(fb[i_start[i]:])
                coords = get_freeboundary_points(tri,ff[i])
                l.append(lineLength(coords))



    return {'connected_freeboundary':ff,
            'length':l}





def drawFreeBoundaries(mesh_case,
                       fb_points,
                       fb_col,
                       fb_width,
                       mesh_surf_color,
                       opacity,
                       smooth_shading,
                       use_transparency,
                       lighting,**kwargs):
    '''
    Draws the boundaries at the edge of a TriSurface mesh with each freeboundary rendered with the colors mentioned in the fb_col 
    Args:
        mesh_case         : openep Case object
        fb_points         : m x 3 coordinate point arrays
        fb_col            : list of RGB colors For eg: fb_col = ['blue','yellow','green','red','orange','brown','magenta']
        fb_width          : Float
        mesh_surf_color   : RGB values between 0 and 1
        opacity           : Float - any value in the range 0-1 
        smooth_shading    : Boolean - True | False  
        use_transparency  : Boolean - True | False
        lighting          : Boolean - True | False
    
    Returns:
        p                 : pyvista plotter handle

    '''

    # Create a pymesh of the openep case
    py_mesh = create_pymesh(mesh_case)
    
    #     Freeboundary Color
    fb_col = ['blue','yellow','green','red','orange','brown','magenta']

    p = pv.Plotter()
    
    p.add_mesh(py_mesh, 
               color=mesh_surf_color,
               opacity=opacity,
               smooth_shading=smooth_shading,
               use_transparency=use_transparency,
               lighting=lighting)
    
    for indx in range(len(fb_points)):
            p.add_lines(fb_points[indx],
                        color=fb_col[indx],
                        width=fb_width)    
    return p

def getAnatomicalStructures(mesh_case, plot, **kwargs):
    
    '''
    GETANATOMICALSTRUCTURES Returns the free boundaries (anatomical 
    structures) described in ep_case and plots them if plot=True with drawFreeBoundaries()
    
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

    a = []
    l = []

    tr = []
    
    x_values = []
    y_values = []
    z_values = []
    
    freeboundary_points = []
    freeboundary_vertex_nodes = []
    
    
    pts = mesh_case.nodes
    face = mesh_case.indices.astype(int)
        
    tm_mesh = tm.Trimesh(vertices=pts,
                        faces=face, 
                        process=False)

       
    freeboundary = tm_mesh.outline().to_dict()
    mesh_outline_entities = tm_mesh.outline().entities

    for i in range(len(mesh_outline_entities)):
        freeboundary_vertex_nodes.append(freeboundary["entities"][i]['points'])
                
    
    # mesh Vertices - finding the coordinates/points of the freeboundaries
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
            

        tr.append(tm.Trimesh(vertices=X,faces=TRI))
        lineLen = lineLength(coords)
        a.append(round(tr[i].area,4))
        l.append(lineLen)
        print('Perimeter is :',l[i],'| Area is:',a[i])
        
    if plot:
        col = ['blue','yellow','green','red','orange','brown','magenta']
        fb_width = 3
        mesh_surf = [0.5,0.5,0.5]
        p = drawFreeBoundaries(mesh_case=mesh_case,
                               fb_points=freeboundary_points,
                               fb_col = col,
                               fb_width = fb_width,
                               mesh_surf_color = mesh_surf,
                               opacity = 0.2,
                               smooth_shading= True,
                               use_transparency=False,
                               lighting=False,
                               **kwargs)
            
        p.background_color='White'
        p.show()
    
        
    return {'FreeboundaryPoints':freeboundary_points,'Lengths':l, 'Area':a, 'tr':tr}

def create_pymesh(mesh_case):

    pts = mesh_case.nodes
    tri = mesh_case.indices.astype(int)
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

    return mesh


    

def DrawMap(ep_case,volt,freeboundary_color,cmap,freeboundary_width,minval,maxval,volt_below_color, volt_above_color, nan_color, plot,**kwargs):
    '''
    DrawMap - plots an OpenEp Map
    Args:
        ep_case,
        volt (str | nx1 array) - 'bip' | interpolated voltagae values 
        freeboundary_color,
        cmap,
        freeboundary_width,
        minval,
        maxval,
        volt_below_color, 
        volt_above_color, 
        nan_color, 
        plot

    Returns:

    '''
    # Create a PyMesh
    py_mesh = create_pymesh(ep_case)
    
    if volt == 'bip':
        volt = ep_case.fields['bip']
    else:
        volt = volt

    # pyvista plotter
    p = pv.Plotter()

    if plot:

        # Plot OpenEp mesh
        sargs = dict(interactive=True, 
                          n_labels=2,
                          label_font_size=18,
                          below_label='  ',
                          above_label='  ')

        freebound = getAnatomicalStructures(ep_case,plot=False)
        p.add_mesh(py_mesh,
               scalar_bar_args=sargs,
               show_edges=False,
               smooth_shading=True,
               scalars=volt,
               nan_color=nan_color,
               clim=[minval,maxval],
               cmap=cmap,
               below_color=volt_below_color,
               above_color=volt_above_color)

    
        # Plot free Boundary - Lines
        for indx in range(len(freebound['FreeboundaryPoints'])):
            p.add_lines((freebound['FreeboundaryPoints'][indx]),
                        color=freeboundary_color,
                        width=freeboundary_width)
        p.show()

    return {'hsurf':p,
            'pyvista-mesh':py_mesh,
            'volt':volt,
            'nan_color':nan_color,
            'minval':minval,
            'maxval':maxval,
            'cmap':cmap,
            'volt_below_color':volt_below_color,
            'volt_above_color':volt_above_color}

    