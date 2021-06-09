from openep import io as openep_io
from openep import case as openep_case
from openep import mesh_routines as openep_mesh
from openep import case_routines as case_routines

import numpy as np
from scipy.interpolate import LinearNDInterpolator as linterp
from scipy.interpolate import NearestNDInterpolator as nearest
from sklearn.neighbors import NearestNeighbors
from matplotlib.cm import jet, rainbow, jet_r, seismic


# # Usecase 039 - Creating a voltage map from electroanatomic mapping data

filename = '/home/jra21/work/source/repos/openep_py/examples/data/new_dataset_1.mat'
distanceThresh = 10


ep_case = openep_io.load_case(filename)
ep_case_mesh = ep_case.create_mesh()
# simplify the mesh
# ep_case_mesh.merge_vertices()

# Load EP CASE - MESH Points (trirep.X)
# pts = ep_case_mesh.vertices
pts = ep_case.nodes



# # Load EGMSurfx and EGM voltage values
coords = ep_case.electric['egmSurfX'].T
data = ep_case.electric['egm'].T


iVp = case_routines.getMappingPointsWithinWoI(ep_case)
# macthing the shape of ivp with data
iVp_data = np.repeat(iVp, repeats=data.shape[1], axis=1)
# macthing the shape of ivp with coords
iVp_coords = np.repeat(iVp, repeats=coords.shape[1],axis=1)

data[~iVp_data] = np.nan
coords[~iVp_coords] = np.nan


# For each mapping point, n, find the voltage amplitude
max_volt = np.amax(a=data,axis=1).reshape(len(data),1)
min_volt = np.amin(a=data,axis=1).reshape(len(data),1)


amplitude_volt = np.subtract(max_volt,min_volt)
# Remove any data with Nans
for indx in range(amplitude_volt.shape[1]):
    tempData = amplitude_volt[:,indx]
    tempCoords = coords
    iNaN = np.isnan(tempData)
    tempData=tempData[~iNaN]
    tempCoords=tempCoords[~iNaN]


    interp = case_routines.openEpDataInterpolator(method='scatteredinterpolant',distanceThreshold=distanceThresh,rbfConstant=1)
    VertexVoltageData = interp.interpolate(x0=tempCoords,d0=tempData,x1=pts)

    print('pts-shape',pts.shape)
    # c = case_routines.LocalSmoothing(x0=tempCoords,x1=pts,smoothingLength=10)
    c = case_routines.rbf(x0=tempCoords,d0=tempData,x1=pts,rbfConstant=1)
    print(c)

ep_case.fields['d1'] = VertexVoltageData
    
# DRAW Map
openep_mesh.compute_field(mesh=ep_case_mesh,fieldname='d1',minval=0,maxval=1,color_map=jet_r)
# ep_case_mesh.show()

