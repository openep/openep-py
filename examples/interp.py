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
# class LinearNDInterpolatorExt(object):
#     def __init__(self, points, values):
#         self.funcinterp = linterp(points, values)
#         self.funcnearest = nearest(points, values)

#     def __call__(self, *args):
#         z = self.funcinterp(*args)
#         chk = np.isnan(z)
#         if chk.any():
#             return np.where(chk, self.funcnearest(*args), z)
#         else:
#             return z



# class openEpDataInterpolator():

#     method = 'scatteredinterpolant'
#     def __init__(self,method):
#         pass
#         # self.method = method


#     def interpolate(self,x0,d0,x1,*args):
#         self.x0 = x0
#         self.d0 = d0
#         self.x1 = x1

#         F = LinearNDInterpolatorExt(self.x0,self.d0,*args)
#         d1 = np.array(F(x1[:,0],x1[:,1],x1[:,2])).reshape(x1.shape[0],1)

#         # Workout if there are any points on the surface that are < 0
#         d1[d1<0] = 0

#         #  work out which points on the surface are too far away from real data
#         # Remove any interpolated values which are outwith the fill threshold 

#         neigh = NearestNeighbors(n_neighbors=1, 
#                                 radius=1, 
#                                 algorithm='auto', 
#                                 leaf_size=30,
#                                 metric='minkowski',
#                                 p=2)

#         neigh.fit(tempCoords)
#         id = neigh.kneighbors(pts,return_distance=False)
#         cPts = tempCoords[id]
#         cPts = np.array(cPts[:,0])
        
#         d = case_routines.distBetweenPoints(cPts,x1)

#         thresholdDistance = np.zeros(shape=d.shape,dtype=np.bool)
#         thresholdDistance[d>distanceThresh] = 1
#         d1[thresholdDistance] = 'nan'
#         d1 = d1.flatten()

#         return d1


        





filename = '/home/jra21/work/source/repos/openep_py/examples/data/new_dataset_2.mat'
interMethod = 'natural'
exterMethod = 'nearest'
distanceThresh = 10


ep_case = openep_io.load_case(filename)
ep_case_mesh = ep_case.create_mesh()

# Load EP CASE - MESH Points (trirep.X)
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

    interp = case_routines.openEpDataInterpolator(method='scatteredinterpolant',distanceThreshold=distanceThresh)
    VertexVoltageData = interp.interpolate(x0=tempCoords,d0=tempData,x1=pts)

ep_case.fields['d1'] = VertexVoltageData
    
# DRAW Map
openep_mesh.compute_field(mesh=ep_case_mesh,fieldname='d1',minval=0,maxval=1,color_map=jet_r)
ep_case_mesh.show()

