
from openep import io as openep_io
from openep import case as openep_case
from openep import mesh_routines as openep_mesh

import numpy as np

# Usecase 039 - Creating a voltage map from electroanatomic mapping data

#Usecase 034 - GETTING MAPPING POINT WITHIN WOI

def getMappingPointsWithinWoI(mesh_case):
    '''
    GETMAPPINGPOINTSWITHINWOI Returns the indices of the mapping points with
    annotated local activation time within the window of interest

    Args:
        mesh_case: Case

    Returns:
        logical array list of valid points; indexes into mesh_case.electric
    '''
    
    # reference annotation
    referenceAnnot = mesh_case.electric['annotations/referenceAnnot'].T

    # woi
    woi = mesh_case.electric['annotations/woi'].T

    # Extract the region of the electrogram of interest for each mapping point, n
    woi = woi + referenceAnnot

    # MapAnnot
    mapAnnot = mesh_case.electric['annotations/mapAnnot'].T

    iPoint = np.ones(shape=mapAnnot.shape, dtype=np.bool)

    # Remove the data points that are outside the region of interest
    for indx in range(mapAnnot.size):
        if (mapAnnot[indx]<woi[indx,0]) or (mapAnnot[indx]>woi[indx,1]):
            iPoint[indx] = False
        else:
            iPoint[indx] = True

    return iPoint




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



filename = '/home/jra21/work/source/repos/openep_py/tests/data/new_dataset_1.mat'
interMethod = 'natural'
exterMethod = 'nearest'
distanceThresh = 10


ep_case = openep_io.load_case(filename)
ep_case_mesh = ep_case.create_mesh()



# # Load EGMSurfx and EGM voltage values
coords = ep_case.electric['egmSurfX'].T
data = ep_case.electric['egm'].T
print('coords',coords.shape)
print('data',data.shape)

iVp = getMappingPointsWithinWoI(ep_case)
# macthing the shape of ivp with data
iVp_data = np.repeat(iVp, repeats=data.shape[1], axis=1)
# macthing the shape of ivp with coords
iVp_coords = np.repeat(iVp, repeats=coords.shape[1],axis=1)
# 
data[~iVp_data] = np.nan
coords[~iVp_coords] = np.nan


# For each mapping point, n, find the voltage amplitude
max_volt = np.amax(a=data,axis=1).reshape(len(data),1)
min_volt = np.amin(a=data,axis=1).reshape(len(data),1)

# print(max_volt.shape)
# print(min_volt.shape)

amplitude_volt = np.subtract(max_volt,min_volt)
# print(amplitude_volt.shape[1])

# Remove any data with Nans
for indx in range(amplitude_volt.shape[1]):
    tempData = amplitude_volt[:,indx]
    tempCoords = coords
    iNaN = np.isnan(tempData)
    tempData=tempData[~iNaN]
    tempCoords=tempCoords[~iNaN]


# Interpolation

















# Interpolate the voltage values across the shell 








# for volt in egm:
#     max_volt.append(np.amax(volt))
#     min_volt.append(np.amin(volt))

# max_volt = np.array(max_volt)
# min_volt = np.array(min_volt)

# #
# amplitude_volt = max_volt - min_volt




# # GENERATEINTERPDATA removes any NaN values in data (and their
# # corresponding location(s) in coords) before calling scatteredInterpolant
# # with the interpolation/extrapolation methods specified. Any values greater
# # than distancethresh are removed.


# # For each mapping point, n, find the voltage amplitude
# #  - Apply the user defined function (for example max(egm) – min(egm))
# print(iPoint)

# iPoint = np.where(iPoint==False,'nan',iPoint)
# data[iPoint] = np.where
# for item in iPoint:
#     print(item)


# print('iPoint-shape',np.where(iPoint=='False','nan',iPoint))






# data[iPoint:,] = np.where(iPoint=='False','nan',data[iPoint:,])

# print(data[np.invert(iPoint)].shape)
# print(data[np.where(~iPoint):])

# for item in c:
#     print(item)

# check for false iPoint and replace the data corresponding to the index of iPoint to False
