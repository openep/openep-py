
from openep import io as openep_io
from openep import case as openep_case
from openep import mesh_routines as openep_mesh

import numpy as np

# Usecase 039 - Creating a voltage map from electroanatomic mapping data

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





# GENERATEINTERPDATA removes any NaN values in data (and their
# corresponding location(s) in coords) before calling scatteredInterpolant
# with the interpolation/extrapolation methods specified. Any values greater
# than distancethresh are removed.

#  Usecase 034 - GETTING MAPPING POINT WITHIN WOI

# reference annotationprint('egm',egm.shape)
referenceAnnot = ep_case.electric['annotations/referenceAnnot'].T

# woi
woi = ep_case.electric['annotations/woi'].T
# print('egm',woi.shape)

# Extract the region of the electrogram of interest for each mapping point, n
woi = woi + referenceAnnot
# print('region of interest\n',woi)

# MapAnnot
mapAnnot = ep_case.electric['annotations/mapAnnot'].T
# print('mapAnnot\n',mapAnnot.shape)

iPoint = np.ones(shape=mapAnnot.shape, dtype=np.bool8)
# print(mapAnnot.size)

# Remove the data points that are outside the region of interest
for indx in range(mapAnnot.size):
    if (mapAnnot[indx]<woi[indx,0]) or (mapAnnot[indx]>woi[indx,1]):
        iPoint[indx] = False
    else:
        iPoint[indx] = True

# print(np.all(iPoint))

# For each mapping point, n, find the voltage amplitude
#  - Apply the user defined function (for example max(egm) – min(egm))

print(iPoint)

iPoint = np.where(iPoint==False,'nan',iPoint)
data[iPoint] = np.where
for item in iPoint:
    print(item)


# print('iPoint-shape',np.where(iPoint=='False','nan',iPoint))






# data[iPoint:,] = np.where(iPoint=='False','nan',data[iPoint:,])

# print(data[np.invert(iPoint)].shape)
# print(data[np.where(~iPoint):])

# for item in c:
#     print(item)

# check for false iPoint and replace the data corresponding to the index of iPoint to False
