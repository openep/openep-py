
import numpy as np
from sklearn.neighbors import NearestNeighbors
from scipy.interpolate import LinearNDInterpolator as linterp
from scipy.interpolate import NearestNDInterpolator as nearest

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


def distBetweenPoints(A, B):
    '''
    DISTBETWEENPOINTS returns the distance from A to B. A and B are specified
    as row vectors [x, y, z] or matrices, with rows representing different
    points. If npoints in A and B are different A must specify one and only 
    one point.

    Args:
        A: row vectors of size mx3 
        B: row vectors of size mx3

    Returns:
        returns an array of distance between the points in the row vectors 

    '''

    diffsq = np.square(np.subtract(A,B))
    d = np.sqrt(np.sum(diffsq,axis=1)).reshape(B.shape[0],1)

    return d


class LinearNDInterpolatorExt(object):
    ''' 
    Interpolation Method - Python implementation of ScatteredInterpolation
    '''

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


class openEpDataInterpolator():

    '''
    OPENEPDATAINTERPOLATOR Creates objects for performing spatial
    interpolation for OpenEP data
    '''

    method = 'scatteredinterpolant'
    def __init__(self,method,distanceThreshold):
        self.method = method
        self.distanceThreshold = distanceThreshold


    def interpolate(self,x0,d0,x1,*args):
        self.x0 = x0
        self.d0 = d0
        self.x1 = x1

        F = LinearNDInterpolatorExt(points=(self.x0[:,0],self.x0[:,1],self.x0[:,2]),values=self.d0)
        d1 = np.array(F(self.x1[:,0],self.x1[:,1],self.x1[:,2])).reshape(self.x1.shape[0],1)

        # Workout if there are any points on the surface that are < 0
        d1[d1<0] = 0

        #  work out which points on the surface are too far away from real data
        # Remove any interpolated values which are outwith the fill threshold 

        neigh = NearestNeighbors(n_neighbors=1, 
                                radius=1, 
                                algorithm='auto', 
                                leaf_size=30,
                                metric='minkowski',
                                p=2)

        neigh.fit(self.x0)

        id = neigh.kneighbors(self.x1,return_distance=False)
        cPts = self.x0[id]
        cPts = np.array(cPts[:,0])
        
        d = distBetweenPoints(cPts,self.x1)

        thresholdDistance = np.zeros(shape=d.shape,dtype=np.bool)
        thresholdDistance[d>self.distanceThreshold] = 1
        d1[thresholdDistance] = 'nan'
        d1 = d1.flatten()

        return d1
        