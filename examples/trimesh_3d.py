# Code written by Jerin Rajan on 06th May 2021
'''
Class to create a Triangular Mesh object

    Parameters:
    -----------
    x -

'''

import trimesh as tm
import numpy as np

'''
Parameters
==========
 x (vertices) : (n,3) float
    Array of vertex locations
 t (faces) : (m, 3) int
    Arry of triangular faces (triangulated)
'''

class tri_mesh:
    def __init__(self,x,t):
        self.x = x
        self.t = t

    def points(self):
        return self.x

    def faces(self):
        return self.t
