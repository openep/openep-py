# Code written by Jerin Rajan on 18th Mar 2021
# read mat file

import os
from pathlib import Path
import scipy.io as sio
import numpy as np
import pandas as pd
import trimesh

# relative path to the sample dataset
# Main MAT File - userdata


# print('cwd:',os.getcwd())
# os.chdir('../repos')

file_path = '../openep/openep-examples-main/openep_dataset_2.mat'
data_tri_X_lst = list()
data_tri_T_lst = list()

# Load the file
main_file = sio.loadmat(file_path)

# Parsing TriRep Data
# Anatomy
# .X & .Triangulation

# Vertices/Points
# data_tri_X = pd.DataFrame(main_file['userdata']['surface_t[riRep_X'])
data_tri_X = main_file['userdata']['surface_triRep_X'][0][0]
for item in data_tri_X:
    data_tri_X_lst.append(item.tolist())

print('size_of_X:',len(data_tri_X_lst))
# Triangulation/faces
# data_tri_T = pd.DataFrame(main_file['userdata']['surface_triRep_Triangulation'][0][0])
data_tri_T = main_file['userdata']['surface_triRep_Triangulation'][0][0]

for item in data_tri_T:
    data_tri_T_lst.append(item.tolist())

print('size_of_X:',len(data_tri_T_lst))

# Convert to a Trimesh object
mesh = trimesh.Trimesh(vertices=data_tri_X_lst,faces=data_tri_T_lst)






# Voltage plots


# Colour Pallette

# Fill Threshold - 0.5mV


# Voltage threshold









# DRAWMAP plots an OpenEP map
# userdata is a Carto data structure
def draw_map(
        userdata, *arg):
    pass
