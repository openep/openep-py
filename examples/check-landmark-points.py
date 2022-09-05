import numpy as np

import openep
from openep._datasets.openep_datasets import DATASET_2


# Add landmark to carto dataset
case = openep.load_openep_mat(DATASET_2)
case.add_landmark('a', 'b', point=np.array([0, 0, 0]))
openep.export_openep_mat(
    case,
    'test-add-landmark.mat',
)
case_with_landmark = openep.load_openep_mat('test-add-landmark.mat')

# Add a landmark to kodex dataset
kodex_case = openep.load_openep_mat(
    '/Users/paul/Downloads/2021-08-05_13-21-51_PL_L465_15_RA_LAT_WT_Abl - RA_Map_Remap1_1_1_1_1_1_1_1.mat'
)
kodex_case.add_landmark('a', 'b', point=np.array([0, 0, 0]))
openep.export_openep_mat(
    kodex_case,
    '/Users/paul/Downloads/2021-08-05_13-21-51_PL_L465_15_RA_LAT_WT_Abl - RA_Map_Remap1_1_1_1_1_1_1_1-landmarks.mat'
)
kodex_case_with_landmark = openep.load_openep_mat(
    '/Users/paul/Downloads/2021-08-05_13-21-51_PL_L465_15_RA_LAT_WT_Abl - RA_Map_Remap1_1_1_1_1_1_1_1-landmarks.mat'
)

# Add landmark to openCARP dataset
carp = openep.load_opencarp(
    points="/Users/paul/github/openep-misc/examples/data/pig21_endo_coarse.pts",
    indices="/Users/paul/github/openep-misc/examples/data/pig21_endo_coarse.elem",
    scale_points=1e-3,
)
unipolar = np.loadtxt("/Users/paul/github/openep-misc/examples/data/pig21_endo_coarse_phie.dat")
carp.add_unipolar_electrograms(unipolar=unipolar)
carp.add_landmark('a', 'b', point=np.array([0, 0, 0]))

openep.export_openep_mat(carp, 'test-add-landmark-carp.mat')
carp_with_landmark = openep.load_openep_mat('test-add-landmark-carp.mat')

print(case, case_with_landmark)
print(kodex_case, kodex_case_with_landmark)
print(carp, carp_with_landmark)
