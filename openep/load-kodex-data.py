import numpy as np
import openep
from openep._datasets.openep_datasets import DATASET_2

case = openep.load_opencarp(
    '../openep-misc/examples/data/pig21_endo_coarse.pts',
    '../openep-misc/examples/data/pig21_endo_coarse.elem'
)
egms = np.loadtxt('../openep-misc/examples/data/pig21_endo_coarse_phie.dat')
case.add_unipolar_electrograms(
    egms,
)
openep.export_openep_mat(case, 'test-export.mat')
case = openep.load_openep_mat('test-export.mat')


#case = openep.load_openep_mat('../openep-misc/examples/data/pig21_endo_coarse.mat')
#case = openep.load_openep_mat('../testdata.mat')

#openep.export_openep_mat(case, "../testdata-exported.mat")
#new_case = openep.load_openep_mat('../testdata-exported.mat')
#print(case)
