import openep
import sys
import os
import math
import pyvista as pv
import numpy as np
import scipy.io
import scipy
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from sklearn.neighbors import KDTree
import pdb
from scipy.interpolate import Rbf
from vedo import *
import pyvistaqt
from sklearn import metrics

from ..case.case_routines import (
    interpolate_general_cloud_points_onto_surface,
    interpolate_activation_time_onto_surface
)

class Divergence():


    def __init__(self, case):
        self.case = case #.copy() change to copy later
        self.include = None,
        self.direction = None,
        self.mesh_original = None,
        self.values = None
        self.divergence_per_point=None

    def calculate(
            self,
            output_binary_field=False,
            method_kws=dict()
    ):
        """Calculate conduction velocity and interpolates each point on a mesh,
            stores in conduction_velocity fields.

        Args:
            method (str): method to use for interpolation. This should be one of:
                - 'triangualtion'
                - 'plane_fitting'
                - 'rbf'
        """

        out = self.calculate_divergence(method_kws)


    def calculate_divergence(self, method_kws, output_binary_field=False):

        if output_binary_field:
            collision_threshold = method_kws.get("collision_threshold",-1)
            focal_threshold = method_kws.get("collision_threshold", 1)

        surface = self.case.create_mesh()
        egmX=self.case.electric.bipolar_egm.points
        LAT=self.case.electric.annotations.local_activation_time - self.case.electric.annotations.reference_activation_time
        points=surface.points
        face=surface.faces

        egmX_point_cloud = pv.PolyData(egmX)
        egmX_point_cloud['LAT[ms]'] = LAT
        mesh_original = pv.PolyData(points, face)

        # interpolated_old = mesh_original.interpolate(egmX_point_cloud, radius=5, n_points=None)
        interpolated_scalar=interpolate_general_cloud_points_onto_surface(
            case=self.case,
            cloud_values=LAT,
            cloud_points=egmX,
        )

        mesh_original['LAT_scalar'] = interpolated_scalar
        deriv = mesh_original.compute_derivative(scalars='LAT_scalar')

        cv_x = deriv['gradient'][:, 0] / np.sum(np.power(deriv['gradient'], 2), axis=1)
        cv_y = deriv['gradient'][:, 1] / np.sum(np.power(deriv['gradient'], 2), axis=1)
        cv_z = deriv['gradient'][:, 2] / np.sum(np.power(deriv['gradient'], 2), axis=1)

        direction = np.ndarray(shape=(len(cv_x), 3))
        cv_direction = np.ndarray(shape=(len(cv_x), 3))
        for i in range(len(cv_x)):
            direction[i][0] = cv_x[i] / np.sqrt(
                np.sum(np.power(cv_x[i], 2) + np.power(cv_y[i], 2) + np.power(cv_z[i], 2)))
            direction[i][1] = cv_y[i] / np.sqrt(
                np.sum(np.power(cv_x[i], 2) + np.power(cv_y[i], 2) + np.power(cv_z[i], 2)))
            direction[i][2] = cv_z[i] / np.sqrt(
                np.sum(np.power(cv_x[i], 2) + np.power(cv_y[i], 2) + np.power(cv_z[i], 2)))
            cv_direction[i][0] = cv_x[i]
            cv_direction[i][1] = cv_y[i]
            cv_direction[i][2] = cv_z[i]

        mesh_original['activation_direction'] = cv_direction
        div = mesh_original.compute_derivative(scalars='activation_direction', divergence=True)
        sargs = dict(
            title_font_size=32,
            label_font_size=32,
            shadow=True,
            n_labels=11,
            italic=True,
            fmt="%.1f",
            font_family="arial",
            color='black'
        )

        if output_binary_field:
            div_value = []
            for i in range(len(div['divergence'])):
                if div['divergence'][i] < collision_threshold or div['divergence'][i] > focal_threshold:
                    div_value.append(1)
                else:
                    div_value.append(0)
        else:
            div_value = div['divergence']

        mesh_divergence = pv.PolyData(points, face)
        mesh_divergence['div'] = div_value
        divergence_per_point = mesh_divergence.cell_data_to_point_data()

        self.case.fields.divergence = div_value
        self.divergence_per_point=divergence_per_point

        self.direction, self.mesh_original, self.values = direction, mesh_original, div_value
