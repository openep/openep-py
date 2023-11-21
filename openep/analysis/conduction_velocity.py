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
    interpolate_general_cloud_points_onto_surface
)

class ConductionVelocity():


    def __init__(self, case):
        self.case = case
        self.prep_egmX, self.prep_LAT = self._preprocess_data()
        self.include = None
        self.values = None
        self.points = None

    def calculate(
            self,
            method='triangulation',
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
        NAME_TO_FUNC = {
            "triangulation": self.calculate_with_triangualtion,
            "plane_fitting": self.calculate_with_plane_fitting,
            "rbf": self.calculate_with_rbf,
        }

        method = method.lower()
        if method not in NAME_TO_FUNC:
            raise ValueError(f"`method` must be one of {NAME_TO_FUNC.keys()}.")

        cv_calc_method = NAME_TO_FUNC[method]
        out = cv_calc_method(method_kws)

        self.case.fields.conduction_velocity = interpolate_general_cloud_points_onto_surface(
            case=self.case,
            cloud_values=self.case.analyse.conduction_velocity.values,
            cloud_points=self.case.analyse.conduction_velocity.points
        )
        return out

    def _preprocess_data(self):
        egmX = self.case.electric.bipolar_egm.points
        voltage = self.case.electric.bipolar_egm.voltage
        LAT = self.case.electric.annotations.local_activation_time - self.case.electric.annotations.reference_activation_time
        surface = self.case.create_mesh()

        valid_id, invalid_id = np.where(LAT != -10000), np.where(LAT == -10000)
        egmX_clean = np.delete(egmX, np.where(LAT == -10000), 0)
        LAT_clean = np.delete(LAT, np.where(LAT == -10000))

        return egmX_clean, LAT_clean

    def calculate_with_triangualtion(self, method_kws=dict()):

        print('---> Starting triangulation')
        egmX, LAT = self.prep_egmX, self.prep_LAT
        cv = []
        min_theta = method_kws.get('min_theta', 30)
        max_elec_distance = method_kws.get('max_elec_dist', 10)
        min_elec_distance = method_kws.get('min_elec_dist', 1.5)
        min_lat_difference = method_kws.get('min_lat_diff', 2)
        cv_centers_x = []
        cv_centers_y = []
        cv_centers_z = []
        cv_centers_id = []
        alpha = method_kws.get('alpha', 5)

        egmX_point_cloud = pv.PolyData(egmX)
        print('---> Triangulation process')
        surf = egmX_point_cloud.delaunay_3d(alpha=alpha)
        print('---> Triangulation process')
        print('---> CV calculation (Triangulation method started ...)')

        for num_point in range(len(surf.cells_dict[5])):

            vtx_id = []
            lat = []

            for i in range(3):
                vtx_id.append(surf.cells_dict[5][num_point][i])
                lat.append(LAT[surf.cells_dict[5][num_point][i]])

            id_lat_sorted = np.argsort(lat)

            O = [egmX[int(vtx_id[id_lat_sorted[0]])], lat[id_lat_sorted[0]]]
            A = [egmX[int(vtx_id[id_lat_sorted[1]])], lat[id_lat_sorted[1]]]
            B = [egmX[int(vtx_id[id_lat_sorted[2]])], lat[id_lat_sorted[2]]]

            OA = np.sqrt(sum(np.power(np.subtract(O[0], A[0]), 2)))
            OB = np.sqrt(sum(np.power(np.subtract(O[0], B[0]), 2)))
            AB = np.sqrt(sum(np.power(np.subtract(A[0], B[0]), 2)))

            tOA = A[1] - O[1]
            tOB = B[1] - O[1]

            theta = np.arccos((np.power(OA, 2) + np.power(OB, 2) - np.power(AB, 2)) / (2 * OA * OB))
            # check if the conditions are meet to accept the triangle set as a viable acceptable one
            if (math.degrees(theta) >= min_theta and OA >= min_elec_distance and OA <= max_elec_distance and
                    OB >= min_elec_distance and OB <= max_elec_distance and tOA >= min_lat_difference and tOB >= min_lat_difference):

                alpha = np.arctan((tOB * OA - tOA * OB * np.cos(theta)) / (tOA * OB * np.sin(theta)))
                cv_temp = (OA / tOA) * np.cos(alpha)
                cv.append(cv_temp)
                cv_centers_x.append((O[0][0] + A[0][0] + B[0][0]) / 3)
                cv_centers_y.append((O[0][1] + A[0][1] + B[0][1]) / 3)
                cv_centers_z.append((O[0][2] + A[0][2] + B[0][2]) / 3)
                cv_centers_id.append(vtx_id[id_lat_sorted[0]])
            else:
                cv.append(np.nan)
                cv_centers_x.append((O[0][0] + A[0][0] + B[0][0]) / 3)
                cv_centers_y.append((O[0][1] + A[0][1] + B[0][1]) / 3)
                cv_centers_z.append((O[0][2] + A[0][2] + B[0][2]) / 3)
                cv_centers_id.append(vtx_id[id_lat_sorted[0]])
            ####### Create a triangulation mesh from egmX pointcloud #####
        cv_centers = [cv_centers_x, cv_centers_y, cv_centers_z]
        print('---> CV calculation ended')
        print(cv)

        egmX_point_cloud = pv.PolyData(np.transpose(cv_centers))
        egmX_point_cloud['CV[m/s]'] = cv
        # interpolated_mesh_cv = self.mesh_original.interpolate(egmX_point_cloud, radius=10, n_points=None, null_value=np.nan)


        egmX_point_cloud = pv.PolyData(np.transpose(cv_centers))
        egmX_point_cloud['CV[m/s]'] = cv


        # interpolated_cv_field = interpolate_conduction_velocity_onto_surface(
        #     self.case,
        #     np.array(cv),
        #     egmX_point_cloud
        # )
        # interpolated_mesh_cv = self.mesh_original.interpolate(egmX_point_cloud, radius=10, n_points=None, null_value=np.nan)

        self.values = np.array(cv)
        self.points = egmX_point_cloud.points

        return np.array(cv), egmX_point_cloud.points

    def calculate_with_plane_fitting(self, method_kws=dict()):
        egmX, LAT = self.prep_egmX, self.prep_LAT

        cv = []
        cv_centers_x = []
        cv_centers_y = []
        cv_centers_z = []
        cv_centers_id = []
        cent = np.random.random((len(egmX), 3))
        direction = np.random.random((len(egmX), 3))

        leaf_size = method_kws.get('leaf_size', 5)
        tree = KDTree(egmX, leaf_size=leaf_size)
        all_nn_indices = tree.query_radius(egmX, r=10)  # 10 closest neighnours
        all_nns = [[egmX[idx] for idx in nn_indices] for nn_indices in all_nn_indices]

        min_n_nearest_neighbours = method_kws.get('min_n_nearest_neighbours', 5)
        print(len(egmX))
        for k in range(len(egmX)):
            print(k, end='\r')
            if len(all_nns[k]) >= min_n_nearest_neighbours:
                cv_data = np.zeros(shape=(len(all_nns[k]) + 1, 4))
                cv_data[0][:3] = egmX[k]
                cv_data[0][3] = LAT[k]
                for i in range(len(all_nns[k])):
                    cv_data[i + 1][:3] = all_nns[k][i]
                    cv_data[i + 1][3] = LAT[all_nn_indices[k][i]]

                x1 = cv_data[:, 0]
                x2 = cv_data[:, 1]
                x3 = cv_data[:, 2]
                m = cv_data[:, 3]

                res = minimize(self.func, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], args=(x1, x2, x3, m))
                a, b, c, d, e, f, g, h, q, l = res.x

                x = egmX[k, 0]
                y = egmX[k, 1]
                z = egmX[k, 2]

                tx = 2 * a * x + d * y + e * z + g
                ty = 2 * b * y + d * x + f * z + h
                tz = 2 * c * z + e * x + f * y + q

                vx = tx / (tx ** 2 + ty ** 2 + tz ** 2)
                vy = ty / (tx ** 2 + ty ** 2 + tz ** 2)
                vz = tz / (tx ** 2 + ty ** 2 + tz ** 2)

                cent[k, 0] = x
                cent[k, 1] = y
                cent[k, 2] = z

                direction[k, 0] = vx
                direction[k, 1] = vy
                direction[k, 2] = vz

                cv_temp = np.sqrt(np.sum(np.power(vx, 2) + np.power(vy, 2) + np.power(vz, 2)))
                direction[k, 0] = vx
                direction[k, 1] = vy
                direction[k, 2] = vz

                if cv_temp >= 0.2 and cv_temp <= 2:
                    cv.append(cv_temp)
                    cv_centers_x.append(x)
                    cv_centers_y.append(y)
                    cv_centers_z.append(z)
                    cv_centers_id.append(k)
                else:

                    cv.append(cv_temp)
                    cv_centers_x.append(x)
                    cv_centers_y.append(y)
                    cv_centers_z.append(z)
            else:
                x = egmX[k, 0]
                y = egmX[k, 1]
                z = egmX[k, 2]
                cv.append(np.nan)
                cv_centers_x.append(x)
                cv_centers_y.append(y)
                cv_centers_z.append(z)
                cv_centers_id.append(k)
            cv_centers = [cv_centers_x, cv_centers_y, cv_centers_z]

        egmX_point_cloud = pv.PolyData(np.transpose(cv_centers))
        egmX_point_cloud['CV[m/s]'] = cv

        # interpolated_cv_field = interpolate_conduction_velocity_onto_surface(
        #     self.case,
        #     np.array(cv),
        #     egmX_point_cloud
        # )
        # interpolated_mesh_cv = self.mesh_original.interpolate(egmX_point_cloud, radius=10, n_points=None, null_value=np.nan)

        self.values = np.array(cv)
        self.points = egmX_point_cloud.points

        return np.array(cv), egmX_point_cloud.points


        # interpolated_mesh_cv = self.mesh_original.interpolate(egmX_point_cloud, radius=10, n_points=None, null_value=np.nan)
        # return interpolated_mesh_cv.point_data.active_scalars

    def func(self, coef, x1, x2, x3, m):
        ''' This the function we use to fit a surface to the activations in plane
        fitting method'''
        a = coef[0]
        b = coef[1]
        c = coef[2]
        d = coef[3]
        e = coef[4]
        f = coef[5]
        g = coef[6]
        h = coef[7]
        q = coef[8]
        l = coef[9]
        # print(np.linalg.norm(m-(a*x1**2 + b*x2**2 + c*x3**2 + d*x1*x2 + e*x1*x3 + f*x2*x3 + g*x1 +h*x2 + q*x3 + l)))
        norm = np.linalg.norm(m - (
                    a * x1 ** 2 + b * x2 ** 2 + c * x3 ** 2 + d * x1 * x2 + e * x1 * x3 + f * x2 * x3 + g * x1 + h * x2 + q * x3 + l))
        return norm

    def calculate_with_rbf(self):
        pass