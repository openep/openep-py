# OpenEP
# Copyright (c) 2021 OpenEP Collaborators
#
# This file is part of OpenEP.
#
# OpenEP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OpenEP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program (LICENSE.txt).  If not, see <http://www.gnu.org/licenses/>

"""Module containing analysis classes"""
from ._conduction_velocity import *
from ..case.case_routines import interpolate_general_cloud_points_onto_surface


class Analyse:
    """
    Analyse class for managing analysis of case data.

    Args:
        case (Case): The case data to be analysed.

    Attributes:
        conduction_velocity (ConductionVelocity): An instance of ConductionVelocity for conducting
                                                  velocity analysis on the case.

        divergence (Divergence): An instance of Divergence for conducting divergence analysis on the case.

    """
    def __init__(self, case):
        self.conduction_velocity = ConductionVelocity(case)
        self.divergence = Divergence(case)


class ConductionVelocity:
    """
    Class for calculating and analyzing conduction velocity.

    This class provides methods to calculate conduction velocity from case data and optionally
    interpolate these values onto a mesh. It handles different methods of calculation and manages
    the resulting data.

    Attributes:
        _case (Case): Internal reference to the case data.
        values (np.ndarray): The calculated conduction velocity values.
        centers (np.ndarray): The centers corresponding to the conduction velocity values.

    Parameters:
        case (Case): The case data from which conduction velocity is to be calculated.
    """
    def __init__(self, case):
        self._case = case
        self.values = None
        self.centers = None

    def calculate_cv(
            self,
            method='triangulation',
            include=None,
            apply_scalar_field=True,
            interpolation_kws=None,
            **method_kwargs
    ):
        """
        Calculate the conduction velocity and interpolate points onto a mesh.

        This method calculates conduction velocities using one of the specified methods and
        optionally interpolates the calculated values onto a mesh, storing the results in the
        'conduction_velocity' field of the case object.

        Args:
            method (str): The method to use for calculating conduction velocity. Options are:
                          - 'triangulation'
                          - 'plane_fitting'
                          - 'rbf'
                          Defaults to 'triangulation'.

            include (np.ndarray, optional): A boolean array specifying which data points to include
                                        in the preprocessing for bipolar EGM data and LAT.
                                        Defaults to "includes" specified in the imported data.

            apply_scalar_field (bool, optional): If True, applies the calculated conduction velocity
                                                 values as a scalar field on the case object's mesh.
                                                 Defaults to True.

            interpolation_kws (dict, optional): Keyword arguments for interpolation function.
                                    > interpolate_general_cloud_points_onto_surface(**interpolation_kws)
                                    Defaults to None.

            **method_kwargs: Additional keyword arguments specific to the chosen calculation method.

        Returns:
            tuple: A tuple containing two elements:
                   - values (np.ndarray): An array of calculated conduction velocity values.
                   - centers (np.ndarray): An array of centers corresponding to the calculated values.

        Raises:
            ValueError: If the specified method is not among the available options.

        Example:
            >> cv = ConductionVelocity(case)
            >> values, centers = cv.calculate(method='plane_fitting', apply_scalar_field=True)
        """

        supported_cv_methods = {
            "triangulation": triangulation,
            "plane_fitting": plane_fitting,
            "rbf": radial_basis_function
        }

        if method.lower() not in supported_cv_methods:
            raise ValueError(f"`method` must be one of {supported_cv_methods.keys()}.")

        if include is None and self._case.electric.include is None:
            raise TypeError(f"include object is of 'NoneType'")
        else:
            include = self._case.electric.include.astype(bool) if include is None else include

        interpolation_kws = dict() if interpolation_kws is None else interpolation_kws

        bipolar_egm_pts = self._case.electric.bipolar_egm.points[include]
        lat = (self._case.electric.annotations.local_activation_time[include]
               - self._case.electric.annotations.reference_activation_time[include])

        cv_method = supported_cv_methods[method]
        self.values, self.centers = cv_method(bipolar_egm_pts, lat, **method_kwargs)

        if apply_scalar_field:
            self._case.fields.conduction_velocity = interpolate_general_cloud_points_onto_surface(
                case=self._case,
                cloud_values=self.values,
                cloud_points=self.centers,
                **interpolation_kws
            )

        return self.values, self.centers


class Divergence:
    """
    Class for calculating and analyzing divergence of conduction velocity.

    This class provides methods to calculate the divergence of conduction velocity from case data.
    It can output binary fields indicating regions of divergence and manage the resulting data.

    Attributes:
        _case (Case): Internal reference to the case data.
        direction (np.ndarray): The calculated divergence direction vectors.
        values (np.ndarray): The calculated divergence values.

    Args:
        case (Case): The case data from which divergence is to be calculated.
    """
    def __init__(self, case):
        self._case = case
        self.direction = None
        self.values = None

    def calculate_divergence(
        self,
        include=None,
        output_binary_field=False,
        apply_scalar_field=True,
        interpolation_kws=None,
    ):
        """
        Calculate the divergence of conduction velocity and optionally apply the result as a scalar field.

        Parameters:
            include (np.ndarray, optional): A boolean array specifying which data points to include in the
                                            calculation. If None, the 'include' array from the case's electric
                                            data is used. Defaults to None.

            output_binary_field (bool, optional): If True, the output will be a binary field where 1 indicates
                                                  regions of divergence and 0 otherwise.
                                                  Defaults to False.

            apply_scalar_field (bool, optional): If True, the calculated divergence values (or binary field, if
                                                 'output_binary_field' is True) are applied as a scalar field to
                                                 the case object.
                                                 Defaults to True.

            interpolation_kws (dict, optional): Keyword arguments for interpolation function.
                                    > interpolate_general_cloud_points_onto_surface(**interpolation_kws)
                                    Defaults to None.

        Returns:
            tuple of (np.ndarray, np.ndarray): A tuple containing the divergence direction and divergence values.
                                              The direction is a 2D array of shape (N, 3), where N is the number
                                              of points, and each row represents the (x, y, z) components of the
                                              divergence direction. The values are a 1D array of divergence values.

        Example:
            case = load_case_data('path/to/case_data')
            cv = ConductionVelocity(case)
            direction, values = cv.calculate_divergence(output_binary_field=True, apply_scalar_field=True)
        """

        if include is None and self._case.electric.include is None:
            raise TypeError(f"include object is of 'NoneType'")
        else:
            include = self._case.electric.include.astype(bool) if include is None else include

        bipolar_egm_pts = self._case.electric.bipolar_egm.points[include]
        lat = (self._case.electric.annotations.local_activation_time[include]
               - self._case.electric.annotations.reference_activation_time[include])

        self.direction, self.values = divergence(
            case=self._case,
            bipolar_egm_pts=bipolar_egm_pts,
            local_activation_time=lat,
            output_binary_field=output_binary_field,
            interpolation_kws=interpolation_kws,
        )

        if apply_scalar_field:
            self._case.fields.cv_divergence = self.values

        return self.direction, self.values
