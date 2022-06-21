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

"""
Interpolators - :mod:`openep.case.interpolators`
================================================

This module provides classes for performing interpolation.
"""

from attr import attrs

import numpy as np
from .case_routines import calculate_distance

__all__ = [
    'LocalSmoothingInterpolator',
]


@attrs(auto_attribs=True, auto_detect=True)
class LocalSMoothingInterpolator:
    """Interpolator for performing local  smoothing.

    Args:
        points (np.ndarray): Data point coordinates.
        field (np.ndarray): Scalar values for each point.
        smoothing_length (int): Cutoff distance between `points` and coordinates
            at which the interpolant is evaluated. Typically, values in the range
            5-10 are reasonable. Defaults to 5.
        fill_value (float): Value used to assign to the field at coordinates that
            fall outside the query radius.
    """

    points: np.ndarray
    field: np.ndarray
    smoothing_length: int = 5
    fill_value: float = np.NaN

    def __call__(self, new_points):
        """Evaluate the interpolant.

        Args:
            new_points (np.ndarray): Coordinates at which to evaluate the interpolant.

        Returns:
            y (np.ndarray): Interpolated field at `new_points`.
        """

        n_points = len(new_points)
        new_field = np.full(n_points, fill_value=self.fill_value, dtype=float)
        new_field_gradient = np.full((n_points, 3), fill_value=self.fill_value, dtype=float)

        distances = calculate_distance(
            origin=new_points,
            destination=self.points,
        )

        within_cutoff = distances < self.smoothing_length
        include = np.any(within_cutoff, axis=1)

        for index, (point, distance) in enumerate(zip(new_points, distances)):
              
            if not include[index]:
                continue

            # Calculate field at new points
            dx = self.points[distance < self.smoothing_length] - point
            exponent = distance[distance < self.smoothing_length] / self.smoothing_length
            weights = np.exp(-exponent**2)
            normalised_weights = weights / weights.sum()

            field_value = sum(self.field[distance < self.smoothing_length] * normalised_weights)
            new_field[index] = field_value

            # Calculate gradient of field at new points
            mat = np.zeros((3, 3), dtype=float)
            for weight, d in zip(normalised_weights, dx):
                mat += np.dot(
                    weight * d[np.newaxis, :].T,
                    d[np.newaxis, :],
                )
            mat_inverse = np.linalg.pinv(mat)

            diff = (self.field[distance < self.smoothing_length] - field_value)[np.newaxis, :].T
            field_gradient =  diff * np.dot(dx, mat_inverse) * normalised_weights[np.newaxis, :].T
            new_field_gradient[index] = np.sum(field_gradient, axis=0)

        return new_field
