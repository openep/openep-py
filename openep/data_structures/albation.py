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

"""Module containing classes for storing ablation data."""

from attr import attrs
import numpy as np

__all__ = []


@attrs(auto_attribs=True, auto_detect=True)
class AblationForce:
    """Class for storing data on ablation force.

    Args:
        times (np.ndarray): array of shape N
        force (np.ndarray): array of shape N
        axial_angle (np.ndarray): array of shape N
        lateral_angle (np.ndarray): array of shape N
        points (np.ndarray): array of shape Nx3

    """
    times: np.ndarray
    force: np.ndarray
    axial_angle: np.ndarray
    lateral_angle: np.ndarray
    points: np.ndarray

    def __repr__(self):
        return f"Ablation forces with {len(self.times)} sites."


@attrs(auto_attribs=True, auto_detect=True)
class Ablation:
    """
    Class for storing ablation data.

    Args:
        times (np.ndarray): array of shape N
        power (np.ndarray): array of shape N
        impedance (np.ndarray): array of shape N
        temperature (np.ndarray): array of shape N
        force (AblationForce): data on the force at each ablation site. Specifically, the
            force applied, the time at which it was applied, and axial and lateral angles,
            and the 3D coordinates of the ablation site.
    """

    times: np.ndarray
    power: np.ndarray
    impedance: np.ndarray
    temperature: np.ndarray
    force: AblationForce

    def __repr__(self):
        return f"Ablations with {len(self.times)} ablation sites."


def extract_ablation_data(ablation_data):
    """Extract surface data from a dictionary.

    Args:
        ablation_data (dict): Dictionary containing numpy arrays that describe the
            ablation sites.

    Returns:
        ablation (Ablation): times, power, impedance and temperature for each ablation site,
            as well as the force applied.
    """

    times = ablation_data['originaldata']['ablparams']['time']
    power = ablation_data['originaldata']['ablparams']['power']
    impedance = ablation_data['originaldata']['ablparams']['impedance']
    temperature = ablation_data['originaldata']['ablparams']['distaltemp']

    force = AblationForce(
        times=ablation_data['originaldata']['force']['time'],
        force=ablation_data['originaldata']['force']['force'],
        axial_angle=ablation_data['originaldata']['force']['axialangle'],
        lateral_angle=ablation_data['originaldata']['force']['lateralangle'],
        points=ablation_data['originaldata']['force']['position'],
    )
    
    ablation = Ablation(
        times=times,
        power=power,
        impedance=impedance,
        temperature=temperature,
        force=force,
    )

    return ablation
