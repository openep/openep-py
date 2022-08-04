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
    times: np.ndarray = None
    force: np.ndarray = None
    axial_angle: np.ndarray = None
    lateral_angle: np.ndarray = None
    points: np.ndarray = None

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

    times: np.ndarray = None
    power: np.ndarray = None
    impedance: np.ndarray = None
    temperature: np.ndarray = None
    force: AblationForce = None

    def __attrs_post_init__(self):

        if self.force is None:
            self.force = AblationForce()

    def __repr__(self):
        n_sites = {len(self.times)} if self.times is not None else 0
        return f"Ablations with {n_sites} ablation sites."


def extract_ablation_data(ablation_data):
    """Extract surface data from a dictionary.

    Args:
        ablation_data (dict): Dictionary containing numpy arrays that describe the
            ablation sites.

    Returns:
        ablation (Ablation): times, power, impedance and temperature for each ablation site,
            as well as the force applied.
    """

    if isinstance(ablation_data, np.ndarray) or ablation_data['originaldata']['ablparams']['time'].size == 0:
        return Ablation()

    times = ablation_data['originaldata']['ablparams']['time'].astype(float)
    power = ablation_data['originaldata']['ablparams']['power'].astype(float)
    impedance = ablation_data['originaldata']['ablparams']['impedance'].astype(float)
    temperature = ablation_data['originaldata']['ablparams']['distaltemp'].astype(float)

    force = AblationForce(
        times=ablation_data['originaldata']['force']['time'].astype(float),
        force=ablation_data['originaldata']['force']['force'].astype(float),
        axial_angle=ablation_data['originaldata']['force']['axialangle'].astype(float),
        lateral_angle=ablation_data['originaldata']['force']['lateralangle'].astype(float),
        points=ablation_data['originaldata']['force']['position'].astype(float),
    )

    ablation = Ablation(
        times=times,
        power=power,
        impedance=impedance,
        temperature=temperature,
        force=force,
    )

    return ablation


def empty_ablation():
    """Create an empty Ablation object with empty numpy arrays.

    Returns:
        ablation (Ablation): times, power, impedance and temperature for each ablation site,
            as well as the force applied.
    """

    force = AblationForce(
        times=np.array([], dtype=float),
        force=np.array([], dtype=float),
        axial_angle=np.array([], dtype=float),
        lateral_angle=np.array([], dtype=float),
        points=np.array([], dtype=float),
    )

    ablation = Ablation(
        times=np.array([], dtype=float),
        power=np.array([], dtype=float),
        impedance=np.array([], dtype=float),
        temperature=np.array([], dtype=float),
        force=force,
    )

    return ablation
