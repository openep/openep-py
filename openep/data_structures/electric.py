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

"""Module containing classes for storing electrogram data."""

from attr import attrs
import numpy as np

__all__ = []


@attrs(auto_attribs=True, auto_detect=True)
class Electrogram:
    """
    Class for storing information about electrograms

    Args:
        egm (np.ndarray): Electrograms
        points (np.ndarray): 3D coordinates of the associated mapping points.
        voltage (np.ndarray): The voltage at each mapping point.
        gain (np.ndarray): gain to apply to each signal
        names (np.ndarray): Names of the associated electrodes.
    """

    egm: np.ndarray
    points: np.ndarray = None
    voltage: np.ndarray = None
    gain: np.ndarray = None
    names: np.ndarray = None

    def __attrs_post_init__(self):
        if self.egm is not None and self.gain is None:
            self.gain = np.zeros(self.egm.shape[0])

    def __repr__(self):
        n_points = len(self.egm) if self.egm is not None else 0
        return f"Electrograms with {n_points} mapping points."


@attrs(auto_attribs=True, auto_detect=True)
class ECG:
    """
    Class for storing information about ECGs

    Args:
        egm (np.ndarray): ECGs
        gain (np.ndarray): gain to apply to each signal
    """

    ecg: np.ndarray
    gain: np.ndarray = None

    def __attrs_post_init__(self):
        if self.ecg is not None and self.gain is None:
            self.gain = np.zeros(self.ecg.shape[0])

    def __repr__(self):
        n_points = len(self.ecg) if self.ecg is not None else 0
        return f"ECGs with {n_points} signals."


@attrs(auto_attribs=True, auto_detect=True)
class Impedance:
    """
    Class for storing information about impedance traces.

    Args:
        times (np.ndarray): Times at which the traces were taken.
        values (np.ndarray): Impedance measurements.
    """

    times: np.ndarray
    values: np.ndarray = None

    def __repr__(self):
        n_traces = len(self.values) if self.values is not None else 0
        return f"Impedance measurements with {n_traces} traces."


@attrs(auto_attribs=True, auto_detect=True)
class ElectricSurface:
    """
    Class for storing information about electroanatomic data projected onto the surface mesh.

    Args:
        nearest_point (np.ndarray): The position of each mapping point projected to the nearest
            point on the 3D surface constructed by the clinical mapping system.
        normals (np.ndarray): Normal to the surface at each of the nearest points.
    """

    nearest_point: np.ndarray
    normals: np.ndarray = None

    def __repr__(self):
        n_points = len(self.nearest_point) if self.nearest_point is not None else 0
        return f"ElectricSurface with {n_points} mapping points."


@attrs(auto_attribs=True, auto_detect=True)
class Annotations:
    """
    Class for storing information about activation times for electrograms.

    Args:
        window_of_interest (np.ndarray): The window of interest for each mapping point
        local_activation_time (np.ndarray): The local activation time for each mapping point
        reference_activation_time (np.ndarray): The reference activation time for each mapping point
    """

    window_of_interest: np.ndarray
    local_activation_time: np.ndarray
    reference_activation_time: np.ndarray

    def __repr__(self):
        n_points = len(self.window_of_interest) if self.window_of_interest is not None else 0
        return f"Annotations with {n_points} mapping points."


@attrs(auto_attribs=True, auto_detect=True)
class Electric:
    """
    Class for storing electrical data obtained during a clinical mapping procedure.

    Args:
        names (np.ndarray): The physician-visible names of the points applied during the
            clinical case.
        internal_names (np.ndarray): The internal names of the points used by the clinical
            electroanatomic mapping system.
        bipolar_egm (Electrogram): The bipolar electrograms, coordinates of mapping points,
            names of the electrodes, and voltages.
        unipolar_egm (Electrogram): The unipolar electrograms, coordinates of the two unipole
            electrodes, names of the two unipole electrodes, and voltages.
        reference_egm (Electrogram): The reference electrograms for each mapping point.
        ecg (np.ndarray): The surface electrocardiogram for each mapping point.
        impedance (Impedance): The time and values of impedance traces for each mapping point
        surface (ElectricSurface): The position of each mapping point projected to the nearest
            point on the 3D surface, as well as the normal to the surface at that point.
        annotations (Annotations): The window of interest, local activation time, and reference
            activation time for each mapping point.
    """

    names: np.ndarray
    internal_names: np.ndarray
    bipolar_egm: Electrogram
    unipolar_egm: Electrogram
    reference_egm: Electrogram
    ecg: np.ndarray
    impedance: Impedance
    surface: ElectricSurface
    annotations: Annotations

    def __repr__(self):
        n_points = len(self.unipolar_egm) if self.unipolar_egm is not None else 0
        return f"Electric data for {n_points} mapping points."


def extract_electric_data(electric_data):
    """Extract electric data from a dictionary.

    Args:
        electric_data (dict): Dictionary containing numpy arrays that describe the
            electric data associated with electrograms taken at various mapping points.

    Returns:
        electric (Electric): object containing electric data associated with electrograms
            taken at various mapping points.
    """

    names = electric_data['tags'].astype(str)
    internal_names = electric_data['names'].astype(str)

    # TODO: check if gain values are stored in electric_data before creating the default values
    bipolar_egm = Electrogram(
        egm=electric_data['egm'],
        points=electric_data['egmX'],
        voltage=electric_data['voltages']['bipolar'],
        gain=np.full_like(electric_data['voltages']['bipolar'], fill_value=1.0),
        names=electric_data['electrodeNames_bip'],
    )
    unipolar_egm = Electrogram(
        egm=electric_data['egmUni'],
        points=electric_data['egmUniX'],
        voltage=electric_data['voltages']['unipolar'],
        gain=None,
        names=electric_data['electrodeNames_uni'],
    )
    reference_egm = Electrogram(
        egm=electric_data['egmRef'],
        gain=np.full(electric_data['egmRef'].shape[0], fill_value=-4.0),
    )

    ecg = ECG(
        ecg=electric_data['ecg'],
        gain=np.full(electric_data['ecg'].shape[0], fill_value=1.0),
    )

    impedance = Impedance(
        times=electric_data['impedances']['time'],
        values=electric_data['impedances']['value'],
    )

    surface = ElectricSurface(
        nearest_point=electric_data['egmSurfX'],
        normals=electric_data['barDirection'],
    )

    annotations = Annotations(
        window_of_interest=electric_data['annotations']['woi'],
        local_activation_time=electric_data['annotations']['mapAnnot'],
        reference_activation_time=electric_data['annotations']['referenceAnnot'],
    )

    electric = Electric(
        names=names,
        internal_names=internal_names,
        bipolar_egm=bipolar_egm,
        unipolar_egm=unipolar_egm,
        reference_egm=reference_egm,
        ecg=ecg,
        impedance=impedance,
        surface=surface,
        annotations=annotations
    )

    return electric


def empty_electric():
    """Create an empty Ablation object with empty numpy arrays.

    Returns:
        electric (Electric): object containing electric data associated with electrograms
            taken at each mapping point.
    """

    names = None
    internal_names = None

    bipolar_egm = Electrogram(
        egm=None,
        points=None,
        voltage=None,
        gain=None,
        names=None,
    )
    unipolar_egm = Electrogram(
        egm=None,
        points=None,
        voltage=None,
        gain=None,
        names=None,
    )
    reference_egm = Electrogram(
        egm=None,
        gain=None,
    )

    ecg = None

    impedance = Impedance(
        times=None,
        values=None,
    )

    surface = ElectricSurface(
        nearest_point=None,
        normals=None,
    )

    annotations = Annotations(
        window_of_interest=None,
        local_activation_time=None,
        reference_activation_time=None,
    )

    electric = Electric(
        names=names,
        internal_names=internal_names,
        bipolar_egm=bipolar_egm,
        unipolar_egm=unipolar_egm,
        reference_egm=reference_egm,
        ecg=ecg,
        impedance=impedance,
        surface=surface,
        annotations=annotations
    )

    return electric
