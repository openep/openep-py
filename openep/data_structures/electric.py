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

from threading import local
from attr import attrs, field
import numpy as np

__all__ = []


class LandmarkPoints:
    """Class for storing information about landmark points.

    Data from **all** mapping points should be passed here. We then filter out
    points that are landmark-only.

    This is to make it easier to make a point a landmark point, or to remove a
    landmark point (i.e. this can be done by changing the flag `LandmarkPoints.is_landmark`)

    Args:
        points (np.ndarray): 3D coordinates of all mapping points.
        names (np.ndarray): User-defined tag for electrode for each mapping point.
        internal_names (np.ndarray): Internal name of each electrode.
    """

    def __init__(
        self,
        points: np.ndarray = None,
        names: np.ndarray = None,
        internal_names: np.ndarray = None,
        is_landmark: np.ndarray = None,
    ):

        self._points = points
        self._names = names
        self._internal_names = internal_names
        self._is_landmark = is_landmark

    @property
    def points(self):
        return self._points[self._is_landmark] if self._points is not None else None

    @property
    def n_points(self):
        return self.points.shape[0] if self.points is not None else 0

    @property
    def names(self):
        return self._names[self._is_landmark] if self._names is not None else None

    @property
    def internal_names(self):
        return self._internal_names[self._is_landmark] if self._internal_names is not None else None

    def __repr__(self):
        return f"Landmarks with {self.n_points} landmark points."


class Electrogram:
    """
    Class for storing information about electrograms

    Args:
        egm (np.ndarray): Electrograms
        points (np.ndarray): 3D coordinates of the associated mapping points.
        voltage (np.ndarray): The voltage at each mapping point.
        gain (np.ndarray): gain to apply to each signal
        names (np.ndarray): Names of the associated electrodes.
        is_electrical (np.ndarray): Flag for whether each point has electrical data
    """

    def __init__(
        self,
        egm: np.ndarray = None,
        points: np.ndarray = None,
        voltage: np.ndarray = None,
        gain: np.ndarray = None,
        names: np.ndarray = None,
        is_electrical: np.ndarray = None,
    ):

        if egm is not None and gain is None:
            gain = np.ones(egm.shape[0])

        if egm is not None and voltage is None:
            voltage = np.full(egm.shape[0], fill_value=np.NaN)

        if egm is not None and names is None:
            names = np.full(egm.shape[0], fill_value='', dtype=str)

        if egm is not None and is_electrical is None:
            is_electrical = np.ones(egm.shape[0], fill_value=True, dtype=bool)

        self._egm = egm
        self._points = points
        self._voltage = voltage
        self._gain = gain
        self._names = names

        # Not all points have electrical signals, we need to ignore these.
        # e.g. landmark points in Kodex have no electrical data.
        self._is_electrical = is_electrical

    @property
    def egm(self):
        return self._egm[self._is_electrical] if self._egm is not None else None

    @property
    def points(self):
        return self._points[self._is_electrical] if self._points is not None else None

    @property
    def n_points(self):
        return self.egm.shape[0] if self.egm is not None else 0

    @property
    def n_samples(self):
        return self.egm.shape[1] if self.egm is not None else 0

    @property
    def voltage(self):
        return self._voltage[self._is_electrical] if self._voltage is not None else None

    @voltage.setter
    def voltage(self, voltage):
        if isinstance(voltage, np.ndarray) and voltage.shape[0] == self.n_points:
            self._voltage[self._is_electrical] = voltage
        else:
            self._voltage = voltage

    @property
    def gain(self):
        return self._gain[self._is_electrical] if self._gain is not None else None

    @gain.setter
    def gain(self, gain):
        if isinstance(gain, np.ndarray) and gain.shape[0] == self.n_points:
            self._gain[self._is_electrical] = gain
        else:
            self._gain = gain

    @property
    def names(self):
        return self._names[self._is_electrical] if self._names is not None else None

    def __repr__(self):
        return f"Electrograms with {self.n_points} mapping points."


class ECG:
    """
    Class for storing information about ECGs

    Args:
        egm (np.ndarray): ECGs
        channel_names (np.ndarray): ECG channel names
        gain (np.ndarray): gain to apply to each signal
        is_electrical: Flag for whether each point has electrical data
    """

    def __init__(
        self,
        ecg: np.ndarray = None,
        channel_names: np.ndarray = None,
        gain: np.ndarray = None,
        is_electrical: np.ndarray = None,
    ):

        # Ensure gain is present is necessary, and that is has the correct shape
        if ecg is not None and gain is None:
            n_points, n_samples, n_channels = ecg.shape
            gain = np.ones((n_points, n_channels))
        elif gain is not None:
            n_points, n_samples, n_channels = ecg.shape
            gain = gain.reshape((n_points, n_channels))

        self._ecg = ecg
        self._channel_names = channel_names
        self._gain = gain

        # Not all points have electrical signals, we need to ignore these.
        # e.g. landmark points in Kodex have no electrical data.
        self._is_electrical = is_electrical

    @property
    def ecg(self):
        return self._ecg[self._is_electrical, :, :] if self._ecg is not None else None

    @property
    def channel_names(self):
        return self._channel_names

    @property
    def n_points(self):
        return self.ecg.shape[0] if self.ecg is not None else 0

    @property
    def n_samples(self):
        return self.ecg.shape[1] if self.ecg is not None else 0

    @property
    def n_channels(self):
        return self.ecg.shape[2] if self.ecg is not None else 0

    @property
    def gain(self):
        return self._gain[self._is_electrical, :] if self._gain is not None else None

    @gain.setter
    def gain(self, gain):
        if isinstance(gain, np.ndarray) and gain.shape[0] == self.n_points:
            self._gain[self._is_electrical] = gain
        else:
            self._gain = gain

    def __repr__(self):
        return f"ECGs with {self.n_points} signals."


@attrs(auto_attribs=True, auto_detect=True)
class Impedance:
    """
    Class for storing information about impedance traces.

    Args:
        times (np.ndarray): Times at which the traces were taken.
        values (np.ndarray): Impedance measurements.
    """

    # TODO: do we need to treat impedances of landmark-only points differently
    times: np.ndarray = None
    values: np.ndarray = None

    def __repr__(self):
        n_traces = len(self.values) if self.values is not None else 0
        return f"Impedance measurements with {n_traces} traces."


class ElectricSurface:
    """
    Class for storing information about electroanatomic data projected onto the surface mesh.

    Args:
        nearest_point (np.ndarray): The position of each mapping point projected to the nearest
            point on the 3D surface constructed by the clinical mapping system.
        normals (np.ndarray): Normal to the surface at each of the nearest points.
        is_electrical (np.ndarray): Flag for whether each point has electrical data.
        """

    def __init__(
        self,
        nearest_point: np.ndarray = None,
        normals: np.ndarray = None,
        is_electrical: np.ndarray = None,
    ):

        self._nearest_point = nearest_point
        self._normals = normals
        self._is_electrical = is_electrical

        if self._nearest_point is None and self._is_electrical is not None:
            self._nearest_point = np.full((is_electrical.size, 3), fill_value=np.NaN, dtype=float)

        if self._normals is None and self._nearest_point is not None:
            self._normals = np.full_like(self._nearest_point, fill_value=np.NaN, dtype=float)
        
    @property
    def nearest_point(self):
        return self._nearest_point[self._is_electrical] if self._nearest_point is not None else None

    @nearest_point.setter
    def nearest_point(self, nearest_point):
        if isinstance(nearest_point, np.ndarray) and nearest_point.shape[0] == self.n_points:
            self._nearest_point[self._is_electrical] = nearest_point
        else:
            self._nearest_point = nearest_point

    @property
    def normals(self):
        return self._normals[self._is_electrical] if self._normals is not None else None

    @normals.setter
    def normals(self, normals):
        if isinstance(normals, np.ndarray) and normals.shape[0] == self.n_points:
            self._normals[self._is_electrical] = normals
        else:
            self._normals = normals

    @property
    def n_points(self):
        return self.nearest_point.shape[0] if self.nearest_point is not None else 0

    def __repr__(self):
        return f"ElectricSurface with {self.n_points} mapping points."


class Annotations:
    """
    Class for storing information about activation times for electrograms.

    Args:
        window_of_interest (np.ndarray): The window of interest for each mapping point
        local_activation_time (np.ndarray): The local activation time for each mapping point
        reference_activation_time (np.ndarray): The reference activation time for each mapping point
        is_electrical (np.ndarray):  Flag for whether each point has electrical data
        frequency (float, optional): Sample frequency (Hz)
    """

    def __init__(
        self,
        window_of_interest: np.ndarray = None,
        local_activation_time: np.ndarray = None,
        reference_activation_time: np.ndarray = None,
        is_electrical: np.ndarray = None,
        frequency: float = 1000,
    ):

        self._window_of_interest_indices = window_of_interest
        self._local_activation_time_indices = local_activation_time
        self._reference_activation_time_indices = reference_activation_time
        self._is_electrical = is_electrical
        self._frequency = frequency

    @property
    def window_of_interest(self):
        try:
            return self._window_of_interest_indices[self._is_electrical] * 1000.0 / self._frequency
        except TypeError as e:
            return None

    @property
    def local_activation_time(self):
        try:
            return self._local_activation_time_indices[self._is_electrical] * 1000.0 / self._frequency
        except TypeError as e:
            return None

    @property
    def reference_activation_time(self):
        try:
            return self._reference_activation_time_indices[self._is_electrical] * 1000.0 / self._frequency
        except TypeError as e:
            return None

    @property
    def n_points(self):
        return len(self.window_of_interest) if self.window_of_interest is not None else 0

    def __repr__(self):
        return f"Annotations with {self.n_points} mapping points."


class Electric:
    """
    Class for storing electrical data obtained during a clinical mapping procedure.

    Args:
        names (np.ndarray): The physician-visible names of the points applied during the
            clinical case.
        internal_names (np.ndarray): The internal names of the points used by the clinical
            electroanatomic mapping system.
        include (np.ndarray): Flag for whether each point should be used for interpolating
            data onto a surface.
        is_electrical (np.ndarray): Flag for whether each point has electrical data.
            Not all points have electrical signals (e.g. landmark points in Kodex), and we need
            to ignore these when e.g. interpolating, plotting.
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
        frequency (float): The sample frequency of electric signals, in Hz. Defaults to 1000 Hz.
    """
    def __init__(
        self,
        names: np.ndarray = None,
        internal_names: np.ndarray = None,
        include: np.ndarray = None,
        is_electrical: np.ndarray = None,
        bipolar_egm: Electrogram = None,
        unipolar_egm: Electrogram = None,
        reference_egm: Electrogram = None,
        ecg: np.ndarray = None,
        impedance: Impedance = None,
        surface: ElectricSurface = None,
        annotations: Annotations = None,
        frequency: float = 1000,
    ):

        self._names = names
        self._internal_names = internal_names
        self._include = include
        self._is_electrical = is_electrical
        self._is_electrical_indices = np.nonzero(is_electrical)[0].ravel()
        self.bipolar_egm = bipolar_egm
        self.unipolar_egm = unipolar_egm
        self.reference_egm = reference_egm
        self.ecg = ecg
        self.impedance = impedance
        self.surface = surface
        self.annotations = annotations
        self.frequency = frequency

        if bipolar_egm is None:
            self.bipolar_egm = Electrogram()

        if self.unipolar_egm is None:
            self.unipolar_egm = Electrogram()

        if self.reference_egm is None:
            self.reference_egm = Electrogram()

        # We need to extract coords of landmark points
        self._is_landmark = None
        if self._names is not None:
            self._is_landmark = np.full_like(
                self._names,
                fill_value=False,
                dtype=bool,
            )
            self._is_landmark[self._names!=''] = True

        self.landmark_points = LandmarkPoints(
            points=self.bipolar_egm._points,
            names=self._names,
            internal_names=self._internal_names,
            is_landmark=self._is_landmark,
        )

        if self.ecg is None:
            self.ecg = ECG()

        if self.impedance is None:
            self.impedance = Impedance()

        if self.surface is None:
            self.surface = ElectricSurface()

        if self.annotations is None:
            self.annotations = Annotations()

        self._time_indices = np.arange(self.n_samples)

    @property
    def names(self):
        return self._names[self._is_electrical] if self._names is not None else None

    @property
    def internal_names(self):
        return self._internal_names[self._is_electrical] if self._internal_names is not None else None

    @property
    def include(self):
        return self._include[self._is_electrical] if self._include is not None else None

    @include.setter
    def include(self, include):
        if isinstance(include, np.ndarray) and include.size == self.n_points:
            self._include[self._is_electrical] = include
        else:
            self._include = include

    @property
    def times(self):
        return self._time_indices * 1000 / self.frequency  # time in ms

    @property
    def n_points(self):
        return self.bipolar_egm.n_points

    @property
    def n_samples(self):
        return self.bipolar_egm.n_samples

    def __repr__(self):
        return f"Electric data for {self.n_points} mapping points."


    def _add_landmark(
        self,
        name: str,
        internal_name: str,
        point: np.ndarray,
    ):
        """Add a landmark point."""

        # We need to add rows to **all** all signals
        # Not necessary for openep-py, but it is for openep-matlab

        if name == '':
            raise ValueError("name cannot be an empty string.")
        if internal_name == '':
            raise ValueError("internal_name cannot be an empty string.")

        if isinstance(internal_name, str):
            internal_name = np.array([internal_name], dtype=str)
        if isinstance(name, str):
            name = np.array([name], dtype=str)

        is_electrical = np.array([False], dtype=bool)
        is_landmark = np.array([True], dtype=bool)
        include = np.array([False], dtype=int)

        if self._names is None:
            self._names = name
        else:
            self._names = np.hstack([self._names, name])

        if self._internal_names is None:
            self._internal_names = internal_name
        else:
            self._internal_names = np.hstack([self._internal_names, internal_name])

        if self._is_electrical is None:
            self._is_electrical = is_electrical
        else:
            self._is_electrical = np.hstack([self._is_electrical, is_electrical])

        self._is_electrical_indices = np.nonzero(self._is_electrical)[0].ravel()

        if self._is_landmark is None:
            self._is_landmark = is_landmark
        else:
            self._is_landmark = np.hstack([self._is_landmark, is_landmark])

        if self._include is None:
            self._include = include
        else:
            self._include = np.hstack([self._include, include])

        point = np.asarray(point)
        if point.ndim == 1:
            point = point[np.newaxis, :]

        # set n_samples to 1 if we have no electrical data
        n_samples = self.bipolar_egm.n_samples or self.unipolar_egm.n_samples or self.ecg.n_samples or 1
        egm = np.full(n_samples, fill_value=np.NaN, dtype=float)

        # We need to create bipolar egm data if it does not exist
        # This is because openep-matlab stores landmark data with the bipolar data
        if self.bipolar_egm._points is None:

            self.bipolar_egm = Electrogram(
                egm=egm,
                points=point,
                is_electrical=self._is_electrical,
            )

        else:
            self.bipolar_egm._egm = np.vstack([self.bipolar_egm._egm, egm])
            self.bipolar_egm._points = np.vstack([self.bipolar_egm._points, point])
            self.bipolar_egm._voltage = np.hstack([self.bipolar_egm._voltage, [np.NaN]])
            self.bipolar_egm._gain = np.hstack([self.bipolar_egm._gain, [1]])
            self.bipolar_egm._names = np.hstack([self.bipolar_egm._names, ['']])
            self.bipolar_egm._is_electrical = self._is_electrical

        # Add all landmarks
        self.landmark_points = LandmarkPoints(
            points=self.bipolar_egm._points,
            names=self._names,
            internal_names=self._internal_names,
            is_landmark=self._is_landmark,
        )

        # Update unipolar data if necessary
        if self.unipolar_egm._egm is not None:
            
            unipolar_egm = np.full((1, n_samples, 2), fill_value=np.NaN, dtype=float)
            unipolar_point = np.full((1, 3, 2), fill_value=np.NaN, dtype=float)
            unipolar_gain = np.zeros((1, 2))
            unipolar_name = np.full((1, 2), fill_value=' ', dtype=str)

            self.unipolar_egm._egm = np.vstack([self.unipolar_egm._egm, unipolar_egm])
            self.unipolar_egm._points = np.vstack([self.unipolar_egm._points, unipolar_point])
            self.unipolar_egm._voltage = np.hstack([self.unipolar_egm._voltage, [np.NaN]])
            self.unipolar_egm._gain = np.vstack([self.unipolar_egm._gain, unipolar_gain])
            self.unipolar_egm._names = np.vstack([self.unipolar_egm._names, unipolar_name])
            self.unipolar_egm._is_electrical = self._is_electrical

        # Update reference data if necessary
        if self.reference_egm._egm is not None:
            self.reference_egm._egm = np.vstack([self.reference_egm._egm, egm])
            self.reference_egm._voltage = np.hstack([self.reference_egm._voltage, [np.NaN]])
            self.reference_egm._gain = np.hstack([self.reference_egm._gain, [1]])
            self.reference_egm._names = np.hstack([self.reference_egm._names, ['']])
            self.reference_egm._is_electrical = self._is_electrical

        # Update ecg data if necessary
        if self.ecg._ecg is not None:
            ecg = np.full((1, self.ecg.n_samples, self.ecg.n_channels), fill_value=np.NaN, dtype=float)
            ecg_gain = np.ones((1, self.ecg.n_channels), dtype=float)

            self.ecg._ecg = np.vstack([self.ecg._ecg, ecg])
            self.ecg._gain = np.vstack([self.ecg._gain, ecg_gain])
            self.ecg._is_electrical = self._is_electrical

        # Update annotations if necesary
        if self.annotations._window_of_interest_indices is not None:

            woi = np.array([[-8000, 8000]], dtype=int)
            reference_activation_time = [0]
            local_activation_time = [0]

            self.annotations._window_of_interest_indices = np.vstack([
                self.annotations._window_of_interest_indices, woi,
            ])
            self.annotations._reference_activation_time_indices = np.hstack([
                self.annotations._reference_activation_time_indices, reference_activation_time,
            ])
            self.annotations._local_activation_time_indices = np.hstack([
                self.annotations._local_activation_time_indices, local_activation_time,
            ])
            self.annotations._is_electrical = self._is_electrical


def extract_electric_data(electric_data):
    """Extract electric data from a dictionary.

    Args:
        electric_data (dict): Dictionary containing numpy arrays that describe the
            electric data associated with electrograms taken at various mapping points.

    Returns:
        electric (Electric): object containing electric data associated with electrograms
            taken at various mapping points.
    """

    if electric_data['egm'].size == 0 and electric_data['egmUni'].size == 0:
        return empty_electric()

    names = electric_data['tags'].astype(str)
    internal_names = electric_data['names'].astype(str)

    # electric.include is a later addition to the openep format
    if 'include' in electric_data:
        include = electric_data['include'].astype(int)
    else:
        include = np.full_like(
            names,
            fill_value=True,
            dtype=int,
        )

    # We need to know which points are landmark only and have no electrical data.
    # Those that have NaN values for 
    is_electrical = ~np.all(np.isnan(electric_data['egm'].astype(float)), axis=1)

    # Older versions of OpenEP datasets did not have unipolar data or electrode names. Add deafult ones here.
    if 'electrodeNames_bip' not in electric_data:
        electric_data['electrodeNames_bip'] = np.full_like(internal_names, fill_value=" ", dtype=str)
    if 'egmUni' not in electric_data:
        electric_data['egmUni'] = None
        electric_data['egmUniX'] = None
        electric_data['voltages']['unipolar'] = None
        electric_data['electrodeNames_uni'] = None
    else:
        electric_data['egmUni'] = electric_data['egmUni'].astype(float)
        electric_data['egmUniX'] = electric_data['egmUniX'].astype(float)
        electric_data['voltages']['unipolar'] = electric_data['voltages']['unipolar'].astype(float)
        electric_data['electrodeNames_uni'] = electric_data['electrodeNames_uni'].astype(str)

    # Make ecgs correct shape
    ecg_dims = electric_data['ecg'].ndim
    if ecg_dims == 1:
        n_ecg_channels = 0
    elif ecg_dims == 2:
        n_ecg_channels = 1
        n_points, n_samples = electric_data['ecg'].shape
        electric_data['ecg'] = electric_data['ecg'].reshape(n_points, n_samples, 1)
    elif ecg_dims == 3:
        n_ecg_channels = electric_data['ecg'].shape[2]

    # Not all datasets have ecgNames
    electric_data['ecgNames'] = electric_data.get('ecgNames', np.array(['ECG' for _ in range(n_ecg_channels)]))
    if isinstance(electric_data['ecgNames'], str):
        electric_data['ecgNames'] = np.array([electric_data['ecgNames']]).astype(str)
    if electric_data['ecgNames'].size == 0:
        electric_data['ecg'] = None
        electric_data['ecgNames'] = None
    else:
        electric_data['ecg'] = electric_data['ecg'].astype(float)

    # Not all datasets have gain values
    egm_types = ['', 'Ref', 'Uni']  # bipolar, reference, unipolar
    default_gain_values = [1.0, -4.0, 0.0]
    for egm_type, default_gain_value in zip(egm_types, default_gain_values):

        if electric_data[f'egm{egm_type}'].size == 0:
            electric_data[f'egm{egm_type}Gain'] = None
        elif f'egm{egm_type}Gain' not in electric_data:
            electric_data[f'egm{egm_type}Gain'] = np.full(
                electric_data[f'egm{egm_type}'].shape[0],
                fill_value=default_gain_value,
                dtype=float,
            )
        else:
            electric_data[f'egm{egm_type}Gain'] = electric_data[f'egm{egm_type}Gain'].astype(float)

    # we need a gain value for each electrode of the unipolar signals
    if electric_data['egmUniGain'] is not None and electric_data['egmUniGain'].ndim == 1:
        electric_data['egmUniGain'] = np.tile(electric_data['egmUniGain'], reps=(2, 1)).T

    if 'ecgGain' not in electric_data or electric_data['ecgGain'].size == 0:
        electric_data['ecgGain'] = None
    else:
        electric_data['ecgGain'] = electric_data['ecgGain'].astype(float)

    # Create objects to pass to Electric
    bipolar_egm = Electrogram(
        egm=electric_data['egm'].astype(float),
        points=electric_data['egmX'].astype(float),
        voltage=electric_data['voltages']['bipolar'].astype(float),
        gain=electric_data['egmGain'],
        names=electric_data['electrodeNames_bip'].astype(str) if electric_data['electrodeNames_bip'].size > 0 else None,
        is_electrical=is_electrical,
    )
    unipolar_egm = Electrogram(
        egm=electric_data['egmUni'] if electric_data['egmUni'].size > 0 else None,
        points=electric_data['egmUniX'] if electric_data['egmUniX'].size > 0 else None,
        voltage=electric_data['voltages']['unipolar'] if electric_data['voltages']['unipolar'].size > 0 else None,
        gain=electric_data['egmUniGain'],
        names=electric_data['electrodeNames_uni'] if electric_data['electrodeNames_uni'].size > 0 else None,
        is_electrical=is_electrical if electric_data['egmUni'].size > 0 else None,
    )
    reference_egm = Electrogram(
        egm=electric_data['egmRef'].astype(float) if electric_data['egmRef'].size > 0 else None,
        gain=electric_data['egmRefGain'],
        is_electrical=is_electrical if electric_data['egmRef'].size > 0 else None,
    )

    ecg = ECG(
        ecg=electric_data['ecg'],
        channel_names=electric_data['ecgNames'],
        gain=electric_data['ecgGain'],
        is_electrical=is_electrical if electric_data['ecg'] is not None and electric_data['ecg'].size > 0 else None,
    )

    try:
        impedance_times = electric_data['impedances']['time'].astype(float)
        impedance_values = electric_data['impedances']['value'].astype(float)
    except AttributeError as e:
        # We have a list of arrays, rather than a single array
        impedance_times = electric_data['impedances']['time']
        impedance_values = electric_data['impedances']['value']

    impedance = Impedance(
        times=impedance_times if len(impedance_times) > 0 else None,
        values=impedance_values if len(impedance_values) > 0 else None,
    )

    surface = ElectricSurface(
        nearest_point=electric_data['egmSurfX'].astype(float) if electric_data['egmSurfX'].astype(float).size > 0 else None,
        normals=electric_data['barDirection'].astype(float) if electric_data['barDirection'].astype(float).size > 0 else None,
        is_electrical=is_electrical if electric_data['egmSurfX'].size > 0 else None,
    )

    # If no sample frequency is specified, assume it's 1000 Hz
    frequency = electric_data.get('sampleFrequency', 1000.0)
    try:
        frequency = float(frequency)
    except TypeError as e:
        frequency = 1000.0

    annotations = Annotations(
        window_of_interest=electric_data['annotations']['woi'].astype(int),
        local_activation_time=electric_data['annotations']['mapAnnot'].astype(int),
        reference_activation_time=electric_data['annotations']['referenceAnnot'].astype(int),
        frequency=frequency,
        is_electrical=is_electrical,
    )

    electric = Electric(
        names=names,
        internal_names=internal_names,
        include=include,
        is_electrical=is_electrical,
        bipolar_egm=bipolar_egm,
        unipolar_egm=unipolar_egm,
        reference_egm=reference_egm,
        ecg=ecg,
        impedance=impedance,
        surface=surface,
        annotations=annotations,
        frequency=frequency,
    )

    return electric


def empty_electric():
    """Create an empty Ablation object with empty numpy arrays.

    Returns:
        electric (Electric): object containing electric data associated with electrograms
            taken at each mapping point.
    """

    names = np.array([], dtype=str)
    internal_names = np.array([], dtype=str)
    include = np.array([], dtype=bool)
    is_electrical = np.array([], dtype=bool)

    bipolar_egm = Electrogram(
        egm=np.array([], dtype=float),
        points=np.array([], dtype=float),
        voltage=np.array([], dtype=float),
        gain=np.array([], dtype=float),
        names=np.array([], dtype=str),
        is_electrical=np.array([], dtype=bool),
    )
    unipolar_egm = Electrogram(
        egm=np.array([], dtype=float),
        points=np.array([], dtype=float),
        voltage=np.array([], dtype=float),
        gain=np.array([], dtype=float),
        names=np.array([], dtype=str),
        is_electrical=np.array([], dtype=bool),
    )
    reference_egm = Electrogram(
        egm=np.array([], dtype=float),
        gain=np.array([], dtype=float),
        is_electrical=np.array([], dtype=bool),
    )

    ecg = ECG(
        ecg=np.array([], dtype=float),
        channel_names=np.array([], dtype=str),
        gain=np.array([], dtype=float),
        is_electrical=np.array([], dtype=bool),
    )

    impedance = Impedance(
        times=np.array([], dtype=float),
        values=np.array([], dtype=float),
    )

    surface = ElectricSurface(
        nearest_point=np.array([], dtype=float),
        normals=np.array([], dtype=float),
        is_electrical=np.array([], dtype=bool),
    )

    frequency = np.NaN

    annotations = Annotations(
        window_of_interest=np.array([], dtype=float),
        local_activation_time=np.array([], dtype=float),
        reference_activation_time=np.array([], dtype=float),
        frequency=frequency,
        is_electrical=np.array([], dtype=bool),
    )

    electric = Electric(
        names=names,
        internal_names=internal_names,
        include=include,
        is_electrical=is_electrical,
        bipolar_egm=bipolar_egm,
        unipolar_egm=unipolar_egm,
        reference_egm=reference_egm,
        ecg=ecg,
        impedance=impedance,
        surface=surface,
        annotations=annotations,
        frequency=frequency,
    )

    return electric
