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

import sys
import time
import operator
import functools
import numpy as np
import matplotlib.pyplot as plt

import scipy.spatial
from scipy.interpolate import Rbf
from scipy.interpolate import LinearNDInterpolator as linterp
from scipy.interpolate import NearestNDInterpolator as nearest

__all__ = [
    'get_mapping_points_within_woi',
    'get_woi_times',
    'get_electrograms_at_points',
    'plot_electrograms',
    'calculate_distance',
    'calculate_points_within_distance',
    'LinearNDInterpolatorExt',
    'LocalSmoothing',
    'rbfAssemble',
    'rbfcreate',
    'rbfinterp',
    'rbfcheck',
    'rbf_scipy',
    'OpenEPDataInterpolator',
    'get_voltage_electroanatomic',
]


def _get_reference_annotation(case, indices=None):
    """
    Get the reference annotations for mapping points.

    Args:
        case (Case): openep case object
        indices (ndarray), optional: indices of reference annotations to select. The default
            is None, in which case all reference annotation will be returned.
    Returns:
        annotations (ndarray): reference annotations
    """

    annotations = case.electric['annotations/referenceAnnot'].ravel()
    annotations = annotations[indices] if indices is not None else annotations
    
    return annotations


def _get_window_of_interest(case, indices=None):

    """
    Gets the window of interest for mapping points.
    
    Args:
        case (Case): openep case object
        indices (ndarray), optional: indices of mapping points for which the woi will
            be extracted. The default is None, in which case the woi will be returned for
            every mapping point.

    Returns:
        woi (ndarray): Nx2 array with the windows of interest for selected mapping points.
            For each point, the two columns specify that start and end time respectively
            of the window of interest.
    """
    
    woi = case.electric["annotations/woi"].T
    woi = woi[indices] if indices is not None else woi.copy()

    return woi


def get_mapping_points_within_woi(case, indices=None, buffer=50):
    """
    Extracts the indices of the mapping points that have
    annotated local activation times within the window of interest (woi).

    Args:
        case (Case): openep case object
        indices (ndarray), optional: indices of mapping points to consider. The default
            is None, in which case all points within the window of interest will be returned.
        buffer (float): points within the woi plus/minus this buffer time will
            be considered to be within the woi.

    Returns:
        within_woi (ndarray): boolean array of valid points; True if
            the annotated local activation time is within the woi, False otherwise.
    """

    reference_annot = _get_reference_annotation(case, indices=indices)
    woi = _get_window_of_interest(case, indices=indices)
    woi += reference_annot[:, np.newaxis]

    map_annot = case.electric["annotations/mapAnnot"].ravel()
    map_annot = map_annot[indices] if indices is not None else map_annot
    
    within_woi = np.logical_and(
        map_annot >= woi[:, 0] - buffer,
        map_annot <= woi[:, 1] + buffer,
    )

    return within_woi


def get_electrograms_at_points(
    case,
    woi=True,
    buffer=50,
    indices=None,
    bipolar=True,
    return_names=True,
    return_lat=True,
):
    """
    Extract the electrogram timeseries for a selection of points.

    Args:
        case (Case): openep case object
        woi (bool): If True, only electrograms within the window of interest will be extracted.
            If False, all electrograms will be extracted.
        buffer (float): If woi is True, points within the woi plus/minus this buffer time will
            be considered to be within the woi. If woi is False, buffer is ignored.
        indices (ndarray), optional: indices of mapping points for which electrograms will
            be extracted. If provided along with `woi=True`, only the electrograms that
            are both within the window of interest and selected by `indices` will be extracted.
            If provided along with `woi=False`, all electrograms of mapping points selected by
            `indices` will be returned.
        bipolar (bool): If True, the bipolar traces will be returned. If False, the unipolar
            traces will be returned.
        return_names (bool): If True, the internal names of the points used by the
            clinical electroanatomic mapping system will also be returned.
        return_lat (bool): If True, the local activation time of each mapping point will also be
            returned.

    Returns:
        traces (ndarray): A timeseries of voltages for each selected mapping point.
        names (ndarray): If `return_names` is True, the internal names of the points used by the
            clinical electroanatomic mapping system will be returned.
            local_activation_time (ndarray): If `return_lat` is True, the local activation time of each
            mapping point will be returned.
    """
    
    electrograms = case.electric['egm'].T if bipolar else case.electric['egmUni'].T
    names = case.electric["names"]
    local_activation_time = case.electric['annotations/mapAnnot'].ravel()

    # Filter by selected indices
    if indices is not None:

        # if we have a single index we need to ensure it is an array
        indices = np.asarray([indices], dtype=int) if isinstance(indices, int) else indices
        print(indices)
        electrograms = electrograms[indices]
        names = names[indices]
        local_activation_time = local_activation_time[indices]

    # Filter by window of interest and buffer
    if woi:
        within_woi = get_mapping_points_within_woi(case, indices=indices, buffer=buffer)
        electrograms = electrograms[within_woi]
        names = names[within_woi]
        local_activation_time = local_activation_time[within_woi]

    if return_names and return_lat:
        return electrograms, names, local_activation_time
    elif return_names:
        return electrograms, names
    elif return_lat:
        return electrograms, local_activation_time
    
    return electrograms


def get_woi_times(case, buffer=50, relative=False):
    """
    Extracts the times within the window of interest.

    Args:
        case (Case): openep case object
        buffer (float): times within the window of interest plus/minus this buffer
            time will be considered to be within the woi.
        relative (bool): if True, the time of the reference annotation will be
            subtracted from the returned times.

    Returns:
        times (ndarray): times within the window of interest
    """
    
    times = np.arange(case.electric['egm'].T.shape[1])
    
    woi = case.electric['annotations/woi'].T[0]
    ref_annotation = int(case.electric['annotations/referenceAnnot'].T[0, 0])
    buffer = 50

    start_time, stop_time = woi + ref_annotation + [-buffer, buffer]

    keep_times = np.logical_and(
        times >= start_time,
        times <= stop_time,
    )

    if relative:
        times -= ref_annotation

    return times[keep_times]


def plot_electrograms(times, electrograms, separation=0.5, names=None, axis=None):
    """
    Plot electrogram traces.

    Args:
        times (ndarray): times at which voltages were measured
        electrograms (ndarray): Electrogram traces. Two-dimensional of size N_points x N_times for bipolar voltage,
            or two-dimensional of shape N_points x N_times x 2 for unipolar dimensional.
        woi (bool): If True, the traces will be plotted only within the window of interest.
        buffer (float): If woi is True, points within the woi plus/minus this buffer time will
            be considered to be within the woi. If woi is False, buffer is ignored.
        axis (matplotlib.axes.Axes): Matplotlib Axes on which to plot the traces. If None, a new figure and axes
            will be created.

    Returns:
        axis (matplotlib.axes.Axes): Axes on which the traces have been plotted.
    """
    
    # TODO: coloring of the lines - make it sensible, both for bipolar and unipolar
    # TODO: add title and labels? Don't make the y-axis invisible?
    # TODO: what default size settings should be used for the figure?
    # TODO: add a horizontal line showing the zero-value voltage for each trace
    # TODO: add a vertical line showing time zero
    
    separations = np.arange(electrograms.shape[0]) * separation
    
    if axis is None:
        _, axis = plt.subplots()
    
    plt.sca(axis)
    if electrograms.ndim == 2:  # bipolar voltage
        plt.plot(times, electrograms.T + separations, label=names)
    else:  # unipolar voltages
        plt.plot(times, electrograms[:, :, 0].T + separations, label=names)
        plt.plot(times, electrograms[:, :, 1].T + separations, label=names)

    axis.get_yaxis().set_visible(False)
    
    return axis.get_figure(), axis


def calculate_distance(origin, destination):
    """
    Returns the distance from a set of origin points to a set of destination
    points.

    Args:
        origin (ndarray): Nx3 matrix of coordinates
        destination (ndarray): Mx3 matrix of coordinates

    Returns:
        distances (ndarray): MxN matrix of distances

    """
    
    origin = origin[np.newaxis, :] if origin.ndim == 1 else origin
    destination = destination[np.newaxis, :] if destination.ndim == 1 else destination
    
    distances = scipy.spatial.distance.cdist(
        origin,
        destination,
    ).T
    
    return distances


def calculate_points_within_distance(origin, destination, max_distance, return_distances=True):
    """
    Calculates whether the distances from a set of origin points to a set of
    destination points are each within a specified distance cutoff.

    Args:
        points (ndarray, (N, 3)): array of 3D points
        test_points (ndarray, (M, 3)): array of 3D test points
        max_distance (float): distance threshold between the origin points and
            destination points

    Returns:
        within_max_dist (ndarray, M x N): Boolean array
            that is equal to True if the points are within the maximum distance
            of one another and equal to False otherwise. 
        distances (ndarray, M x N): distance between each point
            and each test_point.
    """

    distances = calculate_distance(origin, destination)
    within_max_distance = distances <= max_distance
    
    if return_distances:
        return within_max_distance, distances
    
    return within_max_distance


class LinearNDInterpolatorExt(object):
    """
    Interpolation Method - Python implementation of ScatteredInterpolation
    """

    def __init__(self, points, values):
        self.funcinterp = linterp(points, values)
        self.funcnearest = nearest(points, values)

    def __call__(self, *args):
        z = self.funcinterp(*args)
        chk = np.isnan(z)
        if chk.any():
            return np.where(chk, self.funcnearest(*args), z)
        else:
            return z


def LocalSmoothing(x0, x1, smoothingLength):
    f_dash = np.zeros(shape=(len(x1), 1), dtype=np.float64)
    df_dash = np.zeros(shape=(len(x1), 3), dtype=np.float64)
    [idx, dists] = calculate_points_within_distance(x0, x1, smoothingLength)

    return [f_dash, df_dash]


def rbfAssemble(x, phi, const, smooth):
    xShape = x.shape
    if len(xShape) == 2:
        dim = xShape[0]
        n = xShape[1]
    else:
        dim = 1
        n = xShape[0]
    A = np.zeros([n, n])
    for i in range(n):
        for j in range(i + 1):
            r = np.linalg.norm(x[:, i] - x[:, j])
            temp, _ = phi(r, const)
            A[i, j] = temp
            A[j, i] = temp
        A[i, i] = A[i, i] - smooth
    # Polynomial part
    P = np.c_[np.ones(x.shape[1]), x.T]
    A = np.c_[A, P]
    A = np.r_[A, np.c_[P.T, np.zeros([dim + 1, dim + 1])]]
    return A


class rbfcreate:
    def __init__(self, x, y, *args):
        start_time = time.time()
        nargin = 2 + len(args)

        Names = [
            "RBFFunction      ",
            "RBFConstant      ",
            "RBFSmooth        ",
            "Stats            ",
        ]
        m = len(Names)
        names = []
        for i in range(m):
            names.append(Names[i].lower())
        # Check input arrays
        xShape = x.shape
        if len(xShape) == 2:
            nXDim = xShape[0]
            nXCount = xShape[1]
        else:
            nXDim = 1
            nXCount = xShape[0]
        yShape = y.shape
        if len(yShape) == 2:
            nYDim = yShape[0]
            nYCount = yShape[1]
        else:
            nYDim = 1
            nYCount = yShape[0]

        if nXCount != nYCount:
            sys.exit("x and y should have the same number of rows")

        if nYDim != 1:
            sys.exit("y should be n by 1 vector")

        self.x = x
        self.y = y

        # Default values
        self.RBFFunction = "linear"
        self.RBFConstant = (
            functools.reduce(
                operator.mul, (np.amax(x.T, axis=0) - np.amin(x.T, axis=0)), 1
            )
            / nXCount
        ) ** (
            1 / nXDim
        )  # approx. average distance between the nodes
        self.RBFSmooth = 0
        self.Stats = "off"

        i = 0
        if np.remainder(nargin - 2, 2) != 0:
            sys.exit("Arguments must occur in name-value pairs.")

        expectval = 0  # start expecting a name, not a value

        while i < nargin - 2:
            arg = args[i]

            if not expectval:
                if not isinstance(arg, str):
                    sys.exit(
                        "Expected argument " + i + " to be a string property name."
                    )

                lowArg = arg.lower()
                j = []
                for p in range(m):
                    if names[p].startswith(lowArg):
                        j.append(p)
                if not j:  # if no matches
                    sys.exit("Unrecognized property name " + arg)
                elif len(j) > 1:  # if more than one match
                    # Check for any exact matches (in case any names are subsets of others)
                    k = []
                    for p in range(m):
                        if lowArg == names[p]:
                            k.append(1)
                    if len(k) == 1:
                        j = k
                    else:
                        msg = (
                            "Ambiguous property name "
                            + arg
                            + "("
                            + Names[j[0]].replace(" ", "")
                        )
                        for p in range(1, len(k)):
                            msg = msg + ", " + Names[j[p]].replace(" ", "")
                        msg = msg + ")."
                        sys.exit(msg)
                expectval = 1  # we expect a value next
            else:
                if Names[j[0]].replace(" ", "") == "RBFFunction":
                    self.RBFFunction = arg
                elif Names[j[0]].replace(" ", "") == "RBFConstant":
                    self.RBFConstant = arg
                elif Names[j[0]].replace(" ", "") == "RBFSmooth":
                    self.RBFSmooth = arg
                elif Names[j[0]].replace(" ", "") == "Stats":
                    self.stats = arg
                expectval = 0
            i += 1

        if expectval:
            sys.exit("Expected value for property " + arg)

        # **************************************************************************
        # Creating RBF Interpolatin
        # **************************************************************************

        if self.RBFFunction == "linear":
            self.rbfphi = lambda r, const: (r, 1 / r)
        elif self.RBFFunction == "cubic":
            self.rbfphi = lambda r, const: (r * r * r, 3 * r)
        elif self.RBFFunction == "multiquadric":
            self.rbfphi = lambda r, const: (
                np.sqrt(1 + r * r / (const * const)),
                1 / (const * const * np.sqrt(1 + r * r / (const * const))),
            )
        elif self.RBFFunction == "thinplate":
            self.rbfphi = lambda r, const: (
                r * r * np.log(r + 1),
                r / ((r + 1) + 2 * np.log(r + 1)),
            )
        elif self.RBFFunction == "gaussian":
            self.rbfphi = lambda r, const: (
                np.exp(-0.5 * r * r / (const * const)),
                -np.exp(-0.5 * r * r / (const * const)) / (const * const),
            )
        else:
            self.rbfphi = lambda r, const: (r, 1 / r)

        phi = self.rbfphi

        A = rbfAssemble(x, phi, self.RBFConstant, self.RBFSmooth)

        b = np.r_[np.reshape(y, [nYCount, 1], order="F"), np.zeros([nXDim + 1, 1])]

        # inverse
        rbfcoeff = np.linalg.lstsq(A, b, rcond=None)[0]
        self.rbfcoeff = rbfcoeff

        if self.Stats == "on":
            print(
                len(y)
                + "point RBF interpolation was created in "
                + time.time()
                - start_time
                + " sec."
            )


def rbfinterp(x, options):

    start_time = time.time()

    phi = options.rbfphi
    rbfconst = options.RBFConstant
    nodes = options.x
    rbfcoeff = options.rbfcoeff

    dim, n = nodes.shape
    dimPoints, nPoints = x.shape

    if dim != dimPoints:
        sys.exit(
            "x should have the same number of rows as an array used to create RBF interpolation"
        )

    f = np.zeros([nPoints])
    df = np.zeros([dimPoints, nPoints])

    xDim, _ = x.shape

    for i in range(nPoints):
        ds = np.zeros([dimPoints, 1])
        dx = (
            np.matmul(np.reshape(x[:, i], [xDim, 1], order="F"), np.ones([1, n]))
            - nodes
        )
        r = np.sqrt(np.sum(dx * dx, axis=0))

        rbf, d_rbf = phi(r, rbfconst)

        # First term is constant coefficient from monomial terms
        s = rbfcoeff[n] + np.sum(rbfcoeff[0:n] * np.reshape(rbf, [n, 1], order="F"))
        for k in range(dim):
            s = s + rbfcoeff[k + n + 1] * x[k, i]  # add monomial terms
            # Second term here are coefficients from monomial terms
            ds[k] = (
                np.sum(
                    rbfcoeff[0:n]
                    * np.reshape(dx[k, :], [n, 1], order="F")
                    * np.reshape(d_rbf, [n, 1], order="F")
                )
                + rbfcoeff[k + n + 1]
            )
        f[i] = s
        df[:, i] = np.reshape(ds, [1, dimPoints], order="F")

    if options.Stats == "on":
        print(
            "Interpolation at "
            + len(f)
            + " points was computed in "
            + time.time()
            - start_time
            + " sec."
        )
    return f, df


def rbfcheck(options):
    start_time = time.time()

    nodes = options.x
    y = options.y

    s, _ = rbfinterp(nodes, options)

    print("RBF Check")
    print("max|y - yi| = " + str(np.max(abs(s - y))))

    if options.Stats == "on":
        print(len(y) + " points were checked in " + time.time() - start_time + " sec.")


def rbf_scipy(x0, d0, x1, rbfConstant):

    x = x0[:, 0]
    y = x0[:, 1]
    z = x0[:, 2]
    d = d0

    F = Rbf(x, y, z, d, function="multiquadric", smooth=rbfConstant)
    d1 = F(x1[:, 0], x1[:, 1], x1[:, 2])

    return d1


class OpenEPDataInterpolator:
    """
    Creates objects for performing spatial interpolation for OpenEP data
    """

    def __init__(self, method, distanceThreshold, rbfConstant):
        self.method = method
        self.distanceThreshold = distanceThreshold
        self.rbfConstant = rbfConstant

    def interpolate(self, x0, d0, x1, *args):
        self.x0 = x0
        self.d0 = d0
        self.x1 = x1

        if self.method == "scatteredinterpolant":
            F = LinearNDInterpolatorExt(
                points=(self.x0[:, 0], self.x0[:, 1], self.x0[:, 2]), values=self.d0
            )
            d1 = np.array(F(self.x1[:, 0], self.x1[:, 1], self.x1[:, 2])).reshape(
                self.x1.shape[0], 1
            )

        elif self.method == "rbf":
            op = rbfcreate(
                self.x0.T,
                self.d0.T,
                "RBFFunction",
                "multiquadric",
                "RBFConstant",
                self.rbfConstant,
            )
            rbfcheck(op)
            out = rbfinterp(self.x1.T, op)
            d1 = np.array(out[0]).reshape(self.x1.shape[0], 1)

        elif self.method == "localSmoothing":
            pass

        else:
            print("Interpolation Method not Found")

        # Workout if there are any points on the surface that are < 0
        d1[d1 < 0] = 0

        #  work out which points on the surface are too far away from real data
        # Remove any interpolated values which are outwith the fill threshold

        from sklearn.neighbors import NearestNeighbors

        neigh = NearestNeighbors(
            n_neighbors=1,
            radius=1,
            algorithm="auto",
            leaf_size=30,
            metric="minkowski",
            p=2,
        )

        neigh.fit(self.x0)

        id = neigh.kneighbors(self.x1, return_distance=False)
        cPts = self.x0[id]
        cPts = np.array(cPts[:, 0])

        # Is this correct? It's only calculating row-wise distance, not all combinaions
        d = calculate_distance(cPts, self.x1).diagonal()[:, np.newaxis]

        thresholdDistance = np.zeros(shape=d.shape, dtype=np.bool)
        thresholdDistance[d > self.distanceThreshold] = 1
        d1[thresholdDistance] = "nan"
        d1 = d1.flatten()

        return d1


def get_voltage_electroanatomic(mesh_case):
    distance_thresh = 10
    rbf_constant_value = 1

    # Anatomic descriptions (Mesh) - nodes and indices
    pts = mesh_case.nodes

    # Electric data
    # Locations – Cartesian co-ordinates, projected on to the surface
    locations = mesh_case.electric['egmX'].T

    i_egm = mesh_case.electric["egm"].T
    i_vp = get_mapping_points_within_woi(mesh_case)[:, np.newaxis]
    # macthing the shape of ivp with data
    i_vp_egm = np.repeat(i_vp, repeats=i_egm.shape[1], axis=1)
    # macthing the shape of ivp with coords
    i_vp_locations = np.repeat(i_vp, repeats=locations.shape[1], axis=1)

    # Replacing the values outside the window of interest with Nan values
    i_egm[~i_vp_egm] = np.nan
    locations[~i_vp_locations] = np.nan

    # For each mapping point, n, find the voltage amplitude
    max_volt = np.amax(a=i_egm, axis=1).reshape(len(i_egm), 1)
    min_volt = np.amin(a=i_egm, axis=1).reshape(len(i_egm), 1)

    amplitude_volt = np.subtract(max_volt, min_volt)

    for indx in range(amplitude_volt.shape[1]):
        temp_data = amplitude_volt[:, indx]
        temp_coords = locations
        i_nan = np.isnan(temp_data)
        temp_data = temp_data[~i_nan]
        temp_coords = temp_coords[~i_nan]

        interp = OpenEPDataInterpolator(
            method="rbf",
            distanceThreshold=distance_thresh,
            rbfConstant=rbf_constant_value,
        )
        vertex_voltage_data = interp.interpolate(x0=temp_coords, d0=temp_data, x1=pts)

    return vertex_voltage_data
