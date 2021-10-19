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

from scipy.interpolate import Rbf
from scipy.interpolate import LinearNDInterpolator as linterp
from scipy.interpolate import NearestNDInterpolator as nearest

__all__ = [
    'get_reference_annotation',
    'get_mapping_points_within_woi',
    'get_window_of_interest',
    'get_egms_at_points',
    'plot_egm',
    'dist_between_points',
    'calculate_point_distance_max',
    'get_electrogram_coordinates',
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


def get_reference_annotation(mesh_case, *args):
    """
    Returns the value of the reference annotation

    Args:
        mesh_case (obj): openep case object
        'iegm'(str): iegm in string format
        ':' (str) or (list)  or range(start int, stop int): The electrogram point(s) for
            which the reference of annotation is required
        *args: Variable length argument list.

    Returns:
        int: ref, lists of the reference annotation for the corresponding iegm values
    """

    n_standard_args = 1
    i_egm = np.array(":")

    n_arg_in = len(args) + 1
    if n_arg_in > n_standard_args:
        for i in range(0, (n_arg_in - n_standard_args), 2):
            if np.char.lower(args[i]) == "iegm":
                i_egm = args[i + 1]
    ref = []
    if (type(i_egm) == str) and (np.char.compare_chararrays(i_egm, ":", "==", True)):
        ref = mesh_case.electric["annotations/referenceAnnot"].T

    else:
        ref_raw = mesh_case.electric["annotations/referenceAnnot"].T
        ref_raw = ref_raw.astype(int)

        for i in i_egm:
            ref.append(ref_raw[i])

    ref = np.array(ref).astype(int)
    return ref


def get_mapping_points_within_woi(mesh_case):
    """
    Returns the indices of the mapping points with
    annotated local activation time within the window of interest

    Args:
        mesh_case (obj): openep case object

    Returns:
        bool: m x 1 logical array of valid points; True if point is within woi, False otherwise.
        indexes into mesh_case.electric
    """

    # reference annotation
    reference_annot = mesh_case.electric["annotations/referenceAnnot"].T

    # woi
    woi = mesh_case.electric["annotations/woi"].T

    # Extract the region of the electrogram of interest for each mapping point, n
    woi = woi + reference_annot

    # MapAnnot
    map_annot = mesh_case.electric["annotations/mapAnnot"].T

    i_point = np.ones(shape=map_annot.shape, dtype=np.bool)

    # Remove the data points that are outside the region of interest
    for indx in range(map_annot.size):
        if (map_annot[indx] < woi[indx, 0]) or (map_annot[indx] > woi[indx, 1]):
            i_point[indx] = False
        else:
            i_point[indx] = True

    return i_point


def get_window_of_interest(mesh_case, *args):

    """
    Returns the window of interest
    Args:
        mesh_case (obj): openep case object
        'iegm'(str): iegm in string format
        ':' (str) or (list)  or range(start int, stop int): The electrogram point(s) for
            which the window of interst is required
        *args: Variable length argument list.

    Returns:
        int: mx2 arrays of window of interest
    """

    n_standard_args = 1
    i_egm = np.array(":")

    n_arg_in = len(args) + 1
    if n_arg_in > n_standard_args:
        for i in range(0, (n_arg_in - n_standard_args), 2):
            if np.char.lower(args[i]) == "iegm":
                i_egm = args[i + 1]
    woi = []
    if (type(i_egm) == str) and (np.char.compare_chararrays(i_egm, ":", "==", True)):
        woi = mesh_case.electric["annotations/woi"].T
        woi = woi.astype(int)

    else:
        woi_raw = mesh_case.electric["annotations/woi"].T
        woi_raw = woi_raw.astype(int)
        for i in i_egm:
            woi.append(woi_raw[i])

    woi = np.array(woi)

    return woi


def get_egms_at_points(mesh_case, filename, *args):
    """
    Access electrogram from stored in the openep data format

    Args:
        mesh_case (obj): openep case object
        filename(str): file path to the openep dataset .mat file
        'iegm'(str): iegm in string format
        ':' (str) or (list)  or range(start int, stop int): The electrogram point(s) for
            which the window of interst is required
        *args: Variable length argument list.

    Returns:
        float: egm_traces,corresponding list of arrays, containing the 'bip' and 'uni'
            voltages for the requested points
        string: egm_names, corresponding list of the egm names.
        int: sample_range,list of all the sample range within the window-of-interest for
            the requested points
    """

    egm_names = []
    egm_traces = []
    sample_range = []
    electric_egm_bip = []
    electric_egm_uni = []

    n_standard_args = 2
    i_egm = np.array(":")
    n_arg_in = len(args) + 1

    buffer = 50
    buffer = [-buffer, buffer]

    names = mesh_case.electric["names"]
    woi = mesh_case.electric["annotations/woi"].T.astype(np.int64)
    ref_annot = mesh_case.electric["annotations/referenceAnnot"].T.astype(np.int64)

    if n_arg_in > n_standard_args:
        for i in range(0, (n_arg_in - n_standard_args), 2):
            if np.char.lower(args[i]) == "iegm":
                i_egm = args[i + 1]

    if (type(i_egm) == str) and (np.char.compare_chararrays(i_egm, ":", "==", True)):
        electric_egm_bip = mesh_case.electric["egm"].T
        electric_egm_uni = mesh_case.electric["egmUni"].T
        sample_range = (woi[:] + ref_annot[:]) + buffer

        electric_egm_bip = np.array(electric_egm_bip).astype(float)
        electric_egm_uni = np.array(electric_egm_uni).astype(float)

        egm_traces.append(electric_egm_bip)
        egm_traces.append(electric_egm_uni[:, :, 0])
        egm_traces.append(electric_egm_uni[:, :, 1])
        egm_names.append(names)

    else:
        electric_egm_bip_raw = mesh_case.electric["egm"].T
        electric_egm_bip_raw = electric_egm_bip_raw.astype(float)
        electric_egm_uni_raw = mesh_case.electric["egmUni"].T
        electric_egm_uni_raw = electric_egm_uni_raw.astype(float)

        if len(i_egm) == 2:
            start_indx = i_egm[0]
            end_indx = i_egm[1] + 1

            if end_indx > len(electric_egm_bip_raw):
                raise IndexError("egm index out of range")

            else:
                for i in range(start_indx, end_indx):
                    egm_traces.append(electric_egm_bip_raw[i])
                    egm_traces.append(electric_egm_uni_raw[i, :, 0])
                    egm_traces.append(electric_egm_uni_raw[i, :, 1])
                    sample_range.append((woi[i] + ref_annot[i]) + buffer)
                    egm_names.append(names[i])
                egm_traces = egm_traces[::-1]

        else:
            if i_egm[0] > (len(electric_egm_bip_raw) - 1):
                raise IndexError("egm index out of range")

            else:
                egm_traces.append(electric_egm_bip_raw[i_egm])
                egm_traces.append(electric_egm_uni_raw[i_egm, :, 0])
                egm_traces.append(electric_egm_uni_raw[i_egm, :, 1])
                electric_egm_bip.append(electric_egm_bip_raw[i_egm].astype(float))
                electric_egm_uni.append(electric_egm_uni_raw[i_egm].astype(float))
                egm_names.append(names[i_egm])
                sample_range = woi[i_egm] + ref_annot[i_egm]
                sample_range = sample_range + buffer

                egm_traces = egm_traces[::-1]

    return {
        "egm_traces": egm_traces,
        "egm_names": egm_names[0],
        "sample_range": sample_range[0],
    }


def plot_egm(egmtraces, sample_range):
    """
    Plots the electrogram for a single index point

    Args:
        egmtraces(float): list of array of voltages (bip and uni)
        sample_range(int):list of the sample range within the window-of-interest for the requested single point

    Returns:
        matplotlib.pyplot: plt, plot of the electrogram voltages
    """
    # Plotting the egms
    fig, axs = plt.subplots(nrows=1, ncols=1)
    seperation = 7

    for i in range(len(egmtraces)):
        y = egmtraces[i][0][sample_range[0] : sample_range[1]]  # noqa E203
        t = np.arange(sample_range[0], sample_range[1], 1)
        axs.plot(t, y + (seperation * i))
        axs.get_yaxis().set_visible(False)
    plt.show()


def dist_between_points(A, B):
    """
    Returns the distance from A to B. A and B are specified
    as row vectors [x, y, z] or matrices, with rows representing different
    points. If npoints in A and B are different A must specify one and only
    one point.

    Args:
        A(float): row vectors of size mx3
        B(float): row vectors of size mx3

    Returns:
        float: d,returns an array of distance between the points in the row vectors

    """

    diff_sq = np.square(np.subtract(A, B))
    d = np.sqrt(np.sum(diff_sq, axis=1)).reshape(B.shape[0], 1)

    return d


def calculate_point_distance_max(points, test_points, max_distance):
    results = []
    dists = []

    for p in test_points:
        dist = np.linalg.norm(points - p, axis=1)
        inds = np.argwhere(dist <= max_distance).flatten()
        results.append(inds)
        dists.append(dist)

    return results, dists


def get_electrogram_coordinates(mesh_case, *args):
    """
    Function which returns the 3-D Cartesian Coordinates of the electrodes used to record the electrograms
    It returns the "bip" coordinates by defualt if "type" not specified when calling the function

    Args:
        mesh_case (obj): Case object from a given openep file
        "type" (str): "type" in string format
        "bip" or "uni" (str): specify the voltage type in string format
        *args: Variable length argument list.

    Returns:
        float: x, mx3 array,if "type","bip" else mx3x2 array: if "type","uni"

    Raises:
        ValueError: If "uni" voltage field is not in the openep file (mesh_case.electric.egmUniX)

    """
    n_standard_args = 1
    v_type = "bip"

    n_arg_in = len(args) + 1
    if n_arg_in > n_standard_args:
        for i in range(0, (n_arg_in - n_standard_args), 2):
            if np.char.lower(args[i]) == "type":
                v_type = args[i + 1]

    if np.char.lower(v_type) == "bip":
        x = mesh_case.electric["egmX"].T

    elif np.char.lower(v_type) == "uni":
        if not ("egmUniX" in mesh_case.electric):
            raise ValueError(
                "There is no unipolar data associated with this openep case"
            )

        else:
            x = mesh_case.electric["egmUniX"].T

    return x


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
    [idx, dists] = calculate_point_distance_max(x0, x1, smoothingLength)

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

        d = dist_between_points(cPts, self.x1)

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
    locations = get_electrogram_coordinates(mesh_case, "type", "bip")

    i_egm = mesh_case.electric["egm"].T
    i_vp = get_mapping_points_within_woi(mesh_case)
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
