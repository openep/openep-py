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
Create meshes from Circle CVI workspaces --- :mod:`openep.io.readers._circle_cvi`
=================================================================================

This module contains functions for creating a PyVista PolyData mesh
from a Circle CVI workspace and a stack of dicoms.

"""

import glob
import pathlib
from xml.dom import minidom
import h5py

import pandas as pd
import numpy as np
import pyvista
import pydicom

__all__ = []


def load_dicoms(dicoms_directory):
    """Load all dicoms from a given directory.

    Args:
        dicom_directory (pathlib.Path): Path to directory containing dicom files.
    Returns:
        dicoms (dict): Dictionary of dicoms loaded from the directory.
    """

    dicoms = {}
    dicom_files = glob.glob(f'{dicoms_directory.as_posix()}/**/*.dcm', recursive=True)

    for dicom_file in dicom_files:
        dicom = pydicom.read_file(dicom_file)
        dicoms[dicom.SOPInstanceUID] = dicom

    return dicoms


def get_contour_nodes(filename):
    """Extract nodes corresponding to contours"""

    dom = minidom.parse(filename.as_posix())
    contour_nodes = [elem for elem in dom.getElementsByTagName("Hash:item")
        if elem.getAttribute("Hash:key") == "Contours" and elem.hasChildNodes()]

    return contour_nodes


def _get_child_elements(node):
    """Get elements of an XML node."""
    return [child for child in node.childNodes if child.nodeType == child.ELEMENT_NODE]


def _get_contour_info(contour_node):
    """Extract contour xy positions and upsample factors from a node."""

    contours = {}
    contour_elements = _get_child_elements(contour_node)

    for contour_element in contour_elements:
        name = contour_element.getAttribute("Hash:key")
        for child in contour_element.getElementsByTagName('Hash:item'):

            if child.getAttribute("Hash:key") == "Points":
                points = []
                for pnode in _get_child_elements(child):
                    x = int(pnode.getElementsByTagName('Point:x')[0].firstChild.data) -1  # Adjustment for 0 numbering of python
                    y = int(pnode.getElementsByTagName('Point:y')[0].firstChild.data) -1
                    points += [[y, x]]  #Transposed to be consistent with dicomm image format.

            if child.getAttribute('Hash:key') == 'SubpixelResolution':
                upsample_factor = int(child.firstChild.data)		

        points = np.asarray(points, dtype=int)
        contours[name] = points
    return contours, upsample_factor


def get_contours(contour_nodes, dicoms, extract_epi, extract_endo):
    """Extract contacts and corresponding dicoms."""

    contour_data = []
    slice_data = []

    for contour_index, contour_node in enumerate(contour_nodes):

        image_id = contour_node.parentNode.getAttribute("Hash:key")
        dicom = dicoms[image_id]

        sliceloc = round(dicom.SliceLocation)
        phase_encode = dicom.InPlanePhaseEncodingDirection
        lge_datetime = dicom.ContentDate + dicom.ContentTime[:6]

        # TODO: check if there is duplicate data (i.e. data with the same key)
        key = [sliceloc, phase_encode, str(lge_datetime)]
        contours, upsample_factor = _get_contour_info(contour_node)

        # Skip if epi or endo data is missing
        if extract_epi and extract_endo:
            is_closed_contour = "saendocardialContour" in contours and "saepicardialContour" in contours
            is_open_contour = "saendocardialOpenContour" in contours and "saepicardialOpenContour" in contours
        elif extract_epi:
            is_closed_contour = "saepicardialContour" in contours
            is_open_contour = "saepicardialOpenContour" in contours
        elif extract_endo:
            is_closed_contour = "saendocardialContour" in contours
            is_open_contour = "saendocardialOpenContour" in contours

        is_contour = is_closed_contour or is_open_contour
        if not is_contour:
            continue

        # Store slice info as a dataframe
        # Data for each column must be in a list or array, otherwise a ValueError is raised
        current_slice = {}
        current_slice['slice_location'] = [int(sliceloc)]
        current_slice['phase_encode'] = [str(phase_encode)]
        current_slice['datetime'] = [lge_datetime]
        current_slice['dicom_path'] = [pathlib.Path(dicom.filename).resolve().as_posix()]
        current_slice['dicom_id'] = [image_id]
        current_slice['slice_thickness'] = [dicom.SliceThickness.real]
        current_slice['basal_slice'] = is_closed_contour
        pixel_spacing_x, pixel_spacing_y = dicom.PixelSpacing
        current_slice['pixel_spacing_x'] = [pixel_spacing_x]
        current_slice['pixel_spacing_y'] = [pixel_spacing_y]
        current_slice['upsample_factor'] = upsample_factor

        current_slice = pd.DataFrame.from_dict(
            current_slice,
            orient='columns',
        )

        contour_data.append(contours)
        slice_data.append(current_slice)

    slice_data = pd.concat(slice_data)

    sort_indices = np.argsort(slice_data.slice_location.values)
    slice_data = slice_data.iloc[sort_indices]
    slice_data = slice_data.reindex(copy=False)
    contour_data = np.asarray(contour_data)[sort_indices]

    return contour_data, slice_data


def _align_contours(contours_xy):
    """Align a stack of contours so they share the same center of mass.

    contours_xy (list): List of numpy arrays, with each array containing the xy
        position of points in a contour.
    """

    contour_coms = np.asarray([np.mean(contour, axis=0) for contour in contours_xy])
    mesh_com = contour_coms[1]  # why use contour_coms[1] as the com, why not e.g. contour_coms[0]?
    for contour, contour_com in zip(contours_xy, contour_coms):
        contour += (mesh_com - contour_com)

    return contours_xy


def _add_z_locations(contours_xy, z_resolution):
    """Add z values to a list of contours.

    Args:
        contours_xy (list): List of numpy arrays, each containing contours to which z
            positions be added.
        z_resolution (float): Distance in z between each slice.
    """

    contour_sizes = np.asarray([len(contour) for contour in contours_xy])
    z_locations = np.asarray([z_resolution * slice_index for slice_index in range(contour_sizes.size)])
    contour_z = np.concatenate([np.full(size, fill_value=z) for size, z in zip(contour_sizes, z_locations)])
    contours_xy = np.concatenate(contours_xy, axis=0)
    contours = np.concatenate([contours_xy, contour_z[:, np.newaxis]], axis=1)

    return contours


def _add_apex(contours, n_slices=1):
    """Add an apex to a set of contours.

    Args:
        contours (np.ndarray): Contours to which an apex will be added.
        n_slices (int): Number of slices to add. Use 1 for endo
            and 2 for epi. Defaults to 1.
    """

    # determine the height(s) in z of the new slice(s)
    unique_z = np.unique(contours[:, 2])
    z_thickness = unique_z[1]  # unique_z[0] is equal to z
    max_z = unique_z[-1]
    new_z_values = max_z + np.arange(1, n_slices+1) * z_thickness * 0.5  # why * 0.5?

    slice_points = contours[contours[:, 2]==max_z]
    n_points_per_slice = slice_points.size

    # move each contour close to the com of the mesh
    fractions_to_centre = 1 - 0.5**np.arange(1, n_slices+1)
    mesh_centre = np.mean(contours[contours[:, 2]==z_thickness], axis=0)

    for fraction_to_centre, z_value in zip(fractions_to_centre, new_z_values):

        new_points = slice_points.copy()
        new_points[:, 2] = z_value
        translate_xy = mesh_centre[:2] - new_points[:, :2]
        new_points[:, :2] += fraction_to_centre * translate_xy
        contours = np.concatenate([contours, new_points], axis=0)

    return contours


def _generate_surface_mesh(contours):
    """Create a PolyData surface mesh from an array a contour points."""

    points = pyvista.PolyData(contours)
    mesh = points.delaunay_3d()
    surface_mesh = mesh.extract_surface()

    return surface_mesh


def create_mesh(dicoms, contours_xy, align_contours=True, n_apex_slices=0):
    """Create a 3D mesh from a set of contours and info about the associated dicoms.

    Args:
        dicoms (pandas.DataFrame): DataFrame containing info about the stack of dicoms.
        contours_xy (list): List of numpy arrays, each containing the xy position of point
            in a contour.
        align_contours (bool, optional): If True, the contours will be translated to share the
            same center of mass in xy.
        n_apex_slices (int, optional): Add an apex to the mesh using this number of slices.
            Useful for adding an apex to ventricles (1 slice for endo and 2 for epi).

    Returns:
        mesh (pyvista.PolyData)
    """

    z_resolution = np.diff(dicoms.slice_location.values)[0]

    if align_contours:
        contours_xy = _align_contours(contours_xy)

    contours = _add_z_locations(contours_xy=contours_xy, z_resolution=z_resolution)

    if n_apex_slices:
        contours = _add_apex(contours, n_slices=2)

    mesh = _generate_surface_mesh(contours=contours)

    return mesh
