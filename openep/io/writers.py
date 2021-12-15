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
Saving datasets --- :mod:`openep.io.writers`
============================================

This module contains functions to export OpenEP data sets.

Example of saving a dataset
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Mesh data can be exported to openCARP format as follows:

.. code:: python

    import openep
    from openep._datasets.openep_datasets import DATASET_2_V73

    case = openep.load_openep_mat(DATASET_2_V73)
    openep.export_openCARP(
        case=case,
        prefix="dataset_2",
    )

This will save the mesh data from the case dataset into openCARP files, i.e.
the points data will be stored in `dataset_2.pts` and the vertex data will be
stored in `dataset_2.elem`.

Warning
-------
Currently, the fibre orientation is not calculated and written to file. Nor is the mesh
processed to make it suitable for performing openCARP simulations.

.. autofunction:: export_openCARP

"""

import pathlib
import numpy as np

from openep.data_structures.case import Case

__all__ = ["export_openCARP"]

def export_openCARP(
    case: Case,
    prefix: str,
):
    """Export mesh data from an OpenEP data to openCARP format.

    Args:
        case (Case): dataset to be exported
        prefix (str): filename prefix for writing mesh data to files
    """

    output_path = pathlib.Path(prefix).resolve()

    # Save points info
    np.savetxt(
        output_path.with_suffix(".pts"),
        case.points,
        fmt="%.6f",
        header=str(case.points.size),
        comments='',
    )

    # Save elements info
    n_triangles = case.indices.shape[0]
    cell_type = np.full(n_triangles, fill_value="Tr")
    region = np.full(n_triangles, 100, dtype=int)

    elements = np.concatenate(
        [
            cell_type[:, np.newaxis],
            case.indices,
            region[:, np.newaxis],
        ],
        axis=1,
        dtype=object,
    )

    np.savetxt(
        output_path.with_suffix(".elem"),
        elements,
        fmt="%s %d %d %d %d",
        header=str(n_triangles),
        comments='',
    )
