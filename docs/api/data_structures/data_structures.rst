`Case` - the fundamental data structure of `openep-py`
======================================================

A `Case` stores all the information obatained during a clinical mapping
procedure.

.. autoclass:: openep.data_structures.case.Case
    :members: create_mesh, get_surface_data, get_field

Note
----

`Case` contains information on the 3D surface and associated scalar fields,
the electrograms, and the ablation sites. These data are, respectively,
stored as the following data structures: :class:`openep.data_structures.surface.Fields`,
:class:`openep.data_structures.electric.Electric`, and
:class:`openep.data_structures.ablation.Ablation`.


.. autoclass:: openep.data_structures.surface.Fields

.. autoclass:: openep.data_structures.electric.Electric

.. autoclass:: openep.data_structures.electric.Electrogram

.. autoclass:: openep.data_structures.electric.Impedance

.. autoclass:: openep.data_structures.electric.ElectricSurface

.. autoclass:: openep.data_structures.electric.Annotations

.. autoclass:: openep.data_structures.ablation.Ablation

.. autoclass:: openep.data_structures.ablation.AblationForce
