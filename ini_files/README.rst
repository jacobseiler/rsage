************************
Ini Files
************************

Description
====================

This directory contains the ``.ini`` files used to run ``RSAGE``. A description
of the ``SAGE`` parameters are included in the ``.ini`` file itself. For a
complete description of the ``cifog`` parameters see the
`README <https://github.com/jacobseiler/grid-model#parameter-file>`_.

None Fields
====================

Some fields can have their entries set to ``None`` to have their values
automatically detected by ``RSAGE``.  This is used to handle intermediate files
that depend upon model parameters (e.g., the ionizing photon files are named
differently depending upon the value of ``fescPrescription``). The following
``_SAGE.ini`` fields allow this:

* ``GalaxyOutputDir`` -> ``<OutputDir>/galaxies``
* ``GridOutputDir`` -> ``<OutputDir>/grids``
* ``PhotoionDir`` -> ``<OutputDir>/grids/cifog``
* ``PhotoionName`` -> ``<RunPrefix>_photHI``
* ``ReionRedshiftName`` -> ``<RunPrefix>_reionization_redshift``

The following ``_cifog.ini`` fields allow this:

* ``inputNionFile`` -> ``<OutputDir>/grids/nion/<RunPrefix>_<NionPrefix>_nionHI`` where ``<NionPrefix>`` depends upon ``fescPrescription`` and the constants used. See ``get_nion_prefix()`` in ``src/sage/self_consistent/selfcon_grid.c`` for a full description.
* ``output_XHII_file`` -> ``<OutputDir>/grids/cifog/<RunPrefix>_XHII``
* ``output_photHI_file`` -> ``<OutputDir>/grids/cifog/<RunPrefix>_photHI``
* ``output_restart_file`` -> ``<OutputDir>/grids/cifog/<RunPrefix>_restart``

.. You will need to update a number of fields to correctly point to where you have
..  saved your simulation trees and the dark matter density fields. For the
.. ``_SAGE.ini`` file, these fields are all above the **``Recipe Flags``** parameters.
.. For the ``_cifog.ini`` file, these fields are:

.. * ``redshiftFile``
.. * ``inputIgmDensityFile``
.. * ``inputIgmDenstiySuffix``
 
.. Directory fields that contain a value of ``None`` (e.g., ``inputNionFile``),
.. will be determined automatically depending on your selection of the
.. ``fescPrescription`` and corresponding constants.

.. Finally, you will need to carefully go through the ``.ini`` files and update
.. the simulation specific constants such as ``Omega``, the Box Size, etc.
