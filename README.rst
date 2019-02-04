|TRAVIS|

************************
RSAGE
************************

Description
====================

This repository contains the Reionization using Semi-Analytic Galaxy Evolution (``RSAGE``) model.  This model is an augmented version of the Semi-Analytic Galaxy Evolution (``SAGE``) model (found `here <https://github.com/darrencroton/sage>`_) that self-consistently couples galaxy evolution with the evolution of ionized gas during the Epoch of Reionization.  This coupling is powered by the semi-numerical code `cifog <https://github.com/annehutter/grid-model>`_ which follows the evolution of ionized hydrogen and helium during reionization. 

Installation
====================

Pre-Requisites
--------------------

``SAGE`` requires only `GSL <https://www.gnu.org/software/gsl/>`_ and should compile mostly out of the box.
``cifog`` has slightly more complex build requirements and we refer to the ``cifog`` `README <https://github.com/annehutter/grid-model#pre-requisities>`_ 
for full details on the required library and dependencies.

Downloading 
--------------------

.. code::

    $ git clone https://github.com/jacobseiler/rsage 
    $ git submodule init
    $ git submodule update 

``cifog`` is kept as an independant repository and managed through the use of Git Submodules. To clone the version of ``cifog`` that is used within ``RSAGE`` (kept in its own separate `fork <https://github.com/jacobseiler/grid-model>`_), the ``git submodule init`` and ``git submodule update`` commands are used following the initial clone. 

Building
--------------------

.. code::

    $ make 

Will create a single ``rsage`` executable in the ``bin/`` directory.

``RSAGE`` is fully MPI compatible which can be enabled by setting ``USE-MPI = true``
in the ``Makefile``. We also offer the ability to build ``RSAGE`` without the
self-consistent reionization by setting ``BUILD_RSAGE = false`` in the
``Makefile``.  This allows one to use the updates made to base ``SAGE`` to
better model the high redshift Universe. 

Running the code 
====================

Like its parent model, ``RSAGE`` requires halo merger trees to run.  In addition to these trees, the computation of the ionization fields with ``cifog`` require dark matter overdensities grids. For testing purposes, we have included a small set of trees and dark matter density fields (link tbd. Email jseiler@swin.edu.au for quick access!). Example parameter files for running the galaxy model and the semi-numerical reionization model are included in the ``ini_files`` directory.  For a complete description of the ``cifog`` parameters see the `README <https://github.com/jacobseiler/grid-model#parameter-file>`_.

To allow the easy handling of intermediate files (e.g., the ionizing photon and
reionization grids), several of the fields within the ``SAGE`` and ``cifog``
``.ini`` files have inputs of ``None``. Leaving these as ``None`` will force
``RSAGE`` to automatically determine the correct paths based on the recipe
flags and prescriptions chosen. **These fields are listed in the ``ini_files``
directory ``README``**.

The model can be run using

.. code::

   $ bin/rsage <SAGE ini file> <cifog ini file>

Or running in parallel,

.. code::

   $ mpirun -np <NUMPROC> bin/rsage <SAGE ini file> <cifog ini file>

``RSAGE`` was fundamentally written to compare multiple models of reionization.
To facilitate this we have included a utility script to generate multiple
``.ini`` and associated ``.slurm`` files. In the ``utils`` directory, the
``setup_runs.py`` file allows you to specify ``SAGE`` and ``cifog`` parameter
fields to change for any number of runs.  

Authors and Reference Papers
============================

The underlying semi-analytic galaxy evolution model was developed by `Darren Croton <https://github.com/darrencroton/sage>`_ and is described in `Croton et al., 2016 <https://arxiv.org/abs/1601.04709>`_.
The semi-numerical model used to compute ionization fields was developed by `Anne Hutter <https://github.com/annehutter/grid-model>`_ and is described in `Hutter, 2018 <https://arxiv.org/abs/1803.00088>`_.
The self-consistent coupling of these models was developed by Jacob Seiler and is described in (Stay tuned!). 

``RSAGE`` is undergoing constant development and improvement. Absolutely any
and all changes, improvements, questions or comments are encouraged and
requested. Feel free to open an issue for any bug-bears or difficulties or get
in touch directly via jseiler@swin.edu.au

.. |TRAVIS| image:: https://travis-ci.org/jacobseiler/rsage.svg?branch=master
       :target: https://travis-ci.org/jacobseiler/rsage
