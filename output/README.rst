************************
RSAGE Plotting Pipeline
************************

Description
===========

This directory contains the scripts used to generate figures and plots from the
``RSAGE`` output.  The default parameter choices here reflect those used for
the Seiler et. al (2018) (TBD) paper.

At its core, ``RSAGE`` was written to compare multiple different models very
efficiently.  To this end, the plotting pipeline was created to ensure that any
arbitrary number of models could be easily and rapidly plotted.

Using the plotting pipeline
===========================

Pre-Requisites
--------------

The pipeline only supports Python 3. Due to the volume of data that ``RSAGE``
outputs, the pipeline is fully parallelized using
`mpi4py <https://mpi4py.readthedocs.io/en/stable/>`_. We require the package to
be installed to run the pipeline.

Basics
------

Plots are created from the ``RSAGE`` output using the ``paper_plots.py`` file.
From this file you can plot both galaxy (e.g., the stellar mass
function or ``fesc`` as a function of stellar mass) and reionization (e.g.,
evolution of neutral hydrogen or optical depth ``tau``) properties.

Plotting is driven entirely through the ``.ini`` files used to run ``RSAGE``
(see the `ini_files <https://github.com/jacobseiler/rsage/tree/master/ini_files>`_ 
directory for examples).  To this end, we require the directory paths within
the ``.ini`` files to be **absolute** and **not relative**.  Doing so allows us
to easily find where the output data files were written.

Setting up your plotting
------------------------

Two basic inputs you need to specify are ``output_directory`` and
``output_format``.  These control the directory where the plots are written to
and the format in which they are written.  If the specified
``output_directory`` does not exist, it will be created.

Following these are the paths to the ``.ini`` files for both the galaxy
(usually tagged with a ``_SAGE`` suffix) and reionization (usually tagged with
a ``cifog`` suffix) portions of ``RSAGE``. These ``.ini`` files are then placed
into two lists: ``gal_ini_files`` and ``reion_ini_files``. **The paths in these
lists, and only these lists, will be plotted**.

Next are ``model_tags``.  These strings will be placed on the legend of the
plots.

Plot toggles
------------

After specifying the models of interest, you need to select which plots you
want to make.  There are a series of toggles that control which plots to make
for both galaxy and reionization properties. These are followed by a series of
options for which snapshots/fixed neutral fractions etc are used to plot at.

**Galaxy Properties**
- ``gal_nion``, the evolution of the ionizing emissivity as a function of time
  since Big Bang. If you're plotting galaxy properties and want the ionizing
  emissivity use this. If you're plotting reionization properties use
  ``reion_nion``. If you're plotting both, use either!
- ``mstar_fesc``, the escape fraction as a function of galaxy stellar mass.
  Plotted in a multi-panel plot, one for each model at the snapshots specified by
  ``plot_snaps_for_models``.
- ``mstar_fej``, the fraction of baryons in the ejected reservoir as a function
  of galaxy stellar mass. Plotted in a multi-panel plot, one for each model at
  the snapshots specified by ``plot_snaps_for_models``.
- ``mstar_SFR``, the star formation rate as a function of galaxy stellar mass.
  Plotted in a multi-panel plot, one for each model at the snapshots specified
  by ``plot_snaps_for_models``.
- ``SMF``, the stellar mass function.

**Reionization Properties**
- ``history``, the evolution of the neutral hydrogen fraction as a function of
  time since the Big Bang.
- ``reion_nion``, the evolution of the ionizing emissivity as a function of time
  since Big Bang. If you're plotting reionization properties and want the
  ionizing emissivity use this. If you're plotting galaxy properties use
  ``gal_nion``. If you're plotting both, use either!
- ``ps_fixed_HI``, the 21cm power spectrum at specified fixed neutral hydrogen
  fractions (given by ``fixed_HI_values``).
- ``optical_depth``, the evolution of the Thomson optical depth ``tau`` as a
  function of time since the Big Bang.
- ``ps_scale``, the large 21cm power spectrum power as a function of ssmall
  scale power. These scales are defined by ``large_scale_def`` and
  ``small_scale_def``.
- ``slices_fixed_XHI``, slices of the ionization field at specified neutral
  hydrogen fractions (given by ``fixed_HI_values``).

Running
-------

Running the pipeline is as simple as executing the ``paper_plots.py`` driver
module.

.. code::

    $ python paper_plots.py

or in parallel on (e.g.,) 4 processors,

.. code::

    $ mpirun -np 4 python paper_plots.py
