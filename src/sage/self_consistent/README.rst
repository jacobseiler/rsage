This version of ``SAGE`` can incorporate reionization in two ways.  These are
controlled by setting the ``ReionizationOn`` parameter to ``2`` or ``3``.

# `ReionizationOn = 2` 

Under this scheme, reionization is computed `pseudo-self-consistently`.  This
is done by running ``SAGE`` once without reionization (``ReionizationOn = 0``)
and then following the evolution of reionization using ``cifog``. 

The files that relate to this scheme have the prefix ``grid``; e.g., ``grid.c``
and ``grid.h``.

# `ReioniationOn = 3`

Under this scheme, reionization is computed fully self-consistently.

The files that relate to this scheme have the prefix ``selfcon``; e.g.,
``selfcon.c`` and ``selfcon.h``.
