.. highlight:: matlab
.. _aptkernel:

*********************
:mod:`aptkernel`
*********************

Computes the dipolar interaction kernel and elements required for the approximate Pake transformation (APT).

Syntax
=========================================

.. code-block:: matlab

    K = aptkernel(t)
    K = aptkernel(t,'Property',Value)

Parameters
    *   ``t`` - Time axis (N-array)
Returns
    *   ``K`` - APT kernel elements (struct)

Description
=========================================

.. code-block:: matlab

    K = aptkernel(t)

Computes a structure ``K`` containing the following fields:

*   ``Base``: (N/2-2)xN point kernel
*   ``NormalizationFactor``: (N/2-2) point array of normalization factors
*   ``FreqAxis``: (N/2-2) point frequency axis
*   ``TimeAxis``: N-point time axis
*   ``Crosstalk``: (N/2-2)x(N/2-2) crosstalk matrix

This structure can be then be passed directly to the :ref:`apt` function for computing the APT.



Optional Arguments
=========================================
Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.

.. code-block:: matlab

    K = aptkernel(t,'Property1',Value1,'Property2',Value2,...)

ExcitationBandwidth
    The excitation bandwidth in MHz of the experiment.

    *Default:* empty, corresponding to infinite excitation bandwidth

    *Example:*

    .. code-block:: matlab

       K = aptkernel(t,'ExcitationBandwidth',100)   % 100 MHz excitation bandwidth

