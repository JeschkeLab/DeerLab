.. highlight:: matlab
.. _aptkernel:

*********************
:mod:`aptkernel`
*********************

Computes the dipolar interaction kernel and elements required for the approximate Pake transformation (APT).

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`K = aptkernel(t,...)`
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **t** - Time axis (N-array)
Returns
    *   **K** - APT kernel elements (struct)

Usage
=========================================

.. code-block:: matlab

    K = aptkernel(t)

Computes a structure ``K`` containing the (N/2-2)xN point kernel, the (N/2-2) point array of normalization factors, N/2-2) point frequency axis and the (N/2-2)x(N/2-2) crosstalk matrix corresponding to the N-point time axis ``t``. The output structure ``K`` contains the following fields:

*   ``Base``
*   ``NormalizationFactor``
*   ``FreqAxis``
*   ``TimeAxis``
*   ``Crosstalk``

This structure can be then be passed directly to the :ref:`apt` function for computing the APT.



Optional Arguments
=========================================
Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    P = winlowpass(args,'Property1',Value1,'Property2',Value2,...)

.. centered:: **Property Names & Descriptions**

ExitationBandwidth
    The excitation bandwidth in MHz of the experiment can be passed as an option to account for it in the kernel.

    *Default:* [*empty*]

    *Example:*

    .. code-block:: matlab

       P = selregparam(args,'ExitationBandwidth',100)

