.. highlight:: matlab
.. _multigauss:


***********************
:mod:`multigauss`
***********************

Multi-Gauss fitting of a distance distribution

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`P = multigauss(S,K,r,Ngauss,...)`
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **S** - Input signal (N-array)
    *   **K** -  Dipolar kernel (NxM-array)
    *   **r** -  Distance Axis (N-array)
    *   **Ngauss** - Maximum number of Gaussians (scalar)
Returns
    *  **P** - Distance Distribution (M-array)
Usage
=========================================

.. code-block:: matlab

        P = multigauss(S,K,r,Ngauss)

Fits the dipolar signal ``S`` to a distance distribution ``P`` using a multi-gauss parametric model according to the dipolar kernel ``K`` and distance axis ``r``. The function chooses the optimal number of Gaussian distributions up to a maximum number given by ``Ngauss`` by means of the corrected Aikaike information criterion (AICC).


Optional Arguments
=========================================
Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.

.. code-block:: matlab

    F = dipolarsignal(args,'Property1',Value1,'Property2',Value2,...)

.. centered:: **Property Names & Descriptions**

See :ref:`fitparamodel` for a detailed list of the property-value pairs accepted by the function.