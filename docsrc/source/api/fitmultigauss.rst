.. highlight:: matlab
.. _fitmultigauss:


***********************
:mod:`fitmultigauss`
***********************

Multi-Gauss fitting of a distance distribution

Syntax
=========================================

.. code-block:: matlab

    P = fitmultigauss(S,K,r,Ngauss)
    [P,param,Nopt,metrics,Peval] = fitmultigauss(S,K,r,Ngauss)
    [P,param,Nopt,metrics,Peval] = fitmultigauss(S,K,r,Ngauss,'Property',Value)


Parameters
    *   ``S`` - Input signal (N-array)
    *   ``K`` -  Dipolar kernel (NxM-array)
    *   ``r`` -  Distance Axis (N-array)
    *   ``Ngauss`` - Maximum number of Gaussians (scalar)
Returns
    *  ``P`` - Distance Distribution (M-array)
    *  ``param`` - Fitted model parameters (array)
    *  ``Nopt`` - Optimal number of Gaussian (scalar)
    *  ``metrics`` - Evaluated model selection functionals (cell array)
    *  ``Peval`` - Fitted distance distributions for each multigauss model (NgaussxN matrix)
Description
=========================================

.. code-block:: matlab

        P = fitmultigauss(S,K,r,Ngauss)

Fits the dipolar signal ``S`` to a distance distribution ``P`` using a multi-gauss parametric model according to the dipolar kernel ``K`` and distance axis ``r``. The function chooses the optimal number of Gaussian distributions up to a maximum number given by ``Ngauss`` by means of the corrected Aikaike information criterion (AICC).

.. code-block:: matlab

    [P,param,Nopt,metrics] = fitmultigauss(args)

If requested alongside the distribution ``P``, the optimal fit model parameters ``param``, the optimal number of gaussians ``Nopt`` and evaluated selection metrics ``metrics`` are returned.

Optional Arguments
=========================================
Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.

.. code-block:: matlab

    P = fitmultigauss(args,'Property1',Value1,'Property2',Value2,...)

See :ref:`fitparamodel` for a detailed list of the property-value pairs accepted by the function.