.. highlight:: matlab
.. _regparamrange:


***********************
:mod:`regparamrange`
***********************

Regularization parameter range estimator

Syntax
=========================================

.. code-block:: matlab

    alphas = regparamrange(K,L)


Parameters
    *   ``K`` - Dipolar kernel (NxM-array)
    *   ``L`` - Regularization operator ((M-order))xM-array)

Returns
    *   ``alphas`` - Regularization parameter candidates (array)

Description
=========================================

.. code-block:: matlab

    alphas = regparamrange(K,L)

Estimates an array of regularization parameter candidates ``alphas`` from the generalized singular value decomposition (GSVD) of the dipolar kernel ``K`` and regularization operator ``L``.

Optional Arguments
=========================================
Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    alphas = regparamrange(args,'Property1',Value1,'Property2',Value2,...)

NoiseDeviation
    Estimation of the noise standard deviation of the signal for scaling of the singular values.

    *Default:* ``0``

    *Example:*

    .. code-block:: matlab

       alphas = regparamrange(args,'NoiseDeviation',0.05)

logResolution
    Logarithmic scale resolution of the array of alpha candidates.

    *Default:* ``0.1``

    *Example:*

    .. code-block:: matlab

       alphas = regparamrange(args,'logResolution',0.01)