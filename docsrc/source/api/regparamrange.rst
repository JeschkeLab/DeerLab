.. highlight:: matlab
.. _regparamrange:


***********************
:mod:`regparamrange`
***********************

Regularization parameter range estimator

-----------------------------


Syntax
=========================================

.. code-block:: matlab

    alphas = regparamrange(K,L)
    alphas = regparamrange(K,L,'Property',Value)


Parameters
    *   ``K`` - Dipolar kernel (*NxM*-element array)
    *   ``L`` - Regularization operator (*(M-order))xM*-element matrix)

Returns
    *   ``alphas`` - List of regularization parameters (array)

-----------------------------


Description
=========================================

.. code-block:: matlab

    alphas = regparamrange(K,L)

Determines an array of regularization parameters ``alphas`` from the generalized singular value decomposition (GSVD) of the dipolar kernel ``K`` and regularization operator ``L``.


-----------------------------



Additional Settings
=========================================

Additional settings can be specified via name-value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.


.. code-block:: matlab

    alphas = regparamrange(___,'Property1',Value1,'Property2',Value2,___)

- ``'Resolution'`` - Resolution
    Logarithmic scale resolution of the array of alpha candidates.

    *Default:* ``0.1``

    *Example:*

		.. code-block:: matlab

			alphas = regparamrange(args,'Resolution',0.01)

- ``'NoiseLevel'`` - Estimation of the noise level
    Estimation of the noise standard deviation of the signal for scaling of the singular values.

    *Default:* ``0``

    *Example:*

		.. code-block:: matlab

			alphas = regparamrange(args,'NoiseLevel',0.05)

