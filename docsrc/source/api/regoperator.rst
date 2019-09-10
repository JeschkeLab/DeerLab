.. highlight:: matlab
.. _regoperator:

*********************
:mod:`regoperator`
*********************

Computes the discrete approximation to the derivative regularization operators.

Syntax
=========================================

.. code-block:: matlab

   L = regoperator(N,order)
   L = regoperator(r,order)


Parameters
    *   ``N`` -  Distance domain size (scalar)
    *   ``d`` - Derivative order (0, 1, or 2)
    *   ``r`` - Distance axis (N array)
Returns
    *   ``L`` - Regularization operator matrix ((N-d)xN matrix)

Description
=========================================
The function can be called as follows

.. code-block:: matlab

   L = regoperator(N,d)

Computes the discrete approximation ``L`` to the derivative operator of order ``d`` on a regular grid with ``N`` points.

.. code-block:: matlab

   L = regoperator(r,d)

For simplicity, the distance axis ``r`` can also be passed, and the number of points ``N`` is computed from the length of ``r``.