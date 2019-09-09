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
   [L,W] = regoperator(N,order)


Parameters
    *   ``N`` -  Distance domain size (scalar)
    *   ``order`` - Derivative order (scalar)
    *   ``r`` - Distance axis (N array)
Returns
    *   ``L`` - Regularization operator matrix ((N-d)xN matrix)
    *   ``W`` - Orthonormal basis for the null space of ``L`` (Nxd matrix)

Description
=========================================
The function can be called as follows

.. code-block:: matlab

   L = regoperator(N,order)

Computes the discrete approximation ``L`` to the derivative operator of order ``d`` on a regular grid with ``N`` points.

.. code-block:: matlab

   L = regoperator(r,order)

For simplicity, the distance axis ``r`` can also be passed, and the number of points ``N`` is computed from the length of ``r``.