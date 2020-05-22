.. highlight:: matlab
.. _regoperator:

*********************
:mod:`regoperator`
*********************

Computes the discrete approximation to the derivative regularization operators.

Syntax
=========================================

.. code-block:: matlab

   L = regoperator(r,order)


Parameters
    *   ``d`` - Derivative order (integer scalar)
    *   ``r`` - Distance axis (*N*-array)
Returns
    *   ``L`` - Regularization operator matrix (*(N-d)xN*-matrix)

Description
=========================================
The function can be called as follows

.. code-block:: matlab

   L = regoperator(r,d)

Computes the discrete approximation ``L`` to the derivative operator of order ``d`` on an arbitrarily incremented distance vector ``r`` with ``N`` points. All finite difference matrix coefficients are computed via Fernberg's method.