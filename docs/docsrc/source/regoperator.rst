.. highlight:: matlab
.. _regoperator:

*********************
:mod:`regoperator`
*********************

Computes the discrete approximation to the derivative regularization operators.

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`[L,W] = regoperator(N,order)`
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **N** -  Distance domain size (scalar)
    *   **order** - Derivative order (scalar)
Returns
    *   **L** - Regularization operator matrix ((N-d)xN matrix)
    *   **W** - orthonormal basis for the null space of L (Nxd matrix)

Usage
=========================================
The function can be called as follows

.. code-block:: matlab

   L = regoperator(N,order)

Computes the discrete approximation ``L`` to the derivative operator of order ``d`` on a regular grid with ``N`` points.
