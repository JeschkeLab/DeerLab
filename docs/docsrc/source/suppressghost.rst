.. highlight:: matlab
.. _suppressghost:

*********************
:mod:`suppressghost`
*********************

Ghost distance suppression in multi-spin systems

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`cS = suppressghost(S,n)`
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **S** - Signal vector (N-array)
    *   **n** - number of radicals (scalar)
Returns
    *   **cS** Corrected signal (N-array)

Usage
=========================================

.. code-block:: matlab

   cS = suppressghost(S,n)

Suppresses multi-spin contributions to the signal ``S`` by means of the power scaling approximation. The scaling is determined by the number of radicals ``n`` in the system. The function returns the power-scaled signal ``cS`` without further normalization.

