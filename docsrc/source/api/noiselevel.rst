.. highlight:: matlab
.. _noiselevel:

*********************
:mod:`noiselevel`
*********************

Estimate the noise standard deviation on a signal.

-----------------------------


Syntax
=========================================

.. code-block:: matlab

   Level = noiselevel(V)
   Level = noiselevel(V,M)

Parameters
    *   ``V`` - Signal vector (*N*-element array)
    *   ``M`` -  Number of points considered (scalar)
Returns
    *  Noise level (scalar)

-----------------------------


Description
=========================================
The function can be called as follows

.. code-block:: matlab

   Level = noiselevel(V)

Returns the standard deviation estimation of the noise in the signal ``V``. The estimation is done from the last quarter of the signal, i.e. ``M=3/4*N``.

-----------------------------


.. code-block:: matlab

    Level = noiselevel(V,M)

If parameter ``M`` is specified, the noise level estimation is done from the last ``M`` points of the N-point signal.

-----------------------------

.. code-block:: matlab

   Level = noiselevel(V)

If S is a 2D-dataset of different scans, the noise standard deviation is estimated from the deviations between scans. The second dimension of ``V`` must contain the different scans. The function then returns the standard deviation of the averaged signal not of the individual scans.