.. highlight:: matlab

*********************
:mod:`noiselevel`
*********************

Estimate the noise level (standard deviation) on a signal.

.. function:: Level = noiselevel(S,M)

    :param S: Signal vector (N-array)
    :param M: Number of points (scalar)
    :returns: - Noise level (scalar)

Usage
=========================================
The function can be called as follows

.. code-block:: matlab

   Level = noiselevel(S)

Returns the standard deviation estimation of the noise in the signal ``S``. The estimation is done from the last quarter of the signal, i.e. ``M=3/4*N``.

.. code-block:: matlab

    Level = noiselevel(S,M)

If parameter ``M`` is specified, the noise level estimation is done from the last ``M`` points of the N-point signal.

.. image:: ./images/noiselevel1.svg

.. hint:: This function returns the estimation of the noise standard deviation (called noise level in DeerAnalysis2). The value returned by this function is the same type of noise level requested by many other function in DeerAnalysis2.