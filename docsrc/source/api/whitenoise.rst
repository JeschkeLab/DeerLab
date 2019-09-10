.. highlight:: matlab
.. _whitegaussnoise:

*************************
:mod:`whitegaussnoise`
*************************

Generates Gaussian-distributed noise with uniform power spectrum distribution

Syntax
=========================================

.. code-block:: matlab

    x = whitegaussnoise(N,level)

Parameters
    *   ``N`` - Output length (scalar)
    *   ``level`` - Output noise standard deviation (scalar)

Returns
    *   ``x`` - Noise vector array (N-array)

Description
=========================================

.. code-block:: matlab

    x = whitegaussnoise(N,level)

Generates a N-point vector ``x`` of Gaussian distributed random noise. The standard deviation of the output noise is determined by the ``level`` input argument. The pseudorandom numbers are generated with a fixed random number generator seed to ensure reproducible outputs.
