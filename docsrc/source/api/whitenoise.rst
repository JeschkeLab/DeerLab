.. highlight:: matlab
.. _whitegaussnoise:

*************************
:mod:`whitegaussnoise`
*************************

Generates Gaussian-distributed noise with uniform power spectrum distribution

-----------------------------



Syntax
=========================================

.. code-block:: matlab

    noise = whitegaussnoise(N,level)
    noise = whitegaussnoise(t,level)

Parameters
    *   ``N`` - Output length (scalar)
    *   ``t`` - Time Axis (*N*-element array)
    *   ``level`` - Output noise standard deviation (scalar)

Returns
    *   ``noise`` - Noise vector array (*N*-element array)

-----------------------------



Description
=========================================

.. code-block:: matlab

    noise = whitegaussnoise(N,level)
    noise = whitegaussnoise(t,level)


Generates a *N*-element vector ``noise`` of Gaussian distributed random noise. The standard deviation of the output noise is determined by the ``level`` input argument.

.. Important::
   The noise vector is generated from a vector of pseudo-random numbers. Each call of ``whitegaussnoise`` will return a different output noise vector. To set the output to a fixed noise realization, the random number generator must be fixed. In MATLAB this can be accomplished by calling ``rng(k)`` where ``k`` is some integer number.

