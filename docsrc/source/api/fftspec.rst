.. highlight:: matlab
.. _fftspec:


***********************
:mod:`fftspec`
***********************

Fast-Fourier transform spectrum

Syntax
=========================================

.. code-block:: matlab

    spec = fftspec(t,S)
    [nu,spec] = fftspec(t,S)

Parameters
    *   ``t`` - Time axis (N-array)
    *   ``S`` - Signal (N-array)
Returns
    *   ``spec`` - Spectrum (M-array)

Description
=========================================

.. code-block:: matlab

    spec = fftspec(t,S)

Computes the magnitude FFT spectrum of the signal ``S`` on the time axis ``t``.

.. code-block:: matlab

    [nu,spec] = fftspec(t,S)

If two output arguments are requested, the frequency axis ``nu`` is returned as well.

Optional Arguments
=========================================
Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    spec = fftspec(args,'Property1',Value1,'Property2',Value2,...)


Type
    Type of spectrum to be returned.

    *   ``abs`` - Magnitude spectrum
    *   ``real`` - Real (absorptive) spectrum
    *   ``imag`` - Imaginary (dispersive) spectrum

    *Default:* ``abs``

    *Example:*

    .. code-block:: matlab

       spec = fftspec(args,'Type','real')


ZeroFilling
    Number of elements in the output FFT spectrum

    *Default:* ``2*length(S)``

    *Example:*

    .. code-block:: matlab

       spec = fftspec(args,'ZeroFilling',400)
