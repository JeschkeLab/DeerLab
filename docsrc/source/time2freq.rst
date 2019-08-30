.. highlight:: matlab
.. _time2freq:

*********************
:mod:`time2freq`
*********************

Conversion from time-axis to frequency-axis

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`nu = time2freq(t,M)`
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **t** - Signal vector (N-array)
    *   **M** - Output length (Scalar)

Returns
    *   **nu** - Fequency axis (M-array)

Usage
=========================================

.. code-block:: matlab

    nu = time2freq(t)

Computes the N-point frequency axis ``nu`` from the N-point input time axis ``t`` according to the Nyquist criterion.

.. code-block:: matlab

    nu = time2freq(t,M)

The length of the output axis can be specified by the parameter ``M``.
