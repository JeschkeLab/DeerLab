.. highlight:: matlab
.. _time2freq:

*********************
:mod:`time2freq`
*********************

Conversion from time-axis to frequency-axis

-----------------------------


Syntax
=========================================

.. code-block:: matlab

    nu = time2freq(t)
    nu = time2freq(t,M)

Parameters
    *   ``t`` - Signal vector (*N*-element array)
    *   ``M`` - Output length (scalar)

Returns
    *   ``nu`` - Frequency axis (*M*-element array)

-----------------------------


Description
=========================================

.. code-block:: matlab

    nu = time2freq(t)

Computes the N-point frequency axis ``nu`` from the N-point input time axis ``t`` according to the Nyquist criterion.

-----------------------------


.. code-block:: matlab

    nu = time2freq(t,M)

The length of the output axis can be specified by the parameter ``M``.


