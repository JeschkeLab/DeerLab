.. highlight:: matlab
.. _longpass:


***********************
:mod:`longpass`
***********************

Filter short distances by means of a distance-domain lowpass FIR filter (longpass)

-----------------------------


Syntax
=========================================

.. code-block:: matlab

    X = longpass(t,S)
    X = longpass(t,S,rpass)
    X = longpass(t,S,rpass,st)


Parameters
    *   ``t`` - Time axis (*N*-element array)
    *   ``S`` - Signal (*N*-element array)
    *   ``rpass`` - Passband distance (scalar)
    *   ``st`` - Transition band steepness (scalar)

Returns
    *   ``X`` - Filtered signal (*N*-element array)

-----------------------------


Description
=========================================

.. code-block:: matlab

        X = longpass(t,S)

Applies a lowpass filter to the signal ``S`` with a N-point time axis ``t`` to suppress any distances below 1.5 nm.


-----------------------------


.. code-block:: matlab

    X = longpass(t,S,rpass)

The passband distance can be specified in nm by ``rpass``. The transition band is computed according to a steepness of 0.8.


-----------------------------


.. code-block:: matlab

    X = longpass(t,S,rpass,st)

The steepness ``st`` can be passed as well as a value between 0 and 1.