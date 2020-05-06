.. highlight:: matlab
.. _winlowpass:

*********************
:mod:`winlowpass`
*********************

Windowed low-pass FIR filter

-----------------------------


Syntax
=========================================

.. code-block:: matlab

    Sout = winlowpass(S,Fstop,Fpass,Fsamp)
    [Sout,H] = winlowpass(S,Fstop,Fpass,Fsamp)
    [Sout,H] = winlowpass(___,'Property',Value)


Parameters
    *   ``S`` - Signal (*N*-element array)
    *   ``Fstop`` - Stop-band frequency (scalar)
    *   ``Fpass`` - Pass-band frequency (scalar)
    *   ``Fsamp`` - Sampling frequency (scalar)
Returns
    *   ``Sout`` - Filtered signal (*N*-element array)
    *   ``H`` - Filter transfer function parameters (array)

-----------------------------


Description
=========================================

.. code-block:: matlab

    Sout = winlowpass(S,Fstop,Fpass,Fsamp)

Filters the N-point signal ``S`` using a Kaiser-windowed lowpass FIR filter with a transition band given by the pass-band ``Fpass`` and stop-band ``Fstop`` frequencies in Hz. The filter is then constructed according to the sampling rate ``Fsamp`` of the signal in Hz.

-----------------------------


.. code-block:: matlab

    [Sout,H] = winlowpass(S,Fstop,Fpass,Fsamp)

The filter transfer function parameters ``H`` requested as a second output argument.

-----------------------------


Additional Settings
=========================================

Additional settings can be specified via name-value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.


.. code-block:: matlab

    P = winlowpass(___,'Property1',Value1,'Property2',Value2,___)

- ``'MinimalAttenuation'`` - Minimal sidelobe attenuation
    Minimal attenuation level [dB] of the first sidelobe after the stop-band.

    *Default:* ``50``

    *Example:*

		.. code-block:: matlab

			P = selregparam(___,'MinimalAttenuation',60)


- ``ForwardBackward'`` - Forward-backward filtering
    Enables/disables forward-backward filtering of the signal for zero-phase filtering.

    *Default:* ``true``

    *Example:*

		.. code-block:: matlab

			P = selregparam(___,'ForwardBackward',false)
