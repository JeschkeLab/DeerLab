.. highlight:: matlab
.. _winlowpass:

*********************
:mod:`winlowpass`
*********************

Windowed low-pass FIR filter

Syntax
=========================================

.. code-block:: matlab

    Sout = winlowpass(S,Fstop,Fpass,Fsamp)
    [Sout,H] = winlowpass(S,Fstop,Fpass,Fsamp)
    [Sout,H] = winlowpass(S,Fstop,Fpass,Fsamp,'Property',Value)


Parameters
    *   ``S`` - Signal (N-array)
    *   ``Fstop`` - Stopband frequency (scalar)
    *   ``Fpass`` - Passband frequency (scalar)
    *   ``Fsamp`` - Sampling frequency (scalar)
Returns
    *   ``Sout`` - Filtered signal (N-array)
    *   ``H`` - Filter transfer function parameters (array)

Description
=========================================

.. code-block:: matlab

    Sout = winlowpass(S,Fstop,Fpass,Fsamp)

Filters the N-point signal ``S`` using a Kaiser-windowed lowpass FIR filter with a transition band given by the passband ``Fpass`` and stopband ``Fstop`` frequencies in Hz. The filter is then constructed according to the sampling rate ``Fsamp`` of the signal in Hz.

.. code-block:: matlab

    [Sout,H] = winlowpass(S,Fstop,Fpass,Fsamp)

The filter transfer function parameters ``H`` requested as a second ouput argument.

Optional Arguments
=========================================
Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    P = winlowpass(args,'Property1',Value1,'Property2',Value2,...)

MinimalAttenuation
    Minimal attenuation level [dB] of the first sidelobe after the stopband.

    *Default:* ``50``

    *Example:*

    .. code-block:: matlab

       P = selregparam(args,'MinimalAttenuation',60)


ForwardBackward
    Enables/disables forward-backward filtering of the signal for zero-phase filtering.

    *Default:* ``true``

    *Example:*

    .. code-block:: matlab

       P = selregparam(args,'ForwardBackward',false)
