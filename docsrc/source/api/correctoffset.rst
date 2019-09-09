.. highlight:: matlab
.. _correctoffset:

***********************
:mod:`correctoffset`
***********************

Amplitude offset correction of dipolar spectroscopy signals

Syntax
=========================================

.. code-block:: matlab

    Vc = correctoffset(V,t)
    [Vc,V0] = correctoffset(V,t)


Parameters
    *   ``V`` - Signal (N-array)
    *   ``t`` - Time axis (N-array)
Returns
    *   ``Vc`` - Normalized signal (N-array)
    *   ``V0`` - Fitted amplitude offset (scalar)

Description
=========================================

.. code-block:: matlab

        [Vc,V0] = correctoffset(V,t)

Determines the amplitude offset ``V0`` of a dipolar signal ``V`` and time axis ``t`` (in ns/us). The output signal ``Vc`` is normalized by the offset via ``Vc = V/V0``.

.. image:: ../images/correctoffset1.svg
    :align: center