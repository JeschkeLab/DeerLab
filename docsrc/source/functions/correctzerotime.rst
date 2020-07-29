.. highlight:: matlab
.. _correctzerotime:

***********************
:mod:`correctzerotime`
***********************

Zero-time correction of dipolar spectroscopy time axis


-----------------------------


Syntax
=========================================

.. code-block:: matlab

    tc = correctzerotime(V,t)
    [tc,t0,idx] = correctzerotime(V,t)


Parameters
    *   ``V`` - Signal (*N*-element array)
    *   ``t`` - Time axis (*N*-element array)
Returns
    *   ``tc`` - Corrected time axis (*N*-element array)
    *   ``t0`` - Zero time (scalar)
    *   ``idx``  - Zero time index (scalar)


-----------------------------


Description
=========================================

.. code-block:: matlab

        [tc,t0,idx] = correctzerotime(V,t)

Determines the zero time of a dipolar signal ``V`` and corrects the time axis ``t`` (in microseconds) for it. The function returns the corrected time axis ``tc``, the corresponding zero time ``t0`` and the index ``idx`` of the zero-time in the ``t`` array.



