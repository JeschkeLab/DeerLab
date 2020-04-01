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
    tc = correctzerotime(V,t,t0)
    [tc,t0,pos] = correctzerotime(V,t,t0)


Parameters
    *   ``V`` - Signal (*N*-element array)
    *   ``t`` - Time axis (*N*-element array)
    *   ``t0`` - Zero time (scalar)
Returns
    *   ``tc`` - Corrected time axis (*N*-element array)
    *   ``t0`` - Zero time (scalar)
    *   ``pos``  - Zero time index (scalar)


-----------------------------


Description
=========================================

.. code-block:: matlab

        [tc,t0,pos] = correctzerotime(V,t)

Determines the zero time of a dipolar signal ``V`` and corrects the time axis ``t`` in ns/us for it. If input signal ``V`` is in ns/us the zero time ``t0`` and corrected time axis ``tc`` will be returned in ns/us.


-----------------------------


.. code-block:: matlab

    tc = correctzerotime(V,t,t0)

Corrects the time axis ``t`` for a given zero-time ``t0`` in ns/us.



