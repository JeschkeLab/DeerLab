.. highlight:: matlab
.. _correctzerotime:

***********************
:mod:`correctzerotime`
***********************

Zero-time correction of dipolar spectroscopy time axis

Syntax
=========================================

.. code-block:: matlab

    ct = correctzerotime(S,t)
    ct = correctzerotime(S,t,zt)
    [ct,zt,pos,norm] = correctzerotime(S,t,zt)


Parameters
    *   ``S`` - Signal (N-array)
    *   ``t`` - Time axis (N-array)
    *   ``zt`` - Zero time (scalar)
    *   ``norm`` - Normalization factor (scalar)
Returns
    *   ``ct`` - Corrected time axis (N-array)
    *   ``zt`` - Zero time (scalar)
    *   ``pos``  - Zero time index (scalar)

Description
=========================================

.. code-block:: matlab

        [ct,zt,pos,norm] = correctzerotime(S,t)

Determines the zero time of a dipolar signal ``S`` and corrects the time axis ``t`` in ns/us for it. If input signal ``S`` is in ns/us the zero time ``zt`` and corrected time axis ``ct`` will be returned in ns/us. The normalization factor ``norm`` corresponds to the corresponding signal intensity at the optimized zero-time and can be used for proper normalization of the signal.

.. code-block:: matlab

    ct = correctzerotime(S,t,zt)

Corrects the time axis ``t`` for a given zero-time ``zt`` in ns/us.

.. image:: ../images/correctzerotime1.svg
