.. highlight:: matlab
.. _overtones:


***********************
:mod:`overtones`
***********************

Calculates the analytical overtone coefficients of RIDME experiments

Syntax
=========================================

.. code-block:: matlab

    c = overtones(n,Tm,T1)


Parameters
    *   ``N`` - Maximal overtone order (scalar)
    *   ``Tm`` - Phase-memory time (scalar)
    *   ``T1`` - Longitudinal relaxation time (scalar)
Returns
    *   ``c`` - Overtone coefficients (N-array)

Description
=========================================

.. code-block:: matlab

    c = overtones(n,Tm,T1)

Computes the overtone coefficients up to ``n``-th order according to analytical equations for a given mixing time ``Tm`` and longitudinal relxation time ``T1``.
