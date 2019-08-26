.. highlight:: matlab
.. _time2dist:

*********************
:mod:`time2dist`
*********************

Conversion from time-axis to distance-axis

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`r = time2dist(t)`
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **t** - Signal vector (N-array)
    *   **M** - Output length (Scalar)

Returns
    *   **r** - Distance Axis (M-array)

Usage
=========================================

.. code-block:: matlab

    [r,rmin,rmax] = time2dist(t)

Computes the N-point distance axis ``r`` according to the input time axis ``t``. The minimal and maximal distances ``rmin``, ``rmax`` are determined by the empirical approximations derived by Gunnar Jeschke as implemented in the older DeerAnalysis versions.

.. code-block:: matlab

    [r,rmin,rmax] = time2dist(t,M)

The length of the output axis can be specified by the parameter ``M``.

These empirical equation approximate the minimal and maximal detectable distances given a certain timestep :math:`\Delta t` and trace length :math:`t_\text{max}`.

.. math:: r_\text{min} = 4\left( \frac{4\Delta t \nu_0}{0.85} \right)^{1/3}

.. math:: r_\text{max} = 6\left( \frac{t_\text{max}}{2} \right)^{1/3}

where :math:`\nu_0` = 52.04 MHz is the dipolar frequency of between two nitroxide electron spins separated by 1 nm.