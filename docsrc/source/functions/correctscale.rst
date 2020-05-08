.. highlight:: matlab
.. _correctscale:

***********************
:mod:`correctscale`
***********************

Amplitude scale correction of dipolar spectroscopy signals

-----------------------------


Syntax
=========================================

.. code-block:: matlab

    correctscale(V,t)
    [Vc,V0] = correctscale(V,t)
    [Vc,V0] = correctscale(V,t,tmax)
    [Vc,V0] = correctscale(V,t,tmax,model)


Parameters
    *   ``V`` -- Signal (*N*-element array)
    *   ``t`` -- Time axis (*N*-element array), in microseconds
    *   ``tmax`` -- Time cutoff for fit range (scalar), in microseconds; default 0.2
    *   ``model`` -- model to fit: ``'deer'`` simple DEER model (default), ``'gauss'`` raised Gaussian
Returns
    *   ``Vc`` - Normalized signal (*N*-element array)
    *   ``V0`` - Fitted amplitude scale (scalar)

-----------------------------


Description
=========================================

.. code-block:: matlab

        [Vc,V0] = correctscale(V,t,tmax,model)

Determines the amplitude scale ``V0`` of a dipolar signal ``V`` and time axis ``t``. The output signal ``Vc`` is normalized by the scale via ``Vc = V/V0``. The amplitude scale is determined by fitting the model specified in ``model`` over a time interval determined by ``abs(t)<=tmax``.

.. code-block:: matlab

        correctscale(V,t,tmax,model)

If no outputs are specified, the results are plotted.

