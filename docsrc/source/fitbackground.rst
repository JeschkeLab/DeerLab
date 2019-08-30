.. highlight:: matlab
.. _fitbackground:


**********************
:mod:`fitbackground`
**********************

Fit the background function in a signal

Syntax
=========================================

---------

.. code-block:: matlab

    B = fitbackground(S,t,@model)
    [B,lambda] = fitbackground(S,t,@model)
    [B,lambda,param] = fitbackground(S,t,@model)
    [B,lambda,param] = fitbackground(S,t,@model,tstart)
    [B,lambda,param] = fitbackground(S,t,@model,[tstart tend])
    [B,lambda,param] = fitbackground(S,t,@model,[tstart tend],options)

Parameters
    *   ``S`` - Data to fit (M-array)
    *   ``t`` - Time axis (N-array)
    *   ``tfit`` - Time axis to fit (M-array)
    *   ``model`` - Background model (function handle)
    *   ``tstart`` - Time at which fit starts (scalar)
    *   ``tstart`` - Time at which fit end (scalar)

Returns
    *   ``B`` - Spectrum (M-array)
    *   ``lambda`` - Modulation depth (scalar)
    *   ``param`` - Fitted parameter values (array)


Description
=========================================

---------

.. code-block:: matlab

   [B,lambda,param] = fitbackground(S,t,@model)

Fits the background ``B`` and the modulation depth ``lambda`` to a time-domain signal ``S`` and time-axis ``t`` based on a given time-domain parametric model ``@model``. The fitted parameters of the model are returned as a last output argument.

.. code-block:: matlab

    [B,lambda,param] = fitbackground(S,t,@model,tstart)

The time at which the background starts to be fitted can be passed as an additional argument ``tstart``.

.. code-block:: matlab

    [B,lambda,param] = fitbackground(S,t,@model,[tstart tend])

The start and end times of the fitting can be specified by passing a two-element array ``[tstart, tend]`` as the argument. If tend is not specified, the end of the signal is selected as the default.


Optional Arguments
=========================================

---------

Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    B = fitbackground(args,'Property1',Value1,'Property2',Value2,...)

**Property Names & Descriptions**

LogFit
    Specifies the whether the logarithm of the signal is to be fitted.

    *Default:* ``false``

    *Example:*

    .. code-block:: matlab

        B = fitbackground(S,t,@td_exp,tstart,'LogFit',true) %Fit the logarithm of the exponential
