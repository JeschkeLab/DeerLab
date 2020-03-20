.. highlight:: matlab
.. _fitbackground:


**********************
:mod:`fitbackground`
**********************

Isolate a background function fit on a dipolar signal

-----------------------------


Syntax
=========================================

.. code-block:: matlab

    B = fitbackground(V,t,@model)
    [B,lambda] = fitbackground(V,t,@model)
    [B,lambda,param] = fitbackground(V,t,@model)
    [B,lambda,param,tstart] = fitbackground(V,t,@model)
    [B,lambda,param] = fitbackground(V,t,@model,tstart)
    [B,lambda,param] = fitbackground(V,t,@model,[tstart tend])
    [B,lambda,param] = fitbackground(V,t,@model,'Property',Value)
    [B,lambda,param] = fitbackground(V,t,@model,[tstart tend],'Property',Value)

Parameters
    *   ``V`` - Data to fit (*M*-element array)
    *   ``t`` - Time axis (*N*-element array)
    *   ``tfit`` - Time axis to fit (*M*-element array)
    *   ``@model`` - Background model (function handle)
    *   ``tstart`` - Time at which fit starts (scalar)
    *   ``tend`` - Time at which fit end (scalar)

Returns
    *   ``B`` - Background function (*M*-element array)
    *   ``lambda`` - Modulation depth (scalar)
    *   ``param`` - Fitted parameter values (array)


-----------------------------


Description
=========================================

.. code-block:: matlab

   [B,lambda,param,tstart] = fitbackground(V,t,@model)

Fits the background ``B`` and the modulation depth ``lambda`` to a time-domain signal ``V`` and time-axis ``t`` based on a given time-domain parametric model ``@model``. When not specified, the optimal fitting start time ''tstart'' is computed automatically by means of the :ref:`backgroundstart` function and returned as an output. The fitted parameters of the model are returned as a last output argument.

-----------------------------


.. code-block:: matlab

    [B,lambda,param] = fitbackground(V,t,@model,tstart)

The time at which the background starts to be fitted can be passed as an additional argument ``tstart``.

-----------------------------


.. code-block:: matlab

    [B,lambda,param] = fitbackground(S,t,@model,[tstart tend])

The start and end times of the fitting can be specified by passing a two-element array ``[tstart, tend]`` as the argument. If ``tend`` is not specified, the end of the signal is selected as the default.


-----------------------------


Optional Arguments
=========================================

Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    B = fitbackground(args,'Property1',Value1,'Property2',Value2,...)

- ``'ModDepth`` - Modulation depth
    Fixes the modulation depth to a user-defined value instead of fitting it along the background.

    *Default:* [*empty*] (automatically fitted)

    *Example:*

		.. code-block:: matlab

			B = fitbackground(V,t,@bg_exp,tstart,'ModDepth',0.45)


- ``'InitialGuess`` - Initial parameter values
    User-given estimation of the fit parameters, passed as an array. If not specified, the parametric model defaults are employed.

    *Default:* [*empty*]

    *Example:*

		.. code-block:: matlab

			B = fitbackground(V,t,@bg_exp,tstart,'InitialGuess',[0.75 3])


- ``'LogFit`` - Fit in log-scale
    Specifies the whether the logarithm of the signal is to be fitted.

    *Default:* ``false``

    *Example:*

		.. code-block:: matlab

			B = fitbackground(V,t,@bg_exp,tstart,'LogFit',true)

- ``'Solver'`` - Optimization solver
    Specifies the solver used for fitting the background model.

    *Default:* ``'lsqnonlin'`` (Optimization Toolbox installed) or ``'nlsqbnd'`` (Optimization Toolbox not installed)

    *Example:*

		.. code-block:: matlab

			B = fitbackground(V,t,@bg_exp,tstart,'Solver','nlsqbnd')
