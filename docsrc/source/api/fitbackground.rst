.. highlight:: matlab
.. _fitbackground:


**********************
:mod:`fitbackground`
**********************

Fit a parametric background model function to a dipolar time-domain signal

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
    *   ``t`` - Time axis, in microseconds (*N*-element array)
    *   ``@model`` - Background model (function handle)
    *   ``tstart`` - Time at which fit starts, in microseconds (scalar)
    *   ``tend`` - Time at which fit ends, in microseconds (scalar)

Returns
    *   ``B`` - Fitted background function evaluated over ``t`` (*M*-element array)
    *   ``lambda`` - Fitted modulation depth (scalar)
    *   ``param`` - Fitted parameter values (array)
    *   ``tstart`` - Automatically determined starting time (if not given as input), in microseconds (scalar)


-----------------------------


Description
=========================================

.. code-block:: matlab

   [B,lambda,param,tstart] = fitbackground(V,t,@model)

Fits the time-domain parametric background model ``@model`` and the modulation depth ``lambda`` to a time-domain signal ``V`` with time-axis ``t``, resulting in a fitted background ``B`` . When not specified, the optimal fitting start time ''tstart'' is computed automatically by means of the :ref:`backgroundstart` function and returned as an output. The fitted parameters of the model are returned as a last output argument.

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

Additional arguments can be specified by optional name-value pairs. All names are case insensitive and the name-value pairs can be passed in any order after the required input arguments have been passed.

.. code-block:: matlab

    B = fitbackground(args,'Property1',Value1,'Property2',Value2,...)

- ``'ModDepth`` - Modulation depth
    Fixes the modulation depth to a user-defined value instead of fitting it along with the background parameters.

    *Default:* ``[]`` (empty) (automatically fitted)

    *Example:*

		.. code-block:: matlab

			B = fitbackground(V,t,@bg_exp,tstart,'ModDepth',0.45)


- ``'InitialGuess`` - Initial parameter values
    User-given estimation of the background parameters, passed as an array. If not specified, the parametric model defaults are employed.

    *Default:* ``[]`` (empty)

    *Example:*

		.. code-block:: matlab

			B = fitbackground(V,t,@bg_exp,tstart,'InitialGuess',[0.75 3])


- ``'LogFit`` - Fit in log-scale
    Specifies whether to use the signal (``false``) or the logarithm of the signal (``true``) during fitting.

    *Default:* ``false``

    *Example:*

		.. code-block:: matlab

			B = fitbackground(V,t,@bg_exp,tstart,'LogFit',true)

- ``'Solver'`` - Optimization solver
    Specifies the solver used for fitting the background model (``lsqnonlin``, ``fminsearchcon``, ``nlsqbnd``).

    *Default:* ``'lsqnonlin'`` (Optimization Toolbox installed) or ``'fminsearchcon'`` (Optimization Toolbox not installed)

    *Example:*

		.. code-block:: matlab

			B = fitbackground(V,t,@bg_exp,tstart,'Solver','nlsqbnd')
