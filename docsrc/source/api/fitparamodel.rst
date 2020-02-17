.. highlight:: matlab
.. _fitparamodel:

*********************
:mod:`fitparamodel`
*********************

Fits a time or distance-domain parametric model to one (or several) signals.

-----------------------------


Syntax
=========================================

.. code-block:: matlab

    [param,fit] = fitparamodel(V,@model,t)
    [param,fit] = fitparamodel(V,@model,t,param0)
    [param,fit] = fitparamodel(V,@model,r,K)
    [param,fit] = fitparamodel(V,@model,r,K,param0)
    [param,fit] = fitparamodel({V1,V2,...},@model,r,{K1,K2,...},param0)
    [param,fit] = fitparamodel({V1,V2,...},@model,{t1,t2,...},param0)
    [param,fit] = fitparamodel(___,'Property',Value)


Parameters
    *   ``V`` - Input signal (*N*-element array)
    *   ``model`` - Parametric model (function handle)
    *   ``t`` -  Model time axis (*M*-element array)
    *   ``r`` -  Model distance axis (*M*-element array)
    *   ``K`` -  Dipolar kernel (*NxM*-element array)
    *   ``param0`` -  Model parameter inital guess (array)
Returns
    *  ``param`` - Fitted model paramters (array)
    *  ``fit`` - Parametric model fit (*M*-element array)


-----------------------------


Description
=========================================

.. code-block:: matlab

    [param,fit] = fitparamodel(V,@model,t)
    [param,fit] = fitparamodel(V,@model,t,param0)

Fits the **time-domain** parametric model ``@model`` to the input signal ``V`` on a time axis ``t``. User-defined inital guess values can be passed as an additional argument, if not they are automatically determined from the model. If the model is a user-defined function handle, the function will require ``param0`` to be passed.

-----------------------------


.. code-block:: matlab

    [param,fit] = fitparamodel(V,@model,r,K)
    [param,fit] = fitparamodel(V,@model,r,K,param0)

Fits the **distance-domain** parametric model ``@model`` to the input signal ``V`` on a distance axis ``r``. The dipolar kernel ``K`` is required as in input for distance-domain fitting. User-defined inital guess values can be passed as an additional argument, if not they are automatically determined from the model. If the model is a user-defined function handle, the function will require ``param0`` to be passed.

-----------------------------


.. code-block:: matlab

    P = fitparamodel({V1,V2,V3},@model,r,{K1,K2,K3})
    P = fitparamodel({V1,V2,V3},@model,r,{K1,K2,K3},param0)

Passing multiple signals/kernels enables **distance-domain global fitting** of the parametric model to a single distribution. The global fit weights are automatically computed according to their contribution to ill-posedness. The multiple signals are passed as a cell array of arrays of sizes *N1*, *N2*,... and a cell array of Kernel matrices with sizes *N1xM*, *N2xM*, ... must be passed as well.

-----------------------------


.. code-block:: matlab

    P = fitparamodel({V1,V2,V3},@model,{t1,t2,t3})
    P = fitparamodel({V1,V2,V3},@model,{t1,t2,t3},param0)

Similarly, **time-domain global fitting** can be used when passing a time-domain ``@model`` and the model time axes ``{t1,t2,...}`` of the corresponding signals.

-----------------------------


User-defined parametric models must have the following function definition structure:

.. code-block:: matlab

    Vfit = model(t,param)
    Pfit = model(r,param)
	
where the ``r`` and ``t`` depend on whether the parametric model is a distance or time-domain model, respectively. Additionally the parametric model can accept a third input argument ``idx`` as follows

.. code-block:: matlab

    Vfit = model(t,param,idx)
    Pfit = model(r,param,idx)

By doing so, ``fitparamodel`` will automatically pass the index ``idx = (1,2,...,N)`` of the input signal cell array  
``{S1,S2,...,SN}`` being currently processed. This allows for implementation different routines in the parametric model for different signals during global fitting. 


-----------------------------


Optional Arguments
=========================================
Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    param = fitparamodel(args,'Property1',Value1,'Property2',Value2,...)


- ``'CostModel'`` - Optimization cost functional
    Type of fitting cost functional to use.

    * ``'lsq'`` - Least-squares fitting
    * ``'chisquared'`` - :math:`\chi^2`-fitting (as in GLADD or DD)


    *Default:* ``lsq``

    *Example:*

		.. code-block:: matlab

			param = fitparamodel(args,'CostModel','chisquared')

- ``'Upper'`` - Parameters upper bound constraints
    Array of upper bounds for the model parameters.

    *Default:* unbounded or automatically set

    *Example:*

		.. code-block:: matlab

			param = fitparamodel(args,'Upper',[1 100])

- ``'Lower'`` - Parameters lower bound constraints
    Array of lower bounds for the model parameters.

    *Default:* unbounded or automatically set

    *Example:*

		.. code-block:: matlab

			param = fitparamodel(args,'Lower',[0 3])

- ``'Solver'`` - Optimization solver
    Numerical solver employed for the minimization of the regularization functional models.

        *   ``'lsqnonlin'`` - Non-linear least squares (requires Opt. toolbox)
        *   ``'fminsearch'`` - Unconstrained minmization (free)
        *   ``'fmincon'`` - Constrained non-linear minimization solver (requires Opt. toolbox)
        *   ``'fminsearchcon'`` - Constrained non-linear minimization solver (free)
        *   ``'nlsqbnd'`` - Non-linear least squares (free)

    *Default:* ``'lsqnonlin'`` (Optimization Toolbox installed) or ``'nlsqbnd'`` (Optimization Toolbox not installed)

    *Example:*

		.. code-block:: matlab

			param = fitparamodel(args,'Solver','fmincon')

- ``'Algorithm'`` - Numerical solver algorithm
    Algorithm to be used by the solvers (see ``fmincon`` or ``lsqnonlin`` MATLAB documentation)

    *Default:* see MATLAB docs

    *Example:*

		.. code-block:: matlab

			param = fitparamodel(args,'Algorithm','trust-region-reflective')

- ``'GlobalWeights'`` - Global analysis weights
    Array of weighting coefficients for the individual signals in global fitting. If not specified, the global fit weights are automatically computed according to their contribution to ill-posedness. The same number of weights as number of input signals is required. Weight values do not need to be normalized.

    *Default:* [*empty*]

    *Example:*

		.. code-block:: matlab

			fit = fitparamodel({S1,S2,S3},{K1,K2,K3},r,L,'tikhonov',a,'GlobalWeights',[0.1 0.6 0.3]])

- ``'TolFun'`` -  Optimizer tolerance value
    Optimizer function tolerance. The solver stops once the regularization functional evaluation reaches a value lower than this tolerance. Lower values increase the precision of the result, albeit at the cost of longer computation times.

    *Default:* ``1e-9``

    *Example:*

		.. code-block:: matlab

			param = fitparamodel(args,'TolFun',1e-20)

- ``'MaxIter'`` - Maximal solver iterations
    Maximum number of iterations of the solver. After the solver exceeds this number the optimization will stop. This option is only relevant for the ``'fmincon'``  and ``'lsqnonneg'`` solvers.

    *Default:* ``2e7``

    *Example:*

		.. code-block:: matlab

			param = fitparamodel(args,'MaxIter',1e10)

- ``'MaxFunEval'`` -  Maximal solver function evalutions
    Maximum number of function evaluation of the solver. After the solver exceeds this number the optimization will stop. This option is only relevant for the ``'fmincon'``  and ``'lsqnonneg'`` solvers.

    *Default:* ``2e7``

    *Example:*

		.. code-block:: matlab

			param = fitparamodel(args,'MaxFunEval',1e10)

- ``'Verbose'`` -  Information display
    Set the level of detail display for the solvers:

        *   ``'off'`` - No information displayed
        *   ``'final'`` - Display solver exit message
        *   ``'iter-detailed'`` - display state of solver at each iteration


    *Default:* ``'off'``

    *Example:*

		.. code-block:: matlab

			param = fitparamodel(args,'Verbose','iter-detailed')