.. highlight:: matlab
.. _fitregmodel:

*********************
:mod:`fitregmodel`
*********************
Fits a parameter-free distance distribution to one (or several) signals using regularization.

-----------------------------


Syntax
=========================================

.. code-block:: matlab

    P = fitregmodel(V,K,r)
    P = fitregmodel(V,K,r,'regtype',alpha)
    P = fitregmodel(V,K,r,'regtype','method')
    P = fitregmodel(V,K,r,'regtype',alpha,'Property',Value)
    P = fitregmodel(V,K,r,'regtype','method','Property',Value)
    [P,alpha] = fitregmodel(___)

Parameters
    *   ``V`` - Input signal (*N*-element array) or input signals (cell array)
    *   ``K`` -  Dipolar kernel matrix (*NxM* array) or matrices (cell array)
    *   ``r`` -  Distance axis (*M*-element array)
    *   ``regtype`` - Regularization type (string)
    *   ``alpha`` - Regularization parameter (scalar)
    *   ``method`` - Regularization parameter selection method (string)

Returns
    *  ``P`` - Distance distribution (*M*-element array)
    *  ``alpha`` - regularization parameter (scalar)

-----------------------------


Description
=========================================

.. code-block:: matlab

    P = fitregmodel(V,K,r)

Fits a distance distribution ``P`` to the input signal ``V`` given the kernel matrix ``K``, using Tikhonov regularization using an AIC-optimized regularization parameter.

-----------------------------

.. code-block:: matlab

    P = fitregmodel(V,K,r,'regtype',alpha)

Fits a distance distribution ``P`` to the input signal ``V`` using the regularization method specified in ``'regtype'`` and the regularization parameter value given in ``alpha``. The available values for ``'regtype'`` are

    *   ``'tikhonov'`` - Tikhonov regularization
    *   ``'tv'`` - Total variation regularization
    *   ``'huber'`` - Pseudo-Huber regularization

-----------------------------


.. code-block:: matlab

    P = fitregmodel(V,K,r,'regtype','method')

Instead of passing a numerical value for the regularization parameter ``alpha``, the name of a selection method ``method`` (e.g. ``'AIC'``, ``'BIC'``, etc.) can be passed and the regularization parameter will be automatically selected by means of the :ref:`selregparam` function.

-----------------------------


.. code-block:: matlab

    P = fitregmodel({V1,V2,...},{K1,K2,...},r,...)

Passing multiple signals and kernels enables global fitting of a kernel model to a single distribution. The global fit weights are automatically computed according to their contribution to ill-posedness. The multiple signals are passed as a cell array of arrays of sizes N1,N2,... and a cell array of kernel matrices with sizes N1xM,N2xM,... must be passed as well.

-----------------------------



Additional Settings
=========================================

Additional settings can be specified via name-value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.


.. code-block:: matlab

    P = fitregmodel(___,'Property1',Value1,'Property2',Value2,___)

- ``'NonNegConstrained'`` - Non-negativity constraint
    Specifies whether the distance distribution ``P`` is to be computed under the non-negativity constraint. If the constraint is lifted, the distance distribution is computed according to the analytical solution of the inverse problem and does not require any numerical solver.

    *Default:* ``true``

    *Example:*

		.. code-block:: matlab

			P = fitregmodel(___,'NonNegConstrained',false)

- ``'HuberParam'`` - Huber parameter value
    Value of the super-parameter used in pseudo-Huber regularization.

    *Default:* ``1.35``

    *Example:*

		.. code-block:: matlab

				P = fitregmodel(___,'HuberParam',2.5)

- ``'RegOrder'`` - Regularization matrix order
    Order of the regularization operator matrix.

    *Default:* ``2``

    *Example:*

		.. code-block:: matlab

			P = fitregmodel(___,'RegOrder',3)


- ``'GlobalWeights'`` - Weights for global fitting
    Array of weighting coefficients for the individual signals in global fitting regularization. If not specified, the global fit weights are automatically computed according to their contribution to ill-posedness. Weight values do not need to be normalized. The same number of weights as number of input signals is required.

    *Default:* [*empty*]

    *Example:*

		.. code-block:: matlab

			P = fitregmodel({S1,S2,S3},{K1,K2,K3},r,L,'tikhonov',a,'GlobalWeights',[0.1 0.6 0.3]])

- ``'Solver'`` - Optimization solver
    Numerical solver employed for solving the regularized optimization problem.

        *   ``'fnnls'`` - Fast non-negative least squares solver
        *   ``'bppnnls'`` - Block principal pivoting non-negative least-squares solver
        *   ``'lsqnonneg'`` - Non-negative least-squares solver
        *   ``'fmincon'`` - Constrained non-linear minimization solver

    *Default:* ``'fnnls'``

    *Example:*

		.. code-block:: matlab

				P = fitregmodel(___,'Solver','fmincon')

- ``'TolFun'`` - Optimizer tolerance value
    Optimizer function tolerance. The solver stops once the regularization functional evaluation reaches a value lower than this tolerance. Lower values increase the precision of the result, albeit at the cost of longer computation times.

    *Default:* ``1e-9``

    *Example:*

		.. code-block:: matlab

			P = fitregmodel(___,'TolFun',1e-20)

- ``'MaxIter'`` - Maximal solver iterations
    Maximum number of iterations of the solver. After the solver exceeds this number the optimization will stop. This option is only relevant for the ``'fmincon'``  and ``'lsqnonneg'`` solvers.

    *Default:* ``2e7``

    *Example:*

		.. code-block:: matlab

			P = fitregmodel(___,'MaxIter',1e10)

- ``'MaxFunEval'`` - Maximal solver function evaluations
    Maximum number of function evaluation of the solver. After the solver exceeds this number the optimization will stop. This option is only relevant for the ``'fmincon'``  and ``'lsqnonneg'`` solvers.

    *Default:* ``2e7``

    *Example:*

		.. code-block:: matlab

			P = fitregmodel(___,'MaxFunEval',1e10)

- ``'Verbose'`` - Information display
    Set the level of detail display for the solvers:

        *   ``'off'`` - No information displayed
        *   ``'final'`` - Display solver exit message
        *   ``'iter-detailed'`` - Display state of solver at each iteration


    *Default:* ``'off'``

    *Example:*

		.. code-block:: matlab

			P = fitregmodel(___,'Verbose','iter-detailed')

