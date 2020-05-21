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
    [P,Pci,alpha,stats] = fitregmodel(___)

Parameters
    *   ``V`` - Input signal (*N*-element array) or input signals (cell array)
    *   ``K`` -  Dipolar kernel matrix (*NxM* array) or matrices (cell array)
    *   ``r`` -  Distance axis (*M*-element array)
    *   ``regtype`` - Regularization type (string)
    *   ``alpha`` - Regularization parameter (scalar)
    *   ``method`` - Regularization parameter selection method (string)

Returns
    *  ``P`` - Distance distribution (*M*-element array)
    *  ``Pci`` - Estimated confidence intervals (*2xM*-element array)
    *  ``alpha`` - regularization parameter (scalar)
    *  ``stats`` - Goodness-of-fit statistics (structure)

-----------------------------


Description
=========================================

.. code-block:: matlab

    [P,Pci] = fitregmodel(V,K,r)

Fits a distance distribution ``P`` to the input signal ``V`` given the kernel matrix ``K``, using Tikhonov regularization using an AIC-optimized regularization parameter. 

The ``Pci`` output contains the upper and lower 95%-confidence bands of the fitted distance distribution. The upper confidence band can be accessed via  ``Pci(:,1)`` and the lower confidence band via ``Pci(:,2)``.

-----------------------------

.. code-block:: matlab

    [P,Pci] = fitregmodel(V,K,r,'regtype',alpha)

Fits a distance distribution ``P`` to the input signal ``V`` using the regularization method specified in ``'regtype'`` and the regularization parameter value given in ``alpha``. The available values for ``'regtype'`` are

    *   ``'tikhonov'`` - Tikhonov regularization
    *   ``'tv'`` - Total variation regularization
    *   ``'huber'`` - Pseudo-Huber regularization

-----------------------------


.. code-block:: matlab

    [P,Pci,alpha] = fitregmodel(V,K,r,'regtype','method')


.. rst-class:: coderef

Instead of passing a numerical value for the regularization parameter ``alpha``, the name of a selection method ``method`` (e.g. ``'AIC'``, ``'BIC'``, etc.) can be passed and the regularization parameter will be automatically selected by means of the :ref:`selregparam` function. The selected optimal value is returned as the last output ``alpha``. 

-----------------------------


.. code-block:: matlab

    [P,Pci] =  = fitregmodel({V1,V2,___},{K1,K2,___},r,___)

Passing multiple signals and kernels enables global fitting of a kernel model to a single distribution. The global fit weights are automatically computed according to their contribution to ill-posedness. The multiple signals are passed as a cell array of arrays of sizes `N_1`, `N_2`,... and a cell array of kernel matrices with sizes `N_1 \times M`, `N_2 \times M`,... must be passed as well.

-----------------------------

.. code-block:: matlab

    [P,Pci,alpha,stats] = fitregmodel(___)

The ``stats`` structure provides several statistical metric which allow judgment on the quality of the fitted signal on the experimental data ``V`` and allows comparison between fits. The structure contains the following fields: 

         *   ``.chi2red`` - Reduced `\chi^2` test
         *   ``.R2`` - `R^2` test
         *   ``.RMSD`` - Root-mean squared deviation (RMSD)
         *   ``.AIC`` - Akaike information criterion
         *   ``.AICc`` - Corrected Akaike information criterion
         *   ``.BIC`` - Bayesian information criterion

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

- ``'NormP'`` - Normalize distance distribution
    Specifies whether the fitted distance distribution should be normalized (``true`` or ``false``). If set to ``true``, ``Pfit`` is normalized such that ``sum(Pfit)*mean(diff(r))==1``.

    *Default:* ``true``

    *Example:*

		.. code-block:: matlab

				P = fitregmodel(___,'NormP',false)

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

- ``'ConfidenceLevel'`` -  Level for parameter confidence bands
    Confidence level(s) of the confidence intervals computed for each fitted parameter. Must be an array containing values between 0 and 1. If more than one confidence level is requested, ``Pci`` is returned as a cell array containing the confidence intervals at the different requested levels.

    *Default:* ``0.95`` (95% confidence bands)

    *Example:*

		.. code-block:: matlab

			P = fitregmodel(___,'ConfidenceLevel',[0.99 0.5])
			Pci99 = Pci{1};
			Pci50 = Pci{2};

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

- ``'normP'`` -  Renormalization of the distance distribution
    This enables/disables the re-normalization of the fitted distance distribution such that ``sum(P)*dr=1``. 

    *Default:* ``true``

    *Example:*

		.. code-block:: matlab

			P = fitregmodel(___,'normP',false)
