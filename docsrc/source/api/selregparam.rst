.. highlight:: matlab
.. _selregparam:

*********************
:mod:`selregparam`
*********************
Optimal selection of a regularization parameter value according to different model selection criteria.

-----------------------------



Syntax
=========================================

.. code-block:: matlab

    alpha = selregparam(S,K,r,'type','method')
    alpha = selregparam(S,K,r,'type','all')
    alpha = selregparam(S,K,r,'type',{'method1','method2',___})
    alpha = selregparam({S1,S2,___},{K1,K2,___},r,'type','method')
    alpha = selregparam(___,'Property',Value)
    [alpha,selfcn,alphas,res,pen] = selregparam(___)

Parameters
    *   ``S`` - Input signal (*N*-element array)
    *   ``K`` -  Dipolar kernel (*NxM*-element array)
    *   ``r`` -  Distance axis (*M*-element array)
    *   ``type`` - Regularization type (string)
    *   ``method`` - Model selection type (string)
Returns
    *   ``alpha`` - Optimal regularization parameter (scalar)
    *   ``selfcn`` - Model selection functional for given ``alphas`` (scalar)
    *   ``alphas`` - Evaluated candidate regularization parameters  (array)
    *   ``res`` - Evaluated regularization residual term  (array)
    *   ``pen`` - Evaluated regularization penalty term  (array)

-----------------------------



Description
=========================================

.. code-block:: matlab

    [alpha] = selregparam(S,K,r,'type','method')

Returns the optimal regularization parameter ``alpha`` from a range of regularization parameter candidates ``alphas``. The parameter for the regularization type given by ``'type'`` is computed based on the input signal ``S``, the dipolar kernel ``K`` and the regularization operator ``L``. The method employed for the selection of the regularization parameter can be specified as the ``'method'`` input argument. The available regularization models specified by ``'type'`` are

    *   ``'tikhonov'`` - Tikhonov regularization
    *   ``'tv'`` - Total variation regularization
    *   ``'huber'`` - Pseudo-Huber regularization

-----------------------------


.. code-block:: matlab

    [alpha,selfcn,alphas] = selregparam(S,K,r,'type',{'method1','method2',___})

If multiple selection methods are passed as a cell array of strings, the function returns ``alpha`` as an array of optimal regularization parameters corresponding to the input methods. The selection models functionals ``selfcn`` are also returned as a cell array of arrays containing the evaluated functionals of the requested models. The order of the output parameters corresponds to the order of the model strings in the input.

-----------------------------

.. code-block:: matlab

    [alpha,selfcn,alphas,res,pen] = selregparam(S,K,r,'type',{'method1','method2',___})

The vector of evaluated residual ``res`` and penalty ``pen`` terms can be requested as additional outputs. They can be used, e.g. to build the L-curve.


-----------------------------

.. code-block:: matlab

    [alpha,selfcn,alphas] = selregparam(S,K,r,'type','all')

Alternatively, the argument ``'all'`` can be passed, which will compute the optimal regularization parameter based on all the selection methods implemented in the function.


-----------------------------


.. code-block:: matlab

  alpha = selregparam({S1,S2,___},{K1,K2,___},r,'type','method')

Passing multiple signals/kernels enables selection of the regularization parameter for global fitting of the regularization model to a single distribution. The global fit weights are automatically computed according to their contribution to ill-posedness. The multiple signals are passed as a cell array of arrays of sizes *N1*, *N2*,... and a cell array of Kernel matrices with sizes *N1xM*, *N2xM*,... must be passed as well.

============ =============== ======================================================
    Available Model Selection  Criteria
-----------------------------------------------------------------------------------
 String        Acronym                      Model Selection Method
============ =============== ======================================================
``'aic'``         AIC           Akaike information criterion
``'aicc'``        AICc          Corrected Akaike information criterion
``'bic'``         BIC           Bayesian information criterion
``'cv'``          CV            Cross-validation
``'gcv'``         GCV           Generalized cross-validation
``'rgcv'``        rGCV          Robust generalized cross-validation
``'srgcv'``       srGCV         Strong-robust generalized cross-validation
``'dp'``          DP            Discrepancy principle
``'ee'``          EE            Extrapolated error
``'gml'``         GML           Generalized maximum-likelihood
``'lc'``          Lc            L-curve (curvature-based)
``'lr'``          Lr            L-curve (radius-based)
``'mcl'``         MCL           Mallows' :math:`C_L`
``'ncp'``         NCP           Normalized cumulative periodogram
``'rm'``          RM            Residual method
============ =============== ======================================================


-----------------------------




Additional Settings
=========================================

Additional settings can be specified via name-value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.



.. code-block:: matlab

    alpha = selregparam(___,'Property1',Value1,'Property2',Value2,___)


- ``'Range'`` - Regularization parameter search range
    Array of regularization parameter candidates to evaluate.

    *Default:* [*empty*] - Computes an optimal range automatically with :ref:`regparamrange`

    *Example:*

		.. code-block:: matlab

			alpha = selregparam(___,'Range',logspace(-3,4,100))


- ``'Search'`` - Regularization parameter search algorithm
    Specifies the type of algorithm used for searching the optimal regularization parameter. The possible settings are: 


		*   ``'grid'`` - Systematic search over a grid of regularization parameter values, using the grid specified in ``'Range'``.
		*   ``'golden'`` - Golden-ratio search algorithm over the interval specified in ``'Range'``.


    *Default:* ``golden``

    *Example:*

		.. code-block:: matlab

			alpha = selregparam(___,'Search','grid')


- ``'NonNegConstrained'`` - Non-negativity constraint
    Specifies whether the distance distribution ``P`` is to be computed under the non-negativity constraint. If the constraint is lifted, the distance distribution is computed according to the analytical solution of the inverse problem.

    *Default:* ``true``

    *Example:*

		.. code-block:: matlab

			alpha = selregparam(___,'NonNegConstrained',false)

- ``'HuberParam'`` - Huber parameter value
    Value of the superparameter used in pseudo-Huber regularization.

    *Default:* ``1.35``

    *Example:*

		.. code-block:: matlab

			alpha = selregparam(___,'HuberParam',2.5)

- ``'GlobalWeights'`` - Weights for global analysis
    Array of weighting coefficients for the individual signals in global fitting regularization. If not specified, the global fit weights are automatically computed according to their contribution to ill-posedness. Weight values do not need to be normalized. The same number of weights as number of input signals is required.

    *Default:* [*empty*]

    *Example:*

		.. code-block:: matlab

			alpha = selregparam(alphas,{S1,S2,S3},{K1,K2,K3},r,L,'tikhonov','aic','GlobalWeights',[0.1 0.6 0.3]])

- ``'TolFun'`` - Optimizer tolerance value
    Optimizer function tolerance. The solver stops once the regularization functional evaluation reaches a value lower than this tolerance. Lower values increase the precision of the result, albeit at the cost of longer computation times.

    *Default:* ``1e-9``

    *Example:*

		.. code-block:: matlab

			alpha = selregparam(___,'TolFun','1e-20')

- ``'RegOrder'`` - Regularization matrix order
    Order of the regularization operator (0,1, 2 or 3).

    *Default:* ``2``

    *Example:*

		.. code-block:: matlab

			alpha = selregparam(___,'RegOrder',3)

- ``'NoiseLevel'`` - Estimation of the noise level
    Level (standard deviation) of the noise in the input signal(s). If not specified, it is automatically computed via :ref:`noiselevel`. If multiple signals are passed (global fitting), the same number of noise levels must be specified. Required only for the ``'dp'`` and ``'mcl'`` selection methods.

    *Default:* [*empty*]

    *Example:*

		.. code-block:: matlab

			alpha = selregparam(___,'NoiseLevel',0.05)

