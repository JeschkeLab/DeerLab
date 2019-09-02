.. highlight:: matlab
.. _selregparam:

*********************
:mod:`selregparam`
*********************
Optimal selection of the regularization parameter according to different model selection criteria.

Syntax
=========================================

.. code-block:: matlab

    [alpha] = selregparam(alphas,S,K,L,'type','method')
    [alpha,F,alphas] = selregparam(alphas,S,K,L,'type','method')
    alpha = selregparam(alphas,S,K,L,'type','all')
    alpha = selregparam(alphas,S,K,L,'type',{'method1','method2','methodN'})
    alpha = selregparam(alphas,{S1,S2,SM},{K1,K2,KM},r,L,'type','method')
    [alpha] = selregparam(alphas,S,K,L,'type','method','Property',Value)

Parameters
    *   ``alphas`` - Candidate regularization parameters (array)
    *   ``S`` - Input signal (N-array)
    *   ``K`` -  Dipolar kernel (NxM-array)
    *   ``L`` - Regularization operator ((M-order))xM-array)
    *   ``type`` - Regularization type (string)
    *   ``method`` - Model selection type (string)
Returns
    *   ``alpha`` - Optimal regularization parameter (scalar)
    *   ``F`` - Model selection functional for given ``alphas`` (scalar)
    *   ``alphas`` - Evaluated candidate regularization parameters  (array)

Description
=========================================

.. code-block:: matlab

    [alpha] = selregparam(alphas,S,K,L,'type','method')

Returns the optimal regularization parameter ``alpha`` from a range of regularization parameter candidates ``alphas``. The parameter for the regularization type given by ``'type'`` is computed based on the input signal ``S``, the dipolar kernel ``K`` and the regularization operator ``L``. The method employed for the selection of the regularization parameter can be specified as the ``'method'`` input argument. The available regularization models specified by ``'type'`` are

    *   ``'tikhonov'`` - Tikhonov regularization
    *   ``'tv'`` - Total variation regularization
    *   ``'huber'`` - Pseudo-Huber regularization

.. code-block:: matlab

    [alpha,F,alphas] = selregparam(alphas,S,K,L,'type',{'method1','method2','method3',...})

If multiple selection methods are passed as a cell array of strings, the function returns ``alpha`` as an array of optimal regularization parameters corresponding to the input methods. The selection models functionals ``F`` are also returned as a cell array of arrays containing the evaluated functionals of the requested models. The order of the output parameters corresponds to the order of the model strings in the input.

.. code-block:: matlab

    [alpha,F,alphas] = selregparam(alphas,S,K,L,'type','all')

Alternatively, the argument ``'all'`` can be passed, which will compute the optimal regularization parameter based on all the selection methods implemented in the function.

.. code-block:: matlab

  alpha = selregparam(alphas,{S1,S2,..,SM},{K1,K2,..,KM},r,L,'type','method')

Passing multiple signals/kernels enables selection of the regularization parameter for global fitting of the regularization model to a single distribution. The global fit weights are automatically computed according to their contribution to ill-posedness. The multiple signals are passed as a cell array of arrays of sizes N1,N2,... and a cell array of Kernel matrices with sizes N1xM,N2xM,... must be passed as well.

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


Optional Arguments
=========================================
Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    P = selregparam(args,'Property1',Value1,'Property2',Value2,...)


Refine
    Specifies whether to enforce a second search around the optimal regularization parameter value with a finer grid to approach a better value of the optimum. If the refinement step does not find any minima, refinenment will descent the functional until a minima is reached. The refined search grid is included in the output ``alphas`` argument.

    *Default:* ``false``

    *Example:*

    .. code-block:: matlab

       P = selregparam(args,'Refine',true)


NonNegConstrained
    Specifies whether the distance distribution ``P`` is to be computed under the non-negativity constraint. If the constraint is lifted, the distance distribution is computed according to the analytical solution of the inverse problem.

    *Default:* ``true``

    *Example:*

    .. code-block:: matlab

       P = selregparam(args,'NonNegConstrained',false)

HuberParam
    Value of the superparameter used in the pseudo-Huber regularization.

    *Default:* ``1.35``

    *Example:*

    .. code-block:: matlab

        P = selregparam(args,'HuberParam',2.5)

GlobalWeights
    Array of weighting coefficients for the individual signals in global fitting regularization. If not specified, the global fit weights are automatically computed according to their contribution to ill-posedness. The weights must be normalized such that the sum over all weights equals one. The same number of weights as number of input signals is required.

    *Default:* [*empty*]

    *Example:*

    .. code-block:: matlab

        P = selregparam(alphas,{S1,S2,S3},{K1,K2,K3},r,L,'tikhonov','aic','GlobalWeights',[0.1 0.6 0.3]])

TolFun
    Optimizer function tolerance. The solver stops once the regularization functional evaluation reaches a value lower than this tolerance. Lower values increase the precision of the result, albeit at the cost of longer computation times.

    *Default:* ``1e-9``

    *Example:*

    .. code-block:: matlab

        P = selregparam(args,'TolFun','1e-20')

NoiseLevel
    Level (standard deviation) of the noise in the input signal(s). If not specified, it is automatically computed via :ref:`noiselevel`. If multiple signals are passed (global fitting), the same number of noise levels must be specified. Required only for the ``'dp'`` and ``'mcl'`` selection methods.

    *Default:* [*empty*]

    *Example:*

    .. code-block:: matlab

        P = selregparam(args,'mcl','NoiseLevel',0.05)
