.. highlight:: matlab

*********************
:mod:`fitregmodel`
*********************
Fits a distance distribution to one (or several) signals by optimization of a regularization functional model.

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`P = fitregmodel(S,K,r,L,'type',alpha,...)`
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **S** - Input signal (N-array)
    *   **K** -  Dipolar kernel (NxM-array)
    *   **r** -  Distance Axis (N-array)
    *   **L** - Regularization operator ((M-order))xM-array)
    *   **type** - Regularization type (string)
    *   **alpha** - Regularization parameter (scalar)
Returns
    *  **P** - Distance Distribution (M-array)

Usage
=========================================

.. code-block:: matlab

    P = fitregmodel(S,K,r,L,'type',alpha)

Fits a regularized distance distribution ``P``  from the input signal ``S`` according to the regularization model specified by the ``'type'`` argument. The available regularization models are

    *   ``'tikhonov'`` - Tikhonov regularization
    *   ``'tv'`` - Total variation regularization
    *   ``'huber'`` - Pseudo-Huber regularization

See the functions documentation for constructing and computing the dipolar kernel (:ref:`dipolarkernel`), regularization operator (:ref:`regoperator`) and regularization parameter (:ref:`selregparam`).

.. code-block:: matlab

    P = fitregmodel({S1,S2,S3},{K1,K2,S3},r,L,'type',alpha)

Passing multiple signals/kernels enables global fitting of the regularization model to a single distribution. The global fit weights are automatically computed according to their contribution to ill-posedness. The multiple signals are passed as a cell array of arrays of sizes N1,N2,... and a cell array of Kernel matrices with sizes N1xM,N2xM,... must be passed as well.

.. note:: The output distance distribution is already normalized by to unity integral and by the distance axis resolution

    .. math:: \mathbf{P}' = \frac{\mathbf{P}}{\sum\mathbf{P}\Delta\mathbf{r}}

Optional Arguments
=========================================
Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    P = fitregmodel(args,'Property1',Value1,'Property2',Value2,...)

.. centered:: **Property Names & Descriptions**

NonNegConstrained
    Specifies whether the distance distribution ``P`` is to be computed under the non-negativity constraint. If the constraint is lifted, the distance distribution is computed according to the analytical solution of the inverse problem and does not require any numerical solver.

    *Default:* ``true``

    *Example:*

    .. code-block:: matlab

       P = fitregmodel(args,'NonNegConstrained',false)

HuberParam
    Value of the superparameter used in the pseudo-Huber regularization.

    *Default:* ``1.35``

    *Example:*

    .. code-block:: matlab

        P = fitregmodel(args,'HuberParam',2.5)

GlobalWeights
    Array of weighting coefficients for the individual signals in global fitting regularization. If not specified, the global fit weights are automatically computed according to their contribution to ill-posedness. The weights must be normalized such that the sum over all weights equals one. The same number of weights as number of input signals is required.

    *Default:* [*empty*]

    *Example:*

    .. code-block:: matlab

        P = fitregmodel({S1,S2,S3},{K1,K2,K3},r,L,'tikhonov',a,'GlobalWeights',[0.1 0.6 0.3]])

Solver
    Numerical solver employed for the minimization of the regularization functional models.

        *   ``'fnnls'`` - Fast non-negative least squares solver
        *   ``'bppnnls'`` - Block principal pivoting non-negative least-squares solver
        *   ``'lsqnonneg'`` - Non-negative least-squares solver
        *   ``fmincon`` - Constrained non-linear minimization solver

    *Default:* ``'fnnls'``

    *Example:*

    .. code-block:: matlab

        P = fitregmodel(args,'Solver','fmincon')

TolFun
    Optimizer function tolerance. The solver stops once the regularization functional evaluation reaches a value lower than this tolerance. Lower values increase the precision of the result, albeit at the cost of longer computation times.

    *Default:* ``1e-9``

    *Example:*

    .. code-block:: matlab

        P = fitregmodel(args,'TolFun',1e-20)

MaxIter
    Maximum number of iterations of the solver. After the solver exceeds this number the optimization will stop. This option is only relevant for the ``'fmincon'``  and ``'lsqnonneg'`` solvers.

    *Default:* ``2e7``

    *Example:*

    .. code-block:: matlab

        P = fitregmodel(args,'MaxIter',1e10)

MaxFunEval
    Maximum number of function evaluation of the solver. After the solver exceeds this number the optimization will stop. This option is only relevant for the ``'fmincon'``  and ``'lsqnonneg'`` solvers.

    *Default:* ``2e7``

    *Example:*

    .. code-block:: matlab

        P = fitregmodel(args,'MaxFunEval',1e10)


