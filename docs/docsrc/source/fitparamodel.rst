 .. highlight:: matlab
.. _fitparamodel:

*********************
:mod:`fitparamodel`
*********************
Fits a distance distribution to one (or several) signals by fitting of a parametric model.

""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`P = fitparamodel(S,K,r,@model)`
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **S** - Input signal (N-array)
    *   **K** -  Dipolar kernel (NxM-array)
    *   **r** -  Distance Axis (N-array)
    *   **model** - Parametric model (function handle)
Returns
    *  **P** - Distance Distribution (M-array)

Usage
=========================================

.. code-block:: matlab

    P = fitparamodel(S,K,r,@model)

Fitting of the N-point signal ``S`` to a M-point distance distribution ``P`` given a M-point distance axis ``r`` and NxM point kernel ``K``. The fitted distribution corresponds to a parametric model calculated by the passed function handle ``@model``.

See the documentation for constructing the dipolar kernel (:ref:`dipolarkernel`).

.. code-block:: matlab

    P = fitparamodel({S1,S2,S3},{K1,K2,S3},r,@model)

Passing multiple signals/kernels enables global fitting of the parametric model to a single distribution. The global fit weights are automatically computed according to their contribution to ill-posedness. The multiple signals are passed as a cell array of arrays of sizes N1,N2,... and a cell array of Kernel matrices with sizes N1xM,N2xM,... must be passed as well.

.. note:: The output distance distribution is already normalized by to unity integral and by the distance axis resolution

    .. math:: \mathbf{P}' = \frac{\mathbf{P}}{\sum\mathbf{P}\Delta\mathbf{r}}

Optional Arguments
=========================================
Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    P = fitparamodel(args,'Property1',Value1,'Property2',Value2,...)

.. centered:: **Property Names & Descriptions**

CostModel
    Type of fitting cost functional to use.

    * ``'lsq'`` - Least-squares fitting
    * ``'chisquared'`` - :math:`\chi^2`-fitting (as in GLADD or DD)


    *Default:* ``lsq``

    *Example:*

    .. code-block:: matlab

       P = fitparamodel(args,'CostModel','chisquared')

Solver
    Numerical solver employed for the minimization of the regularization functional models.

        *   ``'lsqnonlin'`` - Non-linear least squares
        *   ``'fminsearch'`` - Unconstrained minmization
        *   ``fmincon`` - Constrained non-linear minimization solver

    *Default:* ``'lsqnonlin'``

    *Example:*

    .. code-block:: matlab

        P = fitparamodel(args,'Solver','fmincon')

Algorithm
    Algorithm to be used by the solvers (see ``fmincon`` or ``lsqnonlin`` MATLAB documentation)

    *Default:* see MATLAB docs

    *Example:*

    .. code-block:: matlab

        P = fitparamodel(args,'Algorithm','trust-region-reflective')

GlobalWeights
    Array of weighting coefficients for the individual signals in global fitting. If not specified, the global fit weights are automatically computed according to their contribution to ill-posedness. The weights must be normalized such that the sum over all weights equals one. The same number of weights as number of input signals is required.

    *Default:* [*empty*]

    *Example:*

    .. code-block:: matlab

        P = fitparamodel({S1,S2,S3},{K1,K2,K3},r,L,'tikhonov',a,'GlobalWeights',[0.1 0.6 0.3]])

TolFun
    Optimizer function tolerance. The solver stops once the regularization functional evaluation reaches a value lower than this tolerance. Lower values increase the precision of the result, albeit at the cost of longer computation times.

    *Default:* ``1e-9``

    *Example:*

    .. code-block:: matlab

        P = fitparamodel(args,'TolFun',1e-20)

MaxIter
    Maximum number of iterations of the solver. After the solver exceeds this number the optimization will stop. This option is only relevant for the ``'fmincon'``  and ``'lsqnonneg'`` solvers.

    *Default:* ``2e7``

    *Example:*

    .. code-block:: matlab

        P = fitparamodel(args,'MaxIter',1e10)

MaxFunEval
    Maximum number of function evaluation of the solver. After the solver exceeds this number the optimization will stop. This option is only relevant for the ``'fmincon'``  and ``'lsqnonneg'`` solvers.

    *Default:* ``2e7``

    *Example:*

    .. code-block:: matlab

        P = fitparamodel(args,'MaxFunEval',1e10)
