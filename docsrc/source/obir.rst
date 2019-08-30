.. highlight:: matlab
.. _obir:

*********************
:mod:`obir`
*********************

Osher's Bregman-iterated regularization method

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:mod:`P = obir(S,K,r,'type',L,alpha,...)`
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Parameters
    *   **S** - Input signal (N-array)
    *   **K** -  Dipolar kernel (NxM-array)
    *   **r** -  Distance Axis (N-array)
    *   **type** - Regularization type (string)
    *   **L** - Regularization operator ((M-order))xM-array)
    *   **alpha** - Regularization parameter (scalar)
Returns
    *  **P** - Distance Distribution (M-array)


Usage
=========================================

.. code-block:: matlab

     P = obir(S,K,r,'type',L,alpha)

Fit a distance distribution ``P`` on a distance axis ``r`` to the signal ``S`` according to the kernel ''K''. The (M-2)xM point regularization matrix ``L`` and regularization parameter ``alpha` control the regularization properties.

The type of regularization employed in obir is set by the ``'type'`` input argument. The regularization models implemented in ``obir`` are:

*    ``'tikhonov'`` -   Tikhonov regularization
*    ``'tv'``       -   Total variation regularization
*    ``'huber'``    -   Pseudo-Huber regularization


Optional Arguments
=========================================
Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    P = obir(args,'Property1',Value1,'Property2',Value2,...)

.. centered:: **Property Names & Descriptions**

NoiseLevelAim
    Level (standard deviation) of noise at which Bregman iterations must stop. If not specified the noise is automatically computed using the :ref:`noiselevel` function.

    *Default:* [*empty*]

    *Example:*

    .. code-block:: matlab

       P = obir(args,'NoiseLevelAim',0.05)


DivergenceStop
    Specify whether the Bregman iterations must be stopped if the functional value increases instead of decreasing.

    *Default:* ``false``

    *Example:*

    .. code-block:: matlab

       P = obir(args,'DivergenceStop',true)

MaxOuterIter
   Maximal number of Bregman iterations.

    *Default:* ``5000``

    *Example:*

    .. code-block:: matlab

        P = obir(args,'MaxOuterIter',1e5)

AxisHandle
    Axis handle for plotting. If specified the state of the distance distribution at each Bregman iteration is displayed on the given axis object.

    *Default:* [*empty*]

    *Example:*

    .. code-block:: matlab

        P = obir(args,'AxisHandle',gca)

Solver
    Numerical solver employed for the minimization of the regularization functional models.

        *   ``'fnnls'`` - Fast non-negative least squares solver
        *   ``fmincon`` - Constrained non-linear minimization solver

    *Default:* ``'fnnls'``

    *Example:*

    .. code-block:: matlab

        P = obir(args,'Solver','fmincon')

TolFun
    Optimizer function tolerance. The solver stops once the regularization functional evaluation reaches a value lower than this tolerance. Lower values increase the precision of the result, albeit at the cost of longer computation times.

    *Default:* ``1e-9``

    *Example:*

    .. code-block:: matlab

        P = obir(args,'TolFun',1e-20)

MaxIter
    Maximum number of iterations of the solver. After the solver exceeds this number the optimization will stop. This option is only relevant for the ``'fmincon'`` solver.

    *Default:* ``2e7``

    *Example:*

    .. code-block:: matlab

        P = obir(args,'MaxIter',1e10)

MaxFunEval
    Maximum number of function evaluation of the solver. After the solver exceeds this number the optimization will stop. This option is only relevant for the ``'fmincon'`` solver.

    *Default:* ``2e7``

    *Example:*

    .. code-block:: matlab

        P = obir(args,'MaxFunEval',1e10)