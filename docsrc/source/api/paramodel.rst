.. highlight:: matlab
.. _paramodel:

***********************
:mod:`paramodel`
***********************

Converts a function handle to a DeerLab parametric model

Syntax
=========================================

.. code-block:: matlab

    mod = paramodel(fcn,n)
    mod = paramodel(fcn,param0)
    mod = paramodel(fcn,param0,upper,lower)

Parameters
    *   ``fcn`` - Model function (function handle)
    *   ``n`` - Number of model parameters (scalars)
    *   ``param0`` - Number of model parameters (M-array)
    *   ``upper`` - Upper bounds of parameters (M-array)
    *   ``lower`` - Lower bounds of parameters (M-array)

Returns
    *   ``mod`` - Parametric model (function handle)

Description
=========================================

.. code-block:: matlab

    mod = paramodel(fcn,n)

Converts the input function handle ``fcn`` to a valid DeerLab parametric model compatible with the fit functions, e.g. :ref:`fitparamodel`. The number of parameters in the resulting parametric model must be specified by ``n``. The initial guess values of all parameters are set to zero. The resulting parametric model ``mod`` is returned in as a function handle.

.. code-block:: matlab

        mod = paramodel(fcn,param0)

If an array of initial guess values ``param0`` is passed, these are set into the resulting parametric model. The number of parameters is computed from the length of the array. In this case, the model is unconstrained and no bounds are enforced upon the parameter values.

.. code-block:: matlab

        mod = paramodel(fcn,param0,upper,lower)

Two arrays ``upper`` and ``lower`` containing the bounds on the parameters can be passed as additional arguments. The parametric model will then be constrained by these boundaries. The ``upper``, ``lower`` and ``param0`` arrays must be equally long.


Example
=========================================

.. code-block:: matlab

        K = dipolarkernel(r,t)
        %Define a time-domain model, signal+background
        fcn = @(t,p) bg_exp(t,p(1)).*(K*dd_onegauss(r,p(2:3)))
        %Set initial guess values
        param0 = [0.25,0.5,0.1];
        mod = paramodel(fcn,param0);
        %Fit model
        I = eye(size(K));
        Vfit = fitparamodel(Vexp,I,t,mod);