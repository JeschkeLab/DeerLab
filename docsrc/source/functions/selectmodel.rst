.. highlight:: matlab
.. _selectmodel:


***********************
:mod:`selectmodel`
***********************

Selection of an optimal parametric model

-----------------------------



Syntax
=========================================

.. code-block:: matlab

    opt = selectmodel(models,V,r,K,'method')
    [opt,f,param,paramcis] = selectmodel(models,V,r,K,'method')
    [opt,f,param,paramcis] = selectmodel(models,{V1,V2,___},r,{K1,K2,___},'method')
    [opt,f,param,paramcis] = selectmodel(models,V,r,K,{'method1','method2',___})
    [opt,f,param,paramcis] = selectmodel(models,V,r,K,'method',param0)
    [opt,f,param,paramcis] = selectmodel(models,V,t,'method')
    [opt,f,param,paramcis] = selectmodel(models,{V1,V2,___},{t1,t2,___},'method')
    [opt,f,param,paramcis] = selectmodel(models,V,t,'method',param0)
    [opt,f,param,paramcis] = selectmodel(models,V,t,'method',param0,lb,ub)
    [opt,f,param,paramcis] = selectmodel(___,'Property',Value)


Parameters
    *   ``models`` - Input parametric models (cell array of function handles)
    *   ``V`` - Input signal (*N*-element array)
    *   ``r`` -  Distance axis (*N*-element array), in nanometers
    *   ``K`` -  Dipolar kernel (*NxM*-element array)
    *   ``t`` -  Time axis (*N*-array), in microseconds
    *   ``method`` - Model selection type(s) (string or cell array of strings)
    *   ``param0`` -  Initial parameter values for each model (cell array of numerical vectors)
    *   ``lb`` -  Lower bounds of parameter values for each model (cell array of numerical vectors)
    *   ``ub`` -  Upper bounds of parameter values for each model (cell array of numerical vectors)

Returns
    *  ``opt`` - Index of optimal parametric model (scalar)
    *  ``f`` - Evaluated model selection functionals (cell array)
    *  ``param`` - Fitted parameters for each evaluated model (cell array)
    *  ``paramcis`` - Confidence intervals for the fitted parameters for each evaluated model (cell array)


-----------------------------



Description
=========================================

.. code-block:: matlab

        opt = selectmodel({@model1,@model2,___,@modelN},V,r,K,{'aic',___})

Fits the distance distribution models ``model1``, ..., ``modelN`` to a signal ``V``, using the dipolar kernel ``K`` and distance axis ``r``. The models must be passed as a cell array of function handles. The fits are then evaluated according to the model selection criteria specified in the last input argument:

*   ``'aic'`` - Akaike information criterion
*   ``'aicc'`` - Corrected Akaike information criterion
*   ``'bic'`` - Bayesian information criterion
*   ``'rmsd'`` - Root mean square deviation


The function returns an array containing the optimal model for each selection method.


-----------------------------

.. code-block:: matlab

    opt = selectmodel({@model1,@model2,___,@modelN},V,t,{'aic',___})

Fits the time-domain parametric models ``model1``, ..., ``modelN`` to a signal ``V`` defined over the time axis ``t``.


-----------------------------

.. code-block:: matlab


    opt = selectmodel(models,{V1,V2,___},r,{K1,K2,___})
    opt = selectmodel(models,{V1,V2,___},r,{K1,K2,___},param0)

Passing multiple signals/kernels enables distance-domain global fitting of the parametric models to single distributions. The global fit weights are automatically computed according to their contribution to ill-posedness. The multiple signals are passed as a cell array of arrays of sizes *N1*, *N2*,... and a cell array of kernel matrices with sizes *N1xM*, *N2xM*, ... must be passed as well.


-----------------------------

.. code-block:: matlab


    opt = selectmodel(models,{V1,V2,___},{t1,t2,___})
    opt = selectmodel(models,{V1,V2,___},{t1,t2,___},param0)

Similarly, time-domain global fitting can be used when passing time-domain ``models`` and the model time axes ``{t1,t2,___}`` of the corresponding signals.


-----------------------------


.. code-block:: matlab

    opt = selectmodel(models,V,r,K,{'aic',___},{par1,___,parN})
    opt = selectmodel(models,V,t,{'aic',___},{par1,___,parN})


The initial guess values for the parameters of each model can be passed as a cell array ``{par1,___,parN}`` of value vectors. If not, specified the initial values are taken from the info of each parametric model.

-----------------------------


.. code-block:: matlab

    opt = selectmodel(models,V,r,K,{'aic',___},{par1,___,parN},{lb1,___,lbN},{ub1,___,ubN})
    opt = selectmodel(models,V,t,{'aic',___},{par1,___,parN},[],{ub1,___,ubN})
    opt = selectmodel(models,V,t,{'aic',___},{par1,___,parN},{lb1,___,lbN})


Similarly, the lower and upper bounds for the parameters of each model can be passed as cell arrays ``{lb1,___,lbN}`` and ``{ub1,___,ubN}``, where ``lbi`` and ``ubi`` are the lower and upper boundaries of the parameters for the i-th model, respectively. If not, specified the boundaries are taken from the info of each parametric model.


-----------------------------


.. code-block:: matlab

    [opt,f,param,paramcis] = selectmodel(___)

Additional outputs include: the evaluated method selection functionals ``f`` for the different methods, a cell array ``params`` with the fitted parameters for each of the evaluated models, as well as their confidence intervals ``paramcis`` (see :ref:`cireference`).

-----------------------------



Additional Settings
=========================================

Additional settings can be specified via name-value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.

.. code-block:: matlab

    opt = selectmodel(___,'Property1',Value1,'Property2',Value2,___)

- ``'GlobalWeights'`` - Weights for global fitting
    Array of weighting coefficients for the individual signals in global fitting regularization. If not specified, the global fit weights are automatically computed according to their contribution to ill-posedness. Weight values do not need to be normalized. The same number of weights as number of input signals is required.

    *Default:* [*empty*]

    *Example:*

		.. code-block:: matlab

			___ = selectmodel(___,'GlobalWeights',[0.1 0.6 0.3]])

See :ref:`fitparamodel` for a detailed list of other name-value pairs accepted by the function.
