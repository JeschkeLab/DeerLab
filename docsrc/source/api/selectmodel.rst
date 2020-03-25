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

    opt = selectmodel(Models,S,r,K,'method')
    [opt,f] = selectmodel(Models,S,r,K,'method')
    [opt,f,param,paramcis] = selectmodel(Models,S,r,K,'method')
    [opt,f,param,paramcis] = selectmodel(Models,S,r,K,'method',param0)
    [opt,f,param,paramcis] = selectmodel(Models,S,t,'method')
    [opt,f,param,paramcis] = selectmodel(Models,S,t,'method',param0)
    [opt,f,param,paramcis] = selectmodel(Models,S,t,'method','Property',Value)
    [opt,f,param,paramcis] = selectmodel(Models,S,r,K,'method',param0,'Property',Value)
    [opt,f,param,paramcis] = selectmodel(Models,S,t,'method',param0,'Property',Value)


Parameters
    *   ``Models`` - Input parametric models (cell array of function handles)
    *   ``S`` - Input signal (*N*-element array)
    *   ``r`` -  Distance Axis (*N*-element array)
    *   ``K`` -  Dipolar kernel (*NxM*-element array)
    *   ``t`` -  Time Axis (*N*-array)
    *   ``method`` - Model selection type(s) (string or cell array of strings)
    *   ``param0`` -  Initial parameter values for each model (cell array of numerical vectors)
Returns
    *  ``opt`` - Optimal parametric model index (scalar)
    *  ``f`` - Evaluated model selection functional (cell array)
    *  ``param`` - Fitted parameters for each evaluated model (cell array)
    *  ``paramcis`` - Fit confidence intervals for each evaluated model (cell array)


-----------------------------



Description
=========================================

.. code-block:: matlab

        opt = selectmodel({@model1,@model2,..,@modelN},S,r,K,{'aic',..})

Evaluates the fits of the parametric models ``model1``,..., ``modelN`` to a signal ``S`` according to the dipolar kernel ``K`` and distance axis ``r``. The models must be passed as a cell array of function handles. Each fit is then evaluated according to the model selection criterions specified in the last input argument.

*   ``'aic'`` - Akaike information criterion
*   ``'aicc'`` - Corrected Akaike information criterion
*   ``'bic'`` - Bayesian information criterion
*   ``'rmsd'`` - Root mean square deviation


The function returns an array containing the optimal model for each selection method.


-----------------------------


.. code-block:: matlab

        opt = selectmodel({@model1,@model2,..,@modelN},S,t,{'aic',..})

Evaluates the fits of the  time-domain parametric models ``model1``,..., ``modelN`` to a signal ``S`` according to the time axis ``t``.


-----------------------------


.. code-block:: matlab

        opt = selectmodel({@model1,@model2,..,@modelN},S,r,K,{'aic',..},{par1,..parN})
        opt = selectmodel({@model1,@model2,..,@modelN},S,t,{'aic',..},{par1,..parN})


The initial guess values for the parameters of each model can be passed as a cell array ``{par1,...parN}`` of value vectors.


-----------------------------


.. code-block:: matlab

    [opt,f,param,paramcis] = selectmodel(args)

Additional outputs include, the evaluated method selection functionals ``f`` for the different methods and a cell array ``params`` with the fitted parameters for each of the evaluated models , as well as their confidence intervals ``paramcis``.

-----------------------------



Optional Arguments
=========================================
Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.

.. code-block:: matlab

    opt = selectmodel(args,'Property1',Value1,'Property2',Value2,..)

- ``'Upper'`` - Parameter upper bound constraints
    Cell array containing the upper bound values for the parameters of the evaluated parametric models.

    *Default:* [*empty*] - Uses the model's default upper bound values

    *Example:*

		.. code-block:: matlab

			opt = selectmodel({@dd_onegauss,@dd_onerice},S,r,K,'aicc','Upper',{[10 1],[10 2]})

- ``'Lower'`` - Parameter lower bound constraints
    Cell array containing the lower bound values for the parameters of the evaluated parametric models.

    *Default:* [*empty*] - Uses the model's default lower bound values

    *Example:*

		.. code-block:: matlab

			opt = selectmodel({@dd_onegauss,@dd_onerice},S,r,K,'aicc','Lower',{[1 0.1],[10 0.2]})

See :ref:`fitparamodel` for a detailed list of other property-value pairs accepted by the function.