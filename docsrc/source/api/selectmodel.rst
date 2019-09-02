.. highlight:: matlab
.. _selectmodel:


***********************
:mod:`selectmodel`
***********************

Selection of an optimal parametric model

Syntax
=========================================

.. code-block:: matlab

    opt = selectmodel(Models,S,r,K,'method')
    [opt,f] = selectmodel({@model1,@model2,..,@modelN},S,r,K,{'method1',..})
    [opt,f] = selectmodel({@model1,@model2,..,@modelN},S,r,K,{'aic',..},'Property',Value)

Parameters
    *   ``Models`` - Input parametric models (cell array of function handles)
    *   ``S`` - Input signal (N-array)
    *   ``r`` -  Distance Axis (N-array)
    *   ``K`` -  Dipolar kernel (NxM-array)
    *   ``method`` - Model selection type(s) (string)
Returns
    *  ``opt`` - Optimal parametric model index (scalar)
    *  ``f`` - Evaluated model selection functional (array)


Description
=========================================

.. code-block:: matlab

        opt = selectmodel({@model1,@model2,..,@modelN},S,r,K,{'aic',..})

Evaluates the fits of the parametric models ``model1``,..., ``modelN`` to a signal ``S`` according to the dipolar kernel ``K`` and distance axis (r). The models must be passed as a cell array of function handles. Each fit is then evaluated according to the model selection criterions specified in the last input argument.

*   ``'aic'`` - Akaike information criterion
*   ``'aicc'`` - Corrected Akaike information criterion
*   ``'bic'`` - Bayesian information criterion

 Function returns an array containing the optimal model for each selection method.

.. code-block:: matlab

    [opt,f] = selectmodel(args)

The model selection functionals for the different methods can be requested as a second output.

Optional Arguments
=========================================
Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.

.. code-block:: matlab

    opt = selectmodel(args,'Property1',Value1,'Property2',Value2,...)



See :ref:`fitparamodel` for a detailed list of the property-value pairs accepted by the function.