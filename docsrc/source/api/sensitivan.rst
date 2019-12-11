.. highlight:: matlab
.. _validate:

***********************
:mod:`sensitivan`
***********************

Sensitivity analysis via factorial design

Syntax
=========================================

.. code-block:: matlab

    med = sensitivan(fcn,fact)
    [med,up,lo,main,inter,evals] = sensitivan(fcn,fact)
    [med,up,lo,main,inter,evals] = sensitivan(fcn,fact,'Property',Value)

Parameters
    *   ``fcn`` - Function to validate (function handle)
    *   ``fact`` - Factors (struct)

Returns
    *   ``med`` - Output medians (numerical array or cell array)
    *   ``up`` - Output upper bound (numerical array or cell array)
    *   ``lo`` - Output lower bound (numerical array or cell array)
    *   ``main`` - Factor main effects (struct or cell array)
    *   ``inter`` - Factor interactions (matrix or cell array)
    *   ``evals`` - Evaluated function output arguments (cell array)


Description
=========================================

.. code-block:: matlab

    [med,up,lo] = sensitivan(fcn,fact)

Performs a sensitivity analysis of the output variables of the function ``fcn`` with respect to the parameter ranges defined in ``fact``. The output argument of ``fcn`` is evaluated for all combinations of the factors. The function to be validated must be a function handle accepting the ``fact`` struct as its first argument and only argument. 

.. code-block:: matlab

    [med,up,lo] = sensitivan(@(p)myvalfcn(p,varargin),fact)

    function [out1,out2] = myvalfcn(fact,varargin)
       out1 = process(fact.param1,fact.param2)
       out2 = process2(fact.param1,fact.param3)
    end


The median (50th-percentile), upper boundary (75th-percentile) and lower boundary (25th-percentile) of the resulting statistics are returned as the ``median``, ``up`` and ``lo`` outputs. 


.. code-block:: matlab

    [med,up,lo,main,inter] = sensitivan(fcn,fact)
	

For each factor ``fact.param1``, ``fact.param2``,... evaluated in the sensitivity analysis its main effect and interaction with the other factors will be returned as the ``main`` and ``inter`` output variables.


.. code-block:: matlab

    [med,up,lo] = sensitivan(@(p)myvalfcn(p,varargin),fact)
    med1 = med{1};
    med2 = med{2};
    main1 = main{1};
    main2 = main{2};


If the function ``fcn`` returns multiple output arguments, e.g. ``out1`` and ``out2``, all outputs from ``sensitivan`` will be returned as a cell array. Each cell will contain the results of the statistics on each of the function outputs ``out1`` and ``out2``.


.. code-block:: matlab

    [med,up,lo,main,inter,eval] = sensitivan(fcn,fact)

Additionally, a last output argument ``evals`` can be requested, a cell array, containing the ``fcn`` outputs evaluated at each parameter combination.

Optional Arguments
=========================================

Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed..

.. code-block:: matlab

    [median,iqr] = validate(fcn,valpar,'Property1',Value1,'Property2',Value2)

RandPerm
    Specifies whether to randomly permute the validation parameters combinations.

    *Default:* ``true``

    *Example:*

    .. code-block:: matlab

        [median,iqr] = validate(fcn,valpar,'RandPerm',false)

AxisHandle
    Axis handle to plot the state of the validation results at each parameter combination.

    *Default:* [*empty*]

    *Example:*

    .. code-block:: matlab

        [median,iqr] = validate(fcn,valpar,'AxisHandle',gca)

