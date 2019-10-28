.. highlight:: matlab
.. _validate:

***********************
:mod:`validate`
***********************

Statistical validation of script/function parameters.

Syntax
=========================================

.. code-block:: matlab

    median = validate(fcn,varpar)
    [median,iqr,evals] = validate(fcn,varpar)
    [median,iqr,evals] = validate(fcn,varpar,'Property',Value)

Parameters
    *   ``fcn`` - Function to validate (function handle)
    *   ``varpar`` - Validation parameters (struct)

Returns
    *   ``median`` - Validation median value (cell array)
    *   ``std`` - Validation inter-quartile range (cell array)
    *   ``evals`` - Evaluated function output arguments (cell array)


Description
=========================================

.. code-block:: matlab

    [median,iqr] = validate(fcn,varpar)

Performs a sensitivity analysis of the output of the input function ``fcn`` with respect to the parameter ranges defined in ``varpar``. The output argument of ``fcn`` is evaluated for all combinations of the validation parameters. The function to be validated must be a function handle accepting the ``varpar`` struct as its first argument. 

.. code-block:: matlab

    function [out1,out2] = myvalfcn(varpar,varargin)
       out1 = process(varpar.param1,varpar.param2)
       out2 = process2(varpar.param1,varpar.param3)
    end


The median and inter-quartile range (IQR) of the resulting statistics are returned as the ``median`` and ``iqr`` outputs. If the function ``fcn`` returns multiple arguments, all of them are validated and ``median`` and ``iqr`` are returned as cell arrays with a median and IQR-value for each output.

.. code-block:: matlab

    [median,iqr,evals] = validate(fcn,varpar)

Additionally, a third output argument ``evals`` can be requested, a cell array, containing the ``fcn`` outputs evaluated at each parameter combination.

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

