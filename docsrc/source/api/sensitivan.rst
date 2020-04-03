.. highlight:: matlab
.. _sensitivan:

***********************
:mod:`sensitivan`
***********************

Sensitivity analysis via factorial design


------------------------


Syntax
=========================================

.. code-block:: matlab

    stats = sensitivan(fcn,varpar)
    [stats,factors,evals] = sensitivan(fcn,varpar)
    [stats,factors,evals] = sensitivan(___,'Property',Value)

Parameters
    *   ``fcn`` - Function to validate (function handle)
    *   ``varpar`` - Parameter variation structure (struct)

Returns
    *   ``stats`` - Summary statistics (struct)

         *   ``.median`` - Median of output the variables
         *   ``.mean`` - Mean of the output variables
         *   ``.std`` - Standard deviation of the output variables
         *   ``.p2`` - 2nd percentile of the output variables
         *   ``.p25`` - 25th percentile of the output variables
         *   ``.p75`` - 75th percentile of the output variables
         *   ``.p98`` - 98th percentile of the output variables


    *   ``factors`` - Factor analysis results (struct)

         *   ``.main`` - Main effects of the factors
         *   ``.inter`` - interactions between factors
		 
		 
    *   ``evals`` - Evaluated function output arguments (cell array)


------------------------


Description
=========================================

.. code-block:: matlab

    stats = sensitivan(fcn,varpar)

Performs a sensitivity analysis of the output variables of the function ``fcn`` with respect to the parameter ranges defined in ``varpar``. The output argument of ``fcn`` is evaluated for all level-combinations of the factors. The function to be validated must be a function handle accepting the ``varpar`` struct as its first argument and only argument. 

.. code-block:: matlab

    stats = sensitivan(@(p)myvalfcn(p,varargin),varpar)

    function [out1,out2] = myvalfcn(varpar,varargin)
       out1 = process(varpar.param1,varpar.param2)
       out2 = process2(varpar.param1,varpar.param3)
    end


From the evaluation of all level-combinations an ensemble of outputs is obtained on which statistical estimators are used. The summary of these statistics is returned in the ``stats`` structure. This summary contains the mean, median, standard deviation, 98t, 75th, 25th and 2nd percentile values for all outputs. 

The range defined by the 25th and 75th percentiles is known as the inter-quartile range and contains 50% of all cases around the median. The range between the 2nd and 98th percentiles contains 99% of all cases around the median.


------------------------


.. code-block:: matlab

    [stats,factors] = sensitivan(fcn,varpar)
	

For each variable/factor ``varpar.param1``, ``varpar.param2``,... evaluated in the sensitivity analysis its main effect and interaction with the other factors will be returned in the ``factors`` output structure.


------------------------


.. code-block:: matlab

    [stats,factors,evals] = sensitivan(fcn,varpar)

Additionally, a last output argument ``evals`` can be requested, a cell array, containing the ``fcn`` outputs evaluated at each parameter combination.


------------------------


Additional Settings
=========================================

Additional settings can be specified via name-value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.



.. code-block:: matlab

    [median,iqr] = sensitivan(___,'Property1',Value1,'Property2',Value2,___)

- ``'RandPerm'`` - Randomized level-combination evaluation
    Specifies whether to randomly permute the sensitivity anaysis parameter combinations.

    *Default:* ``true``

    *Example:*

		.. code-block:: matlab

			[median,iqr] = sensitivan(___,'RandPerm',false)

- ``'AxisHandle'`` - Plot intermediate results
    Axis handle to plot the state of the validation results at each level combination. The display is refreshed every 5 level cobinations (for the sake of speed).

    *Default:* [*empty*]

    *Example:*

		.. code-block:: matlab

			[median,iqr] = sensitivan(___,'AxisHandle',gca)


- ``'dynamicStats'`` -  Enable dynamic statistical estimators
	Specifies whether the statistical estimators are computed using the full set of observations (``false``) or approximated dynamically at each iteration (``true``).

    *Default:* ``false`` (if enough memory available) ``true`` (if not enough memory available)

    *Example:*

		.. code-block:: matlab

			[median,iqr] = sensitivan(___,'dynamicStats',true)

