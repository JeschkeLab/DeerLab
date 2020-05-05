.. highlight:: matlab
.. _bootan:

***********************
:mod:`bootan`
***********************

Bootstrap analysis for uncertainty estimation

------------------------


Syntax
=========================================

.. code-block:: matlab

    [bootci,stats] = bootan(fcn,V,Vfit)
    [bootci,stats] = bootan(fcn,V,Vfit,Nsamples)
    [bootci,stats] = bootan(___,'Property',Value)

Parameters
    *   ``fcn`` - Function to analyze (function handle)
    *   ``V`` - Experimental signal (*N*-element array)
    *   ``Vfit`` - Fitted signal (*N*-element array)
    *   ``Nsamples`` - Number of bootstrap samples (scalar)

Returns
    *   ``bootci`` - Bootstrapped confidence intervals (struct)
         *   ``.ci50`` - 50%-confidence intervals of the output variables
         *   ``.ci95`` - 95%-confidence intervals of the output variables
         *   ``.ci99`` - 99%-confidence intervals of the output variables

    *   ``stats`` - Summary statistics (struct)

         *   ``.median`` - Median of output the variables
         *   ``.mean`` - Mean of the output variables
         *   ``.std`` - Standard deviation of the output variables
         *   ``.p1``  - 1st percentile of the output variables
         *   ``.p25`` - 25th percentile of the output variables
         *   ``.p75`` - 75th percentile of the output variables
         *   ``.p99`` - 99th percentile of the output variables
         *   ``.boothist`` - Bootstrap histogram of the output variables

             *   ``.edges`` - Histogram edges
             *   ``.bins`` - Histogram bins

         *   ``.bootdist`` - Bootstrap KDE distribution of the output variables

             *   ``.values`` - Evaluated values
             *   ``.pdf`` - Estimated probability density function


------------------------


Description
=========================================

.. code-block:: matlab

    [bootci,stats] = bootan(fcn,V,Vfit)

Performs a uncertainty analysis of the output variables of the function ``fcn`` from 1000 bootstrap samples. The output argument of ``fcn`` is evaluated for all level-combinations of the factors. The function to be analyzed must be a function handle accepting the ``V`` experimental signal as input. Example:

.. code-block:: matlab

    [bootci,stats] = bootan(@(V)myfcn(p,varargin),V,Vfit)

    function [Pfit1,Pfit2] = myfcn(V,varargin)
       Pfit1 = fitparamodel(V,@dd_gauss,r,K)
       Pfit2 = fitparamodel(V,@dd_randcoil,r,K)
    end


From the evaluation of all level-combinations an ensemble of outputs is obtained on which statistical estimators are used. The 99%, 95% and 50% confidence intervals of all output variables are returned in a structure ``bootci``. The summary of these statistics is returned in the ``stats`` structure. This summary contains the mean, median, standard deviation, 99t, 75th, 25th and 1st percentile values for all outputs.

For non-vectorial variables (e.g. parameter-free distributions, background functions,etc.) the ``stats`` structure will contain an histogram of the distribution of values for the different outputs as well as a corresponding probability density function obtained from a kernel densitiy estimation of the histogram data.

------------------------

.. code-block:: matlab

    stats = bootan(fcn,V,Vfit,Nsamples)


The number of bootstrap samples can be specified in ``Nsamples``. The quality of bootstrapping results improve with the number of boostrap samples evaluated. 



------------------------


.. code-block:: matlab

    stats = bootan(fcn,{V1,V2,___},{Vfit1,Vfit2,___},Nsamples)


If the evaluated function ``fcn`` requries multiple signals ``{V1,V2,___}`` as input, these can be specified aloong the same number of fitted signals ``{Vfit1,Vfit2,___}``. 


------------------------


Additional Settings
=========================================

Additional settings can be specified via name-value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.



.. code-block:: matlab

    stats = bootan(___,'Property1',Value1,'Property2',Value2,___)

- ``'Verbose'`` - Display progress information
    Specifies whether to print the progress of the bootstrap analysis on the command window.

    *Default:* ``false``

    *Example:*

		.. code-block:: matlab

			stats = bootan(___,'Verbose',true)


- ``'Resampling'`` - Re-sampling method
    Specifies the method employed for re-sampling new bootstrap samples.

        *   ``'gaussian'`` - Sample noise from a Gaussian distribution
        *   ``'residual'`` - Sample noise from the fit residuals

    *Default:* ``gaussian``

    *Example:*

		.. code-block:: matlab

			stats = bootan(___,'Resampling',residual)

