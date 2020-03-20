.. highlight:: matlab
.. _fitmultigauss:


***********************
:mod:`fitmultigauss`
***********************

Multi-Gauss fitting of a distance distribution

-----------------------------


Syntax
=========================================

.. code-block:: matlab

    P = fitmultigauss(S,K,r,Ngauss)
    P = fitmultigauss(S,t,r,Ngauss)
    P = fitmultigauss(S,K,r,Ngauss,metric)
    P = fitmultigauss(S,t,r,Ngauss,metric,'Background',model)
    [P,param,Nopt,metrics,Peval] = fitmultigauss(S,K,r,Ngauss,metric,'Property',Value)


Parameters
    *   ``S`` - Input signal (*N*-element array)
    *   ``K`` -  Dipolar kernel (*NxM*-element array)
    *   ``r`` -  Distance Axis (*N*-element array)
    *   ``t`` -  Time Axis (*N*-element array)
    *   ``Ngauss`` - Maximum number of Gaussians (scalar)
    *    ``metric`` - Metric for model selection (string)


Returns
    *  ``P`` - Distance Distribution (*M*-element array)
    *  ``param`` - Fitted model parameters (array)
    *  ``Nopt`` - Optimal number of Gaussian (scalar)
    *  ``metrics`` - Evaluated model selection functionals (cell array)
    *  ``Peval`` - Fitted distance distributions for each multigauss model (*Ngauss x N* matrix)

-----------------------------


Description
=========================================

.. code-block:: matlab

        P = fitmultigauss(S,K,r,Ngauss)

Fits the dipolar signal ``S`` to a distance distribution ``P`` using a multi-gauss parametric model according to the dipolar kernel ``K`` and distance axis ``r``. The function chooses the optimal number of Gaussian distributions up to a maximum number given by ``Ngauss`` by means of the corrected Aikaike information criterion (AICC).

-----------------------------


.. code-block:: matlab

        P = fitmultigauss(S,K,r,Ngauss,metric)

The metric employed for the selection of the optimal multigauss model can be specified as an additional input ``metric``. The accepted inputs are:

	*   ``'aic'`` - Akaike information criterion
	*   ``'aicc'`` - Corrected Akaike information criterion
	*   ``'bic'`` - Bayesian information criterion

-----------------------------


.. code-block:: matlab

        P = fitmultigauss(S,t,r,Ngauss)

If the default kernel is to be used, the time-axis can be passed instead of the kernel.

-----------------------------


.. code-block:: matlab

	P = fitmultigauss(S,t,r,Ngauss,metric,'Background',model)

By passing the ``'Background'`` option, the background function and modulation depth are fitted along the multigauss distribution parameters. 

-----------------------------


.. code-block:: matlab

    [P,param,Nopt,metrics] = fitmultigauss(args)

If requested alongside the distribution ``P``, the optimal fit model parameters ``param``, the optimal number of gaussians ``Nopt`` and evaluated selection metrics ``metrics`` are returned.

-----------------------------


Optional Arguments
=========================================
Optional arguments can be specified by parameter/value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.

.. code-block:: matlab

    P = fitmultigauss(args,'Property1',Value1,'Property2',Value2,...)

- ``'Background'`` - Parametric background model
    Function handle of the corresponding time-domain background model.

    *Default:* [*empty*] - Background and modulation depth are not fitted

    *Example:*

		.. code-block:: matlab

			P = fitmultigauss(S,t,r,Ngauss,metric,'Background',@bg_exp)

- ``'Upper'`` - Parameters upper bound constraints
    Array ``[<r>_max FWHM_max]`` containing the upper bound for the FWHM and mean distance of all the Gaussians.

    *Default:* [*empty*] - Uses the model's default upper bound values

    *Example:*

		.. code-block:: matlab

			P = fitmultigauss(arg,'Upper',[10 0.9])

- ``'Lower'`` - Parameters lower bound constraints
    Array ``[<r>_min FWHM_min]`` containing the lower bound for the FWHM and mean distance of all the Gaussians.

    *Default:* [*empty*] - Uses the model's default lower bound values

    *Example:*

		.. code-block:: matlab

			P = fitmultigauss(arg,'Lower',[1 0.1])

- See :ref:`fitparamodel` for a detailed list of other property-value pairs accepted by the function.