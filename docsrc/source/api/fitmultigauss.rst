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

    P = fitmultigauss(S,K,r,Nmax)
    P = fitmultigauss(S,t,r,Nmax)
    P = fitmultigauss(S,K,r,Nmax,metric)
    P = fitmultigauss(S,t,r,Nmax,metric,'Background',model)
    P = fitmultigauss(___,'Property',Value)
    [P,param,Pci,paramci,Nopt,metrics,Peval] = fitmultigauss(___)


Parameters
    *   ``S`` - Input signal (*N*-element array)
    *   ``K`` -  Dipolar kernel (*NxM*-element array)
    *   ``r`` -  Distance Axis (*N*-element array)
    *   ``t`` -  Time Axis (*N*-element array)
    *   ``Nmax`` - Maximum number of Gaussians (scalar)
    *    ``metric`` - Metric for model selection (string)


Returns
    *  ``P`` - Fitted distance distribution (*M*-element array)
    *  ``param`` - Fitted model parameters (*W*-array)
    *  ``Pci`` - Fitted distribution confidence intervals (*Mx2*-array)
    *  ``paramci`` - Fitted parameters confidence intervals(*Wx2*-array)
    *  ``Nopt`` - Optimal number of Gaussian (scalar)
    *  ``metrics`` - Evaluated model selection functionals (cell array)
    *  ``Peval`` - Fitted distance distributions for each multi-gauss model (*Nmax x N* matrix)

-----------------------------


Description
=========================================

.. code-block:: matlab

        P = fitmultigauss(S,K,r,Nmax)

Fits the dipolar signal ``S`` to a distance distribution ``P`` using a multi-gauss parametric model according to the dipolar kernel ``K`` and distance axis ``r``. The function chooses the optimal number of Gaussian distributions up to a maximum number given by ``Nmax`` by means of the corrected Akaike information criterion (AICC).

-----------------------------


.. code-block:: matlab

        P = fitmultigauss(S,K,r,Nmax,metric)

The metric employed for the selection of the optimal multi-gauss model can be specified as an additional input ``metric``. The accepted inputs are:

	*   ``'aic'`` - Akaike information criterion
	*   ``'aicc'`` - Corrected Akaike information criterion
	*   ``'bic'`` - Bayesian information criterion

-----------------------------


.. code-block:: matlab

        P = fitmultigauss(S,t,r,Nmax)

If the default kernel is to be used, the time-axis can be passed instead of the kernel.

-----------------------------


.. code-block:: matlab

	P = fitmultigauss(S,t,r,Nmax,metric,'Background',model)

By passing the ``'Background'`` option, the background function and modulation depth are fitted along the multi-gauss distribution parameters. 

-----------------------------


.. code-block:: matlab

    [P,param,Nopt,metrics] = fitmultigauss(args)

If requested alongside the distribution ``P``, the optimal fit model parameters ``param``, as well their respective confidence intervals ``Pci`` and ``paramci`` the optimal number of Gaussians ``Nopt`` and evaluated selection metrics ``metrics`` are returned.

-----------------------------


Additional Settings
=========================================

Additional settings can be specified via name-value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.


.. code-block:: matlab

    P = fitmultigauss(___,'Property1',Value1,'Property2',Value2,___)

- ``'Background'`` - Parametric background model
    Function handle of the corresponding time-domain background model.

    *Default:* [*empty*] - Background and modulation depth are not fitted

    *Example:*

		.. code-block:: matlab

			P = fitmultigauss(___,'Background',@bg_exp)

- ``'Upper'`` - Parameters upper bound constraints
    Array ``[<r>_max FWHM_max]`` containing the upper bound for the FWHM and mean distance of all the Gaussians.

    *Default:* [*empty*] - Uses the model's default upper bound values

    *Example:*

		.. code-block:: matlab

			P = fitmultigauss(___,'Upper',[10 0.9])

- ``'Lower'`` - Parameters lower bound constraints
    Array ``[<r>_min FWHM_min]`` containing the lower bound for the FWHM and mean distance of all the Gaussians.

    *Default:* [*empty*] - Uses the model's default lower bound values

    *Example:*

		.. code-block:: matlab

			P = fitmultigauss(___,'Lower',[1 0.1])

- See :ref:`fitparamodel` for a detailed list of other property-value pairs accepted by the function.