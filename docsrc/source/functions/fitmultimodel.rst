.. highlight:: matlab
.. _fitmultimodel:


***********************
:mod:`fitmultimodel`
***********************

Automatic multi-component parametric model fitting of a distance distribution

-----------------------------


Syntax
=========================================

.. code-block:: matlab

    P = fitmultimodel(S,K,r,@dd_model,Nmax)
    P = fitmultimodel(S,t,r,@dd_model,Nmax)
    P = fitmultimodel(S,K,r,@dd_model,Nmax,metric)
    P = fitmultimodel(S,t,r,@dd_model,Nmax,metric,'Background',bg_model)
    P = fitmultimodel({V1,V2,___},{K1,K2,___},r,@dd_model,Nmax,metric)
    P = fitmultimodel({V1,V2,___},{t1,t2,___},r,@dd_model,Nmax,metric)
    P = fitmultimodel({V1,V2,___},{t1,t2,___},r,@dd_model,Nmax,metric,'Background',bg_model)
    P = fitmultimodel(___,'Property',Value)
    [P,param,Pci,paramci,Nopt,metrics,Peval] = fitmultimodel(___)


Parameters
    *   ``S`` - Input signal (*N*-element array)
    *   ``K`` -  Dipolar kernel (*NxM*-element array)
    *   ``r`` -  Distance Axis (*N*-element array)
    *   ``t`` -  Time Axis (*N*-element array)
    *   ``@dd_model`` -  Parametric model basis function (function handle)
    *   ``Nmax`` - Maximum number of components (scalar)
    *    ``metric`` - Metric for model selection (string)


Returns
    *  ``P`` - Fitted distance distribution (*M*-element array)
    *  ``param`` - Fitted model parameters (*W*-array)
    *  ``Pci`` - Fitted distribution confidence intervals (*Mx2*-array)
    *  ``paramci`` - Fitted parameters confidence intervals (*Wx2*-array)
    *  ``Nopt`` - Optimal number of components (scalar)
    *  ``metrics`` - Evaluated model selection functionals (cell array)
    *  ``Peval`` - Fitted distance distributions for each multi-component model (*Nmax x N* matrix)

-----------------------------


Description
=========================================

.. code-block:: matlab

        P = fitmultimodel(S,K,r,@dd_model,Nmax)

Fits the dipolar signal ``S`` to a distance distribution ``P`` using a multi-component parametric model according to the dipolar kernel ``K`` and distance axis ``r``. The multi-component model is constructed from the basis function provided as ``@dd_model``. The function chooses the optimal number of components distributions up to a maximum number given by ``Nmax`` by means of the corrected Akaike information criterion (AICC).

-----------------------------


.. code-block:: matlab

        P = fitmultimodel(S,K,r,@,dd_model,Nmax,metric)

The metric employed for the selection of the optimal number of components can be specified as an additional input ``metric``. The accepted inputs are:

	*   ``'aic'`` - Akaike information criterion
	*   ``'aicc'`` - Corrected Akaike information criterion
	*   ``'bic'`` - Bayesian information criterion

-----------------------------


.. code-block:: matlab

        P = fitmultimodel(S,t,r,@,dd_model,Nmax)

If the default kernel is to be used, the time-axis can be passed instead of the kernel.

-----------------------------


.. code-block:: matlab

	P = fitmultimodel(S,t,r,@dd_model,Nmax,metric,'Background',model)

By passing the ``'Background'`` option, the background function and modulation depth are fitted along the multi-component distribution parameters. 

-----------------------------


.. code-block:: matlab

    P = fitmultimodel({V1,V2,___},{K1,K2,___},r,@dd_model,Nmax,metric)

Passing multiple signals/kernels enables distance-domain global fitting of the parametric model to a single distribution. The global fit weights are automatically computed according to their contribution to ill-posedness. The multiple signals are passed as a cell array of arrays of sizes *N1*, *N2*,... and a cell array of kernel matrices with sizes *N1xM*, *N2xM*, ... must be passed as well.


-----------------------------

.. code-block:: matlab


    P = fitmultimodel({V1,V2,___},{t1,t2,___},r,@dd_model,Nmax,metric)
    P = fitmultimodel({V1,V2,___},{t1,t2,___},r,@dd_model,Nmax,metric,'Background',model)

Similarly, time-domain global fitting can be used when passing time-domain ``models`` and the model time axes ``{t1,t2,___}`` of the corresponding signals. If a background model is specified, it will be applied to all input signals. 



-----------------------------


.. code-block:: matlab

    [P,param,Nopt,metrics] = fitmultimodel(____)

If requested alongside the distribution ``P``, the optimal fit model parameters ``param``, as well their respective confidence intervals ``Pci`` and ``paramci`` the optimal number of components ``Nopt`` and evaluated selection metrics ``metrics`` are returned.

-----------------------------


Additional Settings
=========================================

Additional settings can be specified via name-value pairs. All property names are case insensitive and the property-value pairs can be passed in any order after the required input arguments have been passed.


.. code-block:: matlab

    P = fitmultimodel(___,'Property1',Value1,'Property2',Value2,___)

- ``'Background'`` - Parametric background model
    Function handle of the corresponding time-domain background model.

    *Default:* [*empty*] - Background and modulation depth are not fitted

    *Example:*

		.. code-block:: matlab

			P = fitmultimodel(___,'Background',@bg_exp)

- ``'Upper'`` - Parameters upper bound constraints
    An array of *W*-elements containing the upper bounds for the *W* parameters accepted by the model function ``@dd_model``.

    *Default:* [*empty*] - Uses the model's default upper bound values

    *Example:*

		.. code-block:: matlab

			P = fitmultimodel(___,@dd_gauss,___,'Upper',[rmean_max FWHM_max])

- ``'Lower'`` - Parameters lower bound constraints
    An array of *W*-elements containing the lower bounds for the *W* parameters accepted by the model function ``@dd_model``.

    *Default:* [*empty*] - Uses the model's default lower bound values

    *Example:*

		.. code-block:: matlab

			P = fitmultimodel(___,@dd_gauss,___,'Upper',[rmean_in FWHM_min])

- See :ref:`fitparamodel` for a detailed list of other property-value pairs accepted by the function.