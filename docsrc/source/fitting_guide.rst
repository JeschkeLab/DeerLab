Fitting Guide
=========================================

In this guide, we will focus on the last step of the DeerLab workflow, i.e., the fitting. Due to the separation of non-linear parameters `\theta_\mathrm{nonlin}` and linear parameters `\theta_\mathrm{lin}` in the mdoel structure of DeerLab, the program can use penalized separable non-linear least-squares to fit the model to the data. The aim of the fitting procedure is to find the following optimization problem with the objective function 

.. math::

    \min_{\theta_\mathrm{nonlin}} \left\{ \Vert y - A(\theta_\mathrm{nonlin})\theta_\mathrm{lin}(\theta_\mathrm{nonlin}) \Vert^2  + \sum_i \gamma_i\mathcal{P}_i(\theta_\mathrm{nonlin},\theta_\mathrm{lin}) \right\} \\

with the linear part of the objective given by

.. math::

    \theta_\mathrm{lin}(\theta_\mathrm{nonlin}) = {\arg\!\min}_{\theta} \left\{ \Vert y - A(\theta_\mathrm{nonlin})\theta \Vert^2 +  \alpha^2\mathcal{R}(\theta) \right\} \\

and all subject to the boundary conditions defined in the model   

.. math::

    \theta_\mathrm{lb} \leq \theta_\mathrm{nonlin} \leq \theta_\mathrm{ub} \quad\quad \theta_\mathrm{lbl} \leq \theta_\mathrm{lin} \leq \theta_\mathrm{ubl} 


The term `A(\theta_\mathrm{nonlin})` is the non-linear function of the model. A regularization penalty `\mathcal{R}(\theta_\mathrm{lin})` weighted by a regularization weight or parameter `\alpha`, and multiple additional arbitrary penalty terms `\mathcal{P}_i(\theta_\mathrm{nonlin},\theta_\mathrm{lin})` each weighted by its own penalty weight `\gamma_i`, can be added to this objective function if wanted


The ``fit`` function
--------------------

DeerLab provides a centralized fitting function called ``fit``. Its syntax is generally simple, taking a dataset and a model to be fitted to it. The function will always return a ``FitResult`` object containing all of the fitted quantities including uncertainty estimates, as well as other quantities of interest. Assume that we have some dataset ``y`` described by some model called ``model``. The fit can be executed as follows :: 

    # Fit the model to the data
    result = dl.fit(model, y) 

If the model has parameters without assigned start values, the ``fit`` function will return an error and request them to be specified. This can be done via the ``par0`` keyword argument passing the start values as a list ::

    # Fit the model to the data with other start values 
    result = dl.fit(model, y, par0=par0_list)

However, the list must be ordered according to the model parameter ordering. Therefore, it is recommended to specify the start values on the model (as shown :ref:`here <modelling_modifying_parameters>`).  

Fitting models with constants
*****************************

Models with constants (see :ref:`here <modelling_constants>` for details) can be fitted as shown above, with the addeded requirement that the constants must be specified as well when calling ``fit``. Constants can be specified after the data as positional arguments :: 

    # Fit the model (with two constants) to the data
    result = dl.fit(model, y, constant1, constant2) 


Fitting multi-dataset models
****************************

Models that have been merged using the ``merge`` function (see :ref:`here <modelling_merging>` for details) can describe multiple datasets with a single model object and a common parameter set. To fit such a merged model to multiple datasets, the ``fit`` function can be used as above by passing a list of datasets ``[y1,y2,...,yN]`` instead of a single dataset  ::

    # Fit the model to multiple datasets
    result = dl.fit(model, [y1,y2,y3]

The number of datasets must match the number of reponses returned by the model. Additionally, the ordering in the list of datasets must match the order of responses from the model, i.e. ``response1`` of ``model`` will be fitted to ``y1``, etc. 

Adding regularization
---------------------

DeerLab includes the possibility to impose (Tikhonov) regularization based on smoothness of the linear parameters in the model. By default, DeerLab will automatically check the condition number of the matrix returned by the non-linear of the model and determine whether it is well-conditioned or not. If the matrix is ill-conditioned, the program will automatically enforce regularization upon the linear parameters. 
Regularization can be manually be switched on/off via the ``reg`` param :: 

    # Enforce regularization 
    result = dl.fit(model, y, reg=True)
    # Disable regularization 
    result = dl.fit(model, y, reg=False)

The regularization penalty weight (a.k.a regularization parameter) is optimally selected according to a given criterion (by default the Akaike information criterion, AIC). There are different ways to control this process: 

Changing the selection functional 
*********************************

The regularization functional can be changed from the AIC to any other of the built-in functionals via the `regparam` keyword argument. Changing it to another functional will only change how the regularization parameter is selected, but it still will be optimized. For example to switch the selection functional from the AIC to generalized cross-validation (GCV) ::

    # Fit the model to the data, using the GCV criterion 
    result = dl.fit(model,y, regparam='gcv')

A list of the available selection functionals and their string names are given in the reference for ``fit``.

Changing the optimization range
********************************

When using selection functionals to optimize the regularization weights, a Brent-like algorithm is used to search the value which minimizes the given selection functional within a certain range. This range can be manually specified via the ``regparamrange`` keyword argument. It must be passed as a two-element list ``[regparam_lb, regparam_ub]`` with the search boundaries ::

    # Fit the model to the data, with a constrained regparam search range 
    result = dl.fit(model,y, regparamrange=[1e-5,1e-1])

This can be useful for avoiding unwanted local minima of the selection functional causing potential under- or oversmoothing.

Manual specification
**********************

The value of the regularization penalty weight can also be manually specified and fixed to a value for the whole optimization. This can be done via the aforementioned ``regparam`` keyword argument by specifying a value instead of a selection functional :: 

    # Fit the model to the data, using a fixed regularization weight 
    result = dl.fit(model,y, regparam=0.05) 


Adding penalties
----------------

DeerLab provides a flexible system for defining and adding penalties to the objective function of the ``fit`` function in the form of the ``Penalty`` objects which can be passed to the ``penalties`` keyword argument of the ``fit`` function ::

    # Fit the model to the data with an additional penalty
    result = dl.fit(model,y penalties=penalty)

Penalties are only added to the non-linear part of the separable non-linear least-squares objective function used in ``fit``. For the linear part, only Tikhonov regularization can be imposed (see previous section). 


DeerLab's penalties consist of the following components: 

Penalty function 
    The penalty function takes model parameters and returns a vector of values which are appended to the least-squares residual. The function should ideally be convex, monotonically increasing and defined everywhere. Can be freely constructed and defined. 

Penalty weight 
    As its name indicates, the penalty weight balances the influence of the penalty with respect to the other terms in the objective function. It is treated similarly to model parameters, meaning that it has boundaries defined which can be manipulated freely. 

Selection functional    
    The selection criterion desired for the optimized choice of penalty weight. Must be chosen from a collection of selection functionals.  



Constructing a penalty
**********************



