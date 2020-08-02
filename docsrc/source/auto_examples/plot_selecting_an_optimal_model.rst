.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_plot_selecting_an_optimal_model.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_auto_examples_plot_selecting_an_optimal_model.py:


Selecting an optimal parametric model for fitting a dipolar signal
==================================================================

How to optimally select a parametric model for a given dipolar signal.


.. code-block:: python


    import numpy as np
    import matplotlib.pyplot as plt
    from deerlab import *








Data Generation
----------------

Let's start by constructing a simple dipolar signal with some noise arising 
from a bimodal Gaussian distance distribution.


.. code-block:: python


    # Prepare the signal components
    t = np.linspace(-0.3,6,300)
    r = np.linspace(2,6,200)
    P = dd_gauss2(r,[3.8, 0.7, 0.7, 4.5, 0.3, 0.7])

    # Prepare the dipolar kernel and get the signal
    K = dipolarkernel(t,r)
    V = K@P + whitegaussnoise(t,0.01)








Selecting an optimal model
--------------------------

Even though we know the ground truth, in this example we will cosider the 
following set of potential parametric models: 

* Unimodal Rician distribution
* Bimodal Rician distribution
* Trimodal Rician distribution
* Unimodal Gaussian distribution
* Bimodal Gaussian distribution
* Trimodal Gaussian distribution
* Mixed bimodal Gaussian/Rician distribution

The first six models have built-in parametric models which we can use directly. 
The last model we can construct from built-in models using the ``mixmodels`` function.


.. code-block:: python


    # Prepare the mixed model
    dd_rice_gauss = mixmodels(dd_rice,dd_gauss)
 
    # Prepare list of candidate parametric models
    models = [dd_rice,dd_rice2,dd_rice3,dd_gauss,dd_gauss2,dd_gauss3,dd_rice_gauss]








In order to make an appropiate choice, we need some liklihood estimator. All fit functions is DeerLab returns a stats 
dictionary which contains (amongst other estimators) likelihood estimators such as the Akaike information criterion (AIC).
The model with the lowers AIC value can be considered to most likely to be the optimal model.

To do this, we jus have to evaluate the parametric models with ``fitparamodel`` while looping over all the distribution models
we listed above, and collecting the AIC-values for each model.


.. code-block:: python

 
    aic = []
    for model in models:
        info = model()
        # Prepare the signal model with the new distance model
        Vmodel = lambda par: K@model(r,par)
        # Fit the signal
        fit = fitparamodel(V,Vmodel,par0=info['Start'],lb=info['Lower'],ub=info['Upper'])
        parfit = fit.param
        stats= fit.stats
        # Add current AIC value to the list
        aic.append(stats['aic'])





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    d:\lufa\projects\deerlab\deerlab\deerlab\fitparamodel.py:191: UserWarning: The fitted value of parameter #5, is at the lower bound of the range.
      warnings.warn('The fitted value of parameter #{}, is at the lower bound of the range.'.format(p))
    d:\lufa\projects\deerlab\deerlab\deerlab\utils\gof.py:54: RuntimeWarning: divide by zero encountered in double_scalars
      R2 = 1 - np.sum((x-xfit)**2)/np.sum((xfit-np.mean(xfit))**2)
    d:\lufa\projects\deerlab\deerlab\deerlab\fitparamodel.py:191: UserWarning: The fitted value of parameter #8, is at the lower bound of the range.
      warnings.warn('The fitted value of parameter #{}, is at the lower bound of the range.'.format(p))
    d:\lufa\projects\deerlab\deerlab\deerlab\utils\gof.py:54: RuntimeWarning: divide by zero encountered in double_scalars
      R2 = 1 - np.sum((x-xfit)**2)/np.sum((xfit-np.mean(xfit))**2)
    d:\lufa\projects\deerlab\deerlab\deerlab\fitparamodel.py:191: UserWarning: The fitted value of parameter #2, is at the lower bound of the range.
      warnings.warn('The fitted value of parameter #{}, is at the lower bound of the range.'.format(p))
    d:\lufa\projects\deerlab\deerlab\deerlab\fitparamodel.py:191: UserWarning: The fitted value of parameter #5, is at the lower bound of the range.
      warnings.warn('The fitted value of parameter #{}, is at the lower bound of the range.'.format(p))




Since the absolute AIC values have no meaning, it is standard practice to look at the relative 
changes in AIC values between the evaluated models.


.. code-block:: python


    daic = aic - min(aic)








Akaike Weights
-----------------------------------------------------------------------------
 It is often more useful to look at these results from the perspective of
 Akaike weights, i.e. the probabilities of a model being the most optimal.


.. code-block:: python


    weights = 100*np.exp(-(daic/2))/sum(np.exp(-daic/2))








Plot results
------------


.. code-block:: python


    plt.figure(figsize=(9,8))

    plt.subplot(2,2,1)
    plt.plot(t,V,'k.')
    plt.grid(alpha=0.2)
    plt.xlabel('t [$\mu s$]')
    plt.ylabel('V(t)')
    plt.legend(['data'])

    plt.subplot(2,2,2)
    plt.plot(r,P,'k',linewidth=1.5)
    plt.xlabel('r [nm]')
    plt.ylabel('P(r) [nm$^{-1}$]')
    plt.legend(['Ground truth'])
    plt.grid(alpha=0.2)

    modelnames = [model.__name__ for model in models]

    plt.subplot(2,2,3)
    plt.bar(modelnames,daic,color='b',alpha=0.5)
    plt.ylabel('$\Delta$AIC')
    plt.grid(alpha=0.2)
    plt.xticks(rotation=45)

    # Plot the results
    plt.subplot(2,2,4)
    plt.bar(modelnames,weights,color='b',alpha=0.5)
    plt.ylabel('Akaike Weights [%]')
    plt.xticks(rotation=45)
    plt.grid(alpha=0.2)




.. image:: /auto_examples/images/sphx_glr_plot_selecting_an_optimal_model_001.png
    :alt: plot selecting an optimal model
    :class: sphx-glr-single-img





Typically there is not a single optimal model unless the noise level is very
low. Usually several models have similar probabilities and should therefore be presented together. 


.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  13.259 seconds)


.. _sphx_glr_download_auto_examples_plot_selecting_an_optimal_model.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_selecting_an_optimal_model.py <plot_selecting_an_optimal_model.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_selecting_an_optimal_model.ipynb <plot_selecting_an_optimal_model.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
