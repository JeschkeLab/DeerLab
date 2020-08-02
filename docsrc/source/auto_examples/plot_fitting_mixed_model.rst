.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_plot_fitting_mixed_model.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_auto_examples_plot_fitting_mixed_model.py:


Fitting a mixed distance-distribution model
===========================================

Basic manipulation of parametric models and creating mixed models 
for fitting distance distributions.


.. code-block:: python


    import numpy as np
    import matplotlib.pyplot as plt
    from deerlab import *








Simulate the data
-----------------

Let's start by creating a simple dipolar evolution function (i.e. no background 
and full modulation depth) corresponding to a simple 4-pulse DEER signal.


.. code-block:: python


    #Axis definition
    t = np.linspace(-0.5,4,350)
    r = np.linspace(2,6,200)

    # Distribution parameters
    rmean = 4.5
    width = 0.3
    chain = 4.3
    pers = 10
    amp = 0.35

    # Generate distribution
    P = dd_gauss(r,[rmean, width])
    P = amp*P + (1 - amp)*dd_wormchain(r,[chain, pers])
    # Normalize distribution
    P = P/sum(P)/np.mean(np.diff(r))
    # Generate dipolar evolution function
    np.random.seed(0)
    K = dipolarkernel(t,r)
    V = K@P + whitegaussnoise(t,0.02)








Generating a mixed parametric model
-----------------------------------

Let's say our intuiton (which, since we know the ground truth, is exact) on 
the sample indicates that our distribution is a llinear combination of a Gaussian 
distirbution and a worm-like chain model. While DeerAnalysis provides built-in 
parametric models for both models, we require a combination of both. 

For such cases we can use the ``mixmodels`` function to create a custom mixed 
parametric model. It's syntax is rather simple, we just have to pass the desired 
parametric models as lambda functions. 


.. code-block:: python


    #Mix the models into new one
    gausswlc = mixmodels(dd_gauss,dd_wormchain)








Our new model ``gausswlc`` will now describe our sought linear combination of 
both parametric models. We can check the state of the model by retrieving its 
information


.. code-block:: python


    #Get information on the mixed model
    info = gausswlc()








We can see that the ``mixmodels`` function has introduced an ampitude parameters 
as the first parameter of the model. This parameters weights the contribution 
of each individual parametric model. We see also that this is followed by the 
parameters of the Gaussian model and finally with the parameters of the WLC 
model.

Our model is ready, and since it was generated from built-in models we do 
not need to specify any parameters initial values or boundary constraints. These 
can, however, by re-defined if the built-in defaults are not appropiate (see 
other examples). 

Since we are dealing with a distance-domain model we require a dipolar kernel 
to transform our model into time-domain. Remember that our signal in this example 
is a dipolar evolution function, therefore we do not require anything else than 
a very basic dipolar kernel.


.. code-block:: python


    #Generate the dipolar evolution function kernel
    K = dipolarkernel(t,r)

    #Fit the model to the data
    Vmodel = lambda par: K@gausswlc(r,par)
    info = gausswlc()
    par0 = info['Start'] # built-in start values
    lb = info['Lower'] # built-in lower bounds
    ub = info['Upper'] # built-in upper bounds
    fit = fitparamodel(V,Vmodel,par0,lb,ub,MultiStart=10)
    fitpar = fit.param







From the fitted parameter set ``fitpar`` we can now generate our fitted distance 
distribution and the corresponding time-domain fit.


.. code-block:: python


    #Calculate the fitted model
    Pfit = gausswlc(r,fitpar)
    Vfit = Vmodel(fitpar)








Since we know both the ground truth for the distance distribution and the 
dipolar signal, let's see how our fit turned out.


.. code-block:: python


    #Plot results
    plt.subplot(2,1,1)
    plt.plot(t,V,'k.',t,Vfit,'r',linewidth=1.5)
    plt.xlabel('t [$\mu s$]')
    plt.ylabel('V(t)')
    plt.legend(['data','fit'])

    plt.subplot(2,1,2)
    plt.plot(r,P,'k',r,Pfit,'r',linewidth=1.5)
    plt.xlabel('r [nm]')
    plt.ylabel('P(r) [nm$^{-1}$]')
    plt.legend(['truth','fit'])





.. image:: /auto_examples/images/sphx_glr_plot_fitting_mixed_model_001.png
    :alt: plot fitting mixed model
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    <matplotlib.legend.Legend object at 0x000001ED98887518>




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  26.163 seconds)


.. _sphx_glr_download_auto_examples_plot_fitting_mixed_model.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_fitting_mixed_model.py <plot_fitting_mixed_model.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_fitting_mixed_model.ipynb <plot_fitting_mixed_model.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
