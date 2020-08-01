.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_plot_param_uncertainty_propagation.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_auto_examples_plot_param_uncertainty_propagation.py:


Uncertainty propagation from parameter fits using covariance-based uncertainty quantificaion
===========================================================================================

How to propagate the uncertainty of the fitted parameters to the models which depend on them.


.. code-block:: python


    import numpy as np
    import matplotlib.pyplot as plt
    from deerlab import *








Generate data
-------------


.. code-block:: python


    t = np.linspace(-0.2,4,300)
    r = np.linspace(2,5,400)
    center = 3.5 # [nm] Rician center distance
    width = 0.3 # [nm] Rician width
    lam = 0.27 # Modulation depth
    conc = 150 # [uM] Spin concentration
    P = dd_rice(r,[center, width])
    B = bg_hom3d(t,conc,lam)
    K = dipolarkernel(t,r,lam,B)
    V = K@P + whitegaussnoise(t,0.03)








Fit the data
------------
First we define the models for the different elements in our analysis
(background, distribution and dipolar signal). For simplicity these
models take the full parameter set

``par = [lambda center width conc]``

and select the appropiate elements from the parameter set, i.e.

``Pmodel = f(center,width) -> par[1] & par[2]``
``Bmodel = f(conc,lambda)  -> par[3] & par[0]``
``Vmodel = f(par)          -> par[0] & par[1] & par[2] & par[3]``

By defining the models like this, we can spare then the indexing of the
parameters each time we call one of these model and can pass the full
parameter set directly.


.. code-block:: python


    # Pre-calculate the elemental dipolar kernel (for speed)
    K0 = dipolarkernel(t,r)

    Pmodel = lambda par: dd_rice(r,par[1:3])
    Bmodel = lambda par: bg_hom3d(t,par[3],par[0])
    Vmodel = lambda par: (1 - par[0] + par[0]*K0@Pmodel(par))*Bmodel(par)








Next since we are dealing with a custom-defined model we need to specify
the start values as well as boundaries of the parameter set:


.. code-block:: python


    # Parameters:[lam center width conc]
    par0  =      [0.35, 4.0,  0.4, 500 ] # start values
    lower =      [0.10, 2.0,  0.1, 0.1 ] # lower bounds
    upper =      [0.50, 7.0,  0.5, 1500] # upper bounds

    # Finally we can run the fit and get the fitted parameters and their uncertainties
    parfit,paruq,_ = fitparamodel(V,Vmodel,par0,lower,upper)

    # Forward-calculate the models with the fitted parameters
    Vfit = Vmodel(parfit)
    Pfit = Pmodel(parfit)
    Bfit = Bmodel(parfit)
    lamfit = parfit[0]








Uncertainty propagation
------------------------
 In DeerLab, all uncertainty quantification objects contain a method
 ``.propagate()``, which has all the internal information on the 
 covariance matrices required to propagate the uncertainty from 
 the parameters to the fits. 

 Thus, all we neeed to do is call ``.propagate``` and pass the model function
 which we want to propagate the uncertainty to. It is important that if
 the uncertainty quantification structure is defined for N-parameters (N=4
 in this case) the model function must accept all N parameters. Since we
 defined our model function to accept all N parameters already we do not
 need to worry about it.

1. Uncertainty of the dipolar signal fit: This case is easy, we already have the model and it is unconstrained


.. code-block:: python

    Vuq = paruq.propagate(Vmodel) # Uncertainty quantification for Vfit
    Vci95 = Vuq.ci(95) # 95#-confidence intervals for Vfit








2. Uncertainty of the distance distribution: In this case, the distribution has a non-negativity constraint which we
can specify via the lb input. 


.. code-block:: python

    lb = np.zeros_like(r) # Non-negativity constraint
    Puq = paruq.propagate(Pmodel,lb) # Uncertainty quantification for Pfit
    Pci95 = Puq.ci(95) # 95#-confidence intervals for Pfit








3. Uncertainty of the background: In this case, since we want to use this for plotting we need to evaluate
the function (1-lambda)*Bfit instead of just Bfit in order to plot the\
correct function.


.. code-block:: python

    Buq = paruq.propagate(lambda p:(1-p[0])*Bmodel(p)) # Uncertainty quantification for (1-lam)Bfit
    Bci95 = Buq.ci(95) # 95#-confidence intervals for (1-lam)Bfit








Plots
-----


.. code-block:: python


    plt.figure(figsize=(7,7))

    # Time-domain
    plt.subplot(211)
    plt.plot(t,V,'k.',t,Vfit,'r',t,(1-lamfit)*Bfit,'b',linewidth=1.5)
    plt.fill_between(t,Vci95[:,0],Vci95[:,1],color='r',alpha=0.3,linestyle='None')
    plt.fill_between(t,Bci95[:,0],Bci95[:,1],color='b',alpha=0.3,linestyle='None')
    plt.grid(alpha=0.3)
    plt.xlabel('t [$\mu s$]')
    plt.ylabel('V(t)')
    plt.legend(['data','Vfit','Bfit','Vfit 95%-CI','Bfit 95%-CI'])

    # Distance-domain
    plt.subplot(212)
    plt.plot(r,P,'k',r,Pfit,'r',linewidth=1.5)
    plt.fill_between(r,Pci95[:,0],Pci95[:,1],color='r',alpha=0.3,linestyle='None')
    plt.xlabel('r [nm]')
    plt.ylabel('P(r) [nm$^{-1}$]')
    plt.grid(alpha=0.3)
    plt.legend(['truth','Pfit','Pfit 95%-CI'])





.. image:: /auto_examples/images/sphx_glr_plot_param_uncertainty_propagation_001.png
    :alt: plot param uncertainty propagation
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    <matplotlib.legend.Legend object at 0x0000022605DAF518>




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  2.631 seconds)


.. _sphx_glr_download_auto_examples_plot_param_uncertainty_propagation.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_param_uncertainty_propagation.py <plot_param_uncertainty_propagation.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_param_uncertainty_propagation.ipynb <plot_param_uncertainty_propagation.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
