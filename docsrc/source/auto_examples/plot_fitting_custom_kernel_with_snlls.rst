.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_plot_fitting_custom_kernel_with_snlls.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_auto_examples_plot_fitting_custom_kernel_with_snlls.py:


Fitting a custom kernel model with a parameter-free distribution
=================================================================

How the use of SNLLS to fit a kernel model and a parameter-free 
distribution to a dipolar signal.


.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from deerlab import *








Generating a dataset
-----------------------------------------------------------------------------
 For this example we will simulate a simple 4pDEER signal


.. code-block:: python


    t = np.linspace(-0.5,5,300)
    r = np.linspace(2,6,200)

    # Generate ground truth and input signal
    P = dd_gauss2(r,[3.5, 0.4, 0.4, 4.5, 0.7, 0.6])
    lam = 0.36
    c0 = 250 #uM
    B = bg_hom3d(t,c0,lam)
    K = dipolarkernel(t,r,lam,B)
    V = K@P  + whitegaussnoise(t,0.01)








Fitting via SNLLS
------------------
 Now in order to fit a non-linear dipolar kernel model ``Kmodel`` and a
 linear parameter-free distance distribution ``Pfit`` simultaneously, we
 can use the separable non-linear least squares ``SNLLS`` method. 

 First we define the function that contains the model for the dipolar kernel we want to fit. It 
 is a non-linear functon that accepts the parameter array ``p`` and returns the 
 fitted dipolar kernel ``K``. The linear parameters, in this case ``P``, are
 computed by solving a Tikhonov-regularized linear LSQ problem automatically in the ``snlls`` function. 


.. code-block:: python


    def Kmodel(p):

        # Unpack parameters
        lam,c0 = p
        # Get background
        B = bg_hom3d(t,c0,lam)
        # Generate 4pDEER kernel
        K = dipolarkernel(t,r,lam,B)

        return K








Next, there are two different parameter sets being fitted at the same time:
linear and non-linear parameters. Therefore, the lower/upper bounds for
the two sets need (or can) be specified.


.. code-block:: python


    #--------------------------
    # Non-linear parameters:
    #--------------------------
    #       lam  c0
    #--------------------------
    par0 = [0.5, 50 ] # Start values
    lb   = [ 0, 0.05] # lower bounds
    ub   = [ 1, 1000] # upper bounds

    #--------------------------
    # Linear parameters: 
    #--------------------------
    #          Pfit
    #--------------------------
    lbl = np.zeros_like(r) # Non-negativity constraint of P
    ubl = [] # Unconstrained upper boundary

    # Run SNLLS optimization
    fit = snlls(V,Kmodel,par0,lb,ub,lbl,ubl)
    parfit = fit.nonlin
    Pfit = fit.lin
    paruq = fit.uncertainty

    # Get non-linear parameters uncertainty
    param95 = paruq.ci(95,'nonlin')  #  95#-confidence interval

    # Get linear parameters (distribution) uncertainty
    Pci50 = paruq.ci(50,'lin') #  50#-confidence interval
    Pci95 = paruq.ci(95,'lin') #  95#-confidence interval

    # Print result
    print('lambda = {:.2f}({:.2f}-{:.2f})'.format(parfit[0],param95[0,0],param95[0,1]))
    print('c0 = {:.2f}({:.2f}-{:.2f})uM'.format(parfit[1],param95[1,0],param95[1,1]))

    # Get fitted model
    Kfit = Kmodel(parfit)
    Vfit = Kfit@Pfit





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    lambda = 0.36(0.36-0.37)
    c0 = 254.67(243.34-266.00)uM




Plots
------


.. code-block:: python


    plt.subplot(211)
    plt.plot(t,V,'k.',t,Vfit,'b')
    plt.grid(alpha=0.3)
    plt.xlabel('t [$\mu s$]')
    plt.ylabel('V(t)')
    plt.legend(['data','fit'])

    plt.subplot(212)
    plt.plot(r,P,'k',r,Pfit,'b')
    plt.fill_between(r,Pci50[:,0],Pci50[:,1],color='b',alpha=0.4,linestyle='None')
    plt.fill_between(r,Pci95[:,0],Pci95[:,1],color='b',alpha=0.2,linestyle='None')
    plt.grid(alpha=0.3)
    plt.xlabel('r [nm]')
    plt.ylabel('P(r) [nm$^{-1}$]')
    plt.legend(['truth','fit','50%-CI','95%-CI'])








.. image:: /auto_examples/images/sphx_glr_plot_fitting_custom_kernel_with_snlls_001.png
    :alt: plot fitting custom kernel with snlls
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    <matplotlib.legend.Legend object at 0x000001ED98BB1518>




.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  2.488 seconds)


.. _sphx_glr_download_auto_examples_plot_fitting_custom_kernel_with_snlls.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_fitting_custom_kernel_with_snlls.py <plot_fitting_custom_kernel_with_snlls.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_fitting_custom_kernel_with_snlls.ipynb <plot_fitting_custom_kernel_with_snlls.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
