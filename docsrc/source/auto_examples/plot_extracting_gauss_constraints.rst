.. only:: html

    .. note::
        :class: sphx-glr-download-link-note

        Click :ref:`here <sphx_glr_download_auto_examples_plot_extracting_gauss_constraints.py>`     to download the full example code
    .. rst-class:: sphx-glr-example-title

    .. _sphx_glr_auto_examples_plot_extracting_gauss_constraints.py:


Extracting Gaussian constraints from a parameter-free distribution fit
=======================================================================

How to extract Gaussian constraints from a parameter-free fit.

While parameter-free distance distributions are the most robust way to
analyze dipolar signals, many structural biology modelling programs
accept only estimators such as mean distances or Gaussian constraints. 

This example shows how to extract Gaussian constraints from a
parameter-free fit of a dipolar signal and how to calculate the
corresponding uncertainty. 


.. code-block:: python


    import numpy as np
    import matplotlib.pyplot as plt
    from deerlab import *









.. code-block:: python


    # Generating a dataset
    # For this example we will simulate a simple 4pDEER signal

    # Parameters
    t = np.linspace(0,5,250)
    r = np.linspace(1,7,200)
    P = dd_gauss3(r,[4.5, 0.6, 0.4, 3, 0.4, 0.3, 4, 0.7, 0.5])
    lam = 0.3
    conc = 80; #uM

    # Simulate the signal
    Bmodel = lambda t,lam: bg_hom3d(t,conc,lam)
    K = dipolarkernel(t,r,lam,Bmodel)
    np.random.seed(0)
    V = K@P + whitegaussnoise(t,0.01)








Fit the dipolar signal
----------------------
First, we need to fit the parameter-free distance distribution using ``fitsignal()``
We are only interested right now on the fitted distribution and the
corresponding uncertainty quantification, so we will ignore the rest of
the outputs.
%%


.. code-block:: python

    fit = fitsignal(V,t,r,'P',bg_exp,ex_4pdeer,display=True)
    Pfit = fit.P



.. image:: /auto_examples/images/sphx_glr_plot_extracting_gauss_constraints_001.png
    :alt: plot extracting gauss constraints
    :class: sphx-glr-single-img


.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none

    ----------------------------------------------------------------------------
    Goodness of fit
      Vexp[0]: chi2 = 1.008065  RMSD  = 0.009641
    ----------------------------------------------------------------------------
    Fitted parameters and 95%-confidence intervals
      parfit['bg'][0][0]:   0.0768264  (0.0223348, 0.1313180)  Decay Rate (us-1)
      parfit['ex'][0][0]:   0.3129042  (0.2837829, 0.3420256)  Modulation depth ()
    ----------------------------------------------------------------------------




Extract Gaussian constraints from the fit
-----------------------------------------
Next, we will fit a multi-Gauss distribution to the fitted parameter-free
distribution. We can do this by using the ``fitparamodel()`` function (in
this example, fitting a two-Gauss model). 

However, in order to get the correct uncertainty quantification, we need
to specify the covariance matrix of the fitted distribution.
``fitparamodel()`` can then use that information to propagate the error in
``Pfit`` to the Gauss constraints that we then fit.

Extract the uncertainty quantification of the fitted distribution...


.. code-block:: python

    Pfit_uq = fit.Puncert
    # ...specifically its covariance matrix
    Pfit_covmat = Pfit_uq.covmat








Fit a 2-Gauss model to the fitted parameter-free distribution:

    - ``parfit```: will contain the Gaussian constraints
    - ``PGauss```: the corresponding distribution
    - ``paruq```: the uncertainty quantification of our constraints


.. code-block:: python

    Pmodel = lambda p: dd_gauss2(r,p)
    # Get information on the model
    info = dd_gauss2()
    par0 = info['Start']
    lb = info['Lower']
    ub = info['Upper']
    fit = fitparamodel(Pfit,Pmodel,par0,lb,ub,covmatrix=Pfit_covmat)
    parfit = fit.param
    paruq = fit.uncertainty
    PGauss = dd_gauss2(r,parfit)

    # Extract the 95#-confidence intervals...
    par95 = paruq.ci(95)
    # ... and print the results of the constraints 
    print('\nGaussian constraints:')
    info = dd_gauss2()
    for i in range(len(parfit)):
        print('  parfit[{}] = {:2.2f} ({:2.2f}, {:2.2f}) {}'.format(i,parfit[i],par95[i,0],par95[i,1],info['Parameters'][i]))

    # Now propagate the error of the constraints on the model
    lb = np.zeros_like(r) # Non-negativity constraint
    PGauss_uq = paruq.propagate(lambda par: dd_gauss2(r,par),lb)
    PGauss95 = PGauss_uq.ci(95)





.. rst-class:: sphx-glr-script-out

 Out:

 .. code-block:: none


    Gaussian constraints:
      parfit[0] = 3.02 (2.87, 3.17) Center of 1st Gaussian
      parfit[1] = 0.55 (0.29, 0.81) FWHM of 1st Gaussian
      parfit[2] = 0.25 (0.18, 0.33) Amplitude of 1st Gaussian
      parfit[3] = 4.29 (4.22, 4.36) Center of 2nd Gaussian
      parfit[4] = 0.94 (0.75, 1.13) FWHM of 2nd Gaussian
      parfit[5] = 0.67 (0.64, 0.70) Amplitude of 2nd Gaussian





.. code-block:: python


    # Plot the fitted constraints model on top of the parameter-free case
    plt.plot(r,Pfit,'r',linewidth=1.5)
    plt.fill_between(r,Pfit_uq.ci(95)[:,0], Pfit_uq.ci(95)[:,1],facecolor='r',linestyle='None',alpha=0.2)

    plt.plot(r,PGauss,'b',linewidth=1.5)
    plt.fill_between(r,PGauss95[:,0], PGauss95[:,1],facecolor='b',linestyle='None',alpha=0.2)

    plt.xlabel('Distance [nm]')
    plt.ylabel('P [nm$^{-1}$]')
    plt.tight_layout()
    plt.grid(alpha=0.3)
    plt.legend(['Fit','95%-CI','2G-constraints','95%-CI'])
    plt.show()





.. image:: /auto_examples/images/sphx_glr_plot_extracting_gauss_constraints_002.png
    :alt: plot extracting gauss constraints
    :class: sphx-glr-single-img






.. rst-class:: sphx-glr-timing

   **Total running time of the script:** ( 0 minutes  5.704 seconds)


.. _sphx_glr_download_auto_examples_plot_extracting_gauss_constraints.py:


.. only :: html

 .. container:: sphx-glr-footer
    :class: sphx-glr-footer-example



  .. container:: sphx-glr-download sphx-glr-download-python

     :download:`Download Python source code: plot_extracting_gauss_constraints.py <plot_extracting_gauss_constraints.py>`



  .. container:: sphx-glr-download sphx-glr-download-jupyter

     :download:`Download Jupyter notebook: plot_extracting_gauss_constraints.ipynb <plot_extracting_gauss_constraints.ipynb>`


.. only:: html

 .. rst-class:: sphx-glr-signature

    `Gallery generated by Sphinx-Gallery <https://sphinx-gallery.github.io>`_
