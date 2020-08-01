Uncertainty
=========================================

After fitting experimental data with a distance distribution, a background, and a modulation depth, it is important to assess the uncertainty in the fitted parameters.

DeerLab provides uncertainty estimates for all parameters in the form of confidence intervals (CIs). They can be calculated in two ways: either using the covariance matrix, or using bootstrap. The first method is fast, but is not entirely accurate, whereas bootstrap is much slower, but more robust. All confidence intervals provided by DeerLab functions describe the range of values which might contain the ground truth with a ceratin probability.

Covariance Uncertainty Quantification
------------------------------------------

Along with every fit, functions like ``fitsignal`` return confidence intervals for the fitted parameters and the fitted distributions. For example,

.. code-block:: python

   Vfit,Pfit,Bfit,parfit,paruq,_,_ = fitsignal(Vexp,t,r);

The variable ``paruq`` is a uncertainty quantification object :ref:`UncertQuant` which contains the uncertainty information for all fitted parameters, calculated using the standard method based on the variance-covariance matrix.

The confidence intervals (at any confidence level) can be calculated by using the ``paruq.ci()`` method, e.g. the 50%, 75% and 95% confidence intervals: 

.. code-block:: python

    parfit_ci50 = paruq.ci(50)
    parfit_ci75 = paruq.ci(75)
    parfit_ci95 = paruq.ci(95)

Uncertainty can be propagated to dependent models, e.g. propagating the uncertainty in the fit of ``rmean`` and ``fwhm`` to the resulting Gaussian distance distribution. This can be done via the ``paruq.propagate()`` method: 

.. code-block:: python

    # parfit = [rmean fwhm]
    ddmodel = lambda parfit: dd_gauss(r,parfit)
    
    # Get the fitted model
    Pfit = ddmodel(parfit)
    
    # Propagate the error in the parameters to the model
    lb = zeros(numel(Pfit),1)
    ub = []
    Pci = paruq.propagate(ddmodel,lb,ub)

    # Get the 95%-confidence intervals of the fitted distribution
    Pci95 = Pci.ci(95)


Assumptions:
   - The uncertainty in the fitted parameters is described by a Gaussian distribution.
   - The mean of this Gaussian is assumed to be the fitted value and its width by the diagonal elements of the covariance matrix.
   - All parameters are assumed to be unconstrained.


Bootstrap Uncertainty Quantification
------------------------------------------

A more thorough way of assessing parameter uncertainty is bootstrap. In this method, many additional synthetic datasets are generated from the given experimental data and the fitted model and are fitted individually. This yields an ensemble of parameter fits that is analyzed statistically to provide information about the scatter.

Here is an example for a parametric model:

.. code-block:: python

    def fitfcn(V):
        Vmodel = lambda par: K@dd_gauss(r,par)
        parfit,_,_ = fitparamodel(V,Vmodel,r,K)
        return parfit

    bootuq = bootan(fitfcn,Vexp,Vfit,samples=1000,verbose=True);

The output ``bootuq`` is again a :ref:`UncertQuant` object that can be used as described above to evaluate confidence intervals at different confidence levels, e.g the 50%, 75% and 95% confidence intervals: 

.. code-block:: python

    parfit_ci50 = bootuq.ci(50)
    parfit_ci75 = bootuq.ci(75)
    parfit_ci95 = bootuq.ci(95)

The bootstrapped distributions for each parameter can be accessed by using the ``paruq.pardist()`` method, e.g.if the modulation depth is the second fit parameter:

.. code-block:: python

    moddepth_dist = bootuq.pardist(2);


Here is an example for a model with a non-parametric distribution:

.. code-block:: python


    def fitfcn(V):
           _,Pfit,_,parfit,_,_,_ = fitsignal(V,t,r,'P',bg_hom3d,ex_4pdeer)
        return Pfit, parfit.bg, parfit.ex

    bootuq = bootan(fitfcn,Vexp,Vfit,samples=100,verbose=True)

To plot the resulting 95% and 50% confidence interval for the non-parametric distance distribution, use

.. code-block:: python
    
    Pci50 = bootuq.ci(50)
    Pci95 = bootuq.ci(95)
    
    import matplotlib.pyplot as plt
    plt.plot(r,Pfit,'k')
    plt.fill_between(r,Pci50[:,0]; Pci50[:,1],color='r',alpha=,0.5)
    plt.fill_between(r,Pci95[:,0]; Pci95[:,1],color='r',alpha=0.2)

Assumptions:
   - ``Vfit`` is a good fit of the experimental data ``Vexp``.