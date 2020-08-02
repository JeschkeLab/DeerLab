Uncertainty
=========================================

After fitting experimental data with a distance distribution, a background, and a modulation depth, it is important to assess the uncertainty in the fitted parameters.

DeerLab provides uncertainty estimates for all parameters in the form of confidence intervals (CIs). They can be calculated in two ways: either using the covariance matrix, or using bootstrap. The first method is fast, but is not entirely accurate, whereas bootstrap is much slower, but more robust. All confidence intervals provided by DeerLab functions describe the range of values which might contain the ground truth with a ceratin probability.

Covariance Uncertainty Quantification
------------------------------------------

Along with every fit, functions like ``fitsignal`` return confidence intervals for the fitted parameters and the fitted distributions. For example,

.. code-block:: python

   fit = fitsignal(Vexp,t,r)
   Puq = fit.Puncert            # Uncertainty quantification of the distance distribution
   Buq = fit.Buncert            # Uncertainty quantification of the background 
   lamuq = fit.exparamUncert    # Uncertainty quantification of the modulation depth
   

The variables ``Puq``, ``Buq`` and ``lamuq`` are uncertainty quantification objects :ref:`UncertQuant` which contain the full uncertainty information of the corresponding variables, calculated using the standard method based on the variance-covariance matrices.

The confidence intervals (at any confidence level) can be calculated by using the ``ci()`` method, e.g. for the 50% and 95% confidence intervals of the fitted distance distribution are calculated as: 

.. code-block:: python

    Pfit_ci50 = paruq.ci(50)    # 50%-confidence intervals of Pfit
    Pfit_ci95 = paruq.ci(95)    # 95%-confidence intervals of Pfit

Uncertainty can also be propagated to dependent models. For example, assume that we have fitted a single Gaussian distance distribution with ``rmean`` and ``fwhm`` as parameters. Now, we can propagate the uncertainty in the fit of ``rmean`` and ``fwhm`` to the resulting Gaussian distance distribution. This can be done via the ``propagate()`` method, this will create the uncertainty quantification for the fitted distribution: 

.. code-block:: python

    # parfit = [rmean, fwhm]
    ddmodel = lambda parfit: dd_gauss(r,parfit)
    
    # Get the fitted model
    Pfit = ddmodel(parfit)
    
    # Propagate the error in the parameters to the model
    lb = zeros(numel(Pfit),1)   # Non-negativity constraint of the distance distribution
    ub = []]                    # No upper bounds
    Puq = paruq.propagate(ddmodel,lb) 

    # Get the 95%-confidence intervals of the fitted distribution
    Pci95 = Puq.ci(95)


Theoretical assumptions of covariance-based uncertainty analysis:
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
        fit = fitparamodel(V,Vmodel,r,K)
        return fit.param

    bootuq = bootan(fitfcn,Vexp,Vfit,samples=1000,verbose=True);

The output ``bootuq`` is again a :ref:`UncertQuant` object that can be used as described above to evaluate confidence intervals at different confidence levels, e.g the 50% and 95% confidence intervals: 

.. code-block:: python

    parfit_ci50 = bootuq.ci(50)
    parfit_ci95 = bootuq.ci(95)

The bootstrapped distributions for each parameter can be accessed by using the ``pardist()`` method, e.g.if the modulation depth is the second fit parameter:

.. code-block:: python

    moddepth_dist = bootuq.pardist(2);


Here is an example for a model with a non-parametric distribution:

.. code-block:: python


    def fitfcn(V):
           fit = fitsignal(V,t,r,'P',bg_hom3d,ex_4pdeer)
        return fit.P, fit.bgparam, fit.exparam

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