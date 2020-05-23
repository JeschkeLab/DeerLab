Uncertainty
=========================================

After fitting experimental data with a distance distribution, a background, and a modulation depth, it is important to assess the uncertainty (i.e. the error) in the fitted parameters.

DeerLab provides uncertainty estimates for all parameters in the form of confidence intervals (CIs). They can be calculated in two ways: either using the covariance matrix, or using bootstrap. The first method is fast, but is not entirely accurate, whereas bootstrap is much slower, but more robust.

Covariance CIs
------------------------------------------

Along with every fit, functions like ``fitsignal`` return confidence intervals for the fitted parameters and the fitted distributions. For example,

.. code-block:: matlab

   [Vfit,Pfit,Bfit,parfit,parCI] = fitsignal(Vexp,t,r);

The output ``parCI`` contains the 95% confidence interval for all fitted parameters, calculated using the standard method based on the variance-covariance matrix.

Bootstrap CIs
------------------------------------------

A more thorough way of assessing parameter uncertainty is bootstrap. In this method, many additional synthetic datasets are generated from the given experimental data and the fitted model and are fitted individually. This yields an ensemble of parameter fits that is analyzed statistically to provide information about the scatter.

Here is an example for a parametric model:

.. code-block:: matlab

    bootci = bootan(@(V)fitfcn(V,r,K),Vexp,Vfit,1000,'verbose',true);
    
    function parfit = fitfcn(Vin,r,K)
           parfit = fitparamodel(Vin,@dd_gauss,r,K);
    end

The output ``bootci`` contains calculated bootstrap 50%, 90%, and 95% confidence intervals for all parameters.

Here is an example for a model with a non-parametric distribution:

.. code-block:: matlab

    bootci = bootan(@(V)fitfcn(V,t,r),Vexp,Vfit,100,'verbose',true);

    function [Pfit, parfit.bg, parfit.ex] = fitfcn(Vin,t,r)
           [~,Pfit,~,parfit] = fitsignal(Vin,t,r,'P',@bg_hom3d,@ex_4pdeer,[],'RegParam',1);
    end

To plot the resulting CIs for the non-parametric distance distribution, use

.. code-block:: matlab
    
    Pci = bootci{1}.ci95;
    plot(r,Pfit,'k',r,Pci(:,1),'r',r,Pci(:,2),'r')
