Uncertainty
=========================================

After fitting experimental data with a distance distribution, a background, and a modulation depth, it is important to assess the uncertainty (i.e. the error) in the fitted parameters.

DeerLab provides uncertainty estimates for all parameters in the form of confidence intervals (CIs). They can be calculated in two ways: either using the covariance matrix, or using bootstrap. The first method is fast, but is not entirely accurate, whereas bootstrap is much slower, but more robust. All confidence intervals provided by DeerLab functions describe the range of values which might contain the ground truth with a ceratin probability.

Covariance Uncertainty Quantification
------------------------------------------

Along with every fit, functions like ``fitsignal`` return confidence intervals for the fitted parameters and the fitted distributions. For example,

.. code-block:: matlab

   [Vfit,Pfit,Bfit,parfit,paruq] = fitsignal(Vexp,t,r);

The output ``paruq`` is a uncertainty quantification structure which contains the uncertainty information for all fitted parameters, calculated using the standard method based on the variance-covariance matrix.

The confidence intervals (at any confidence level) can be calculated by using the ``paruq.ci()`` field, e.g. the 50%, 75% and 95% confidence intervals: 

.. code-block:: matlab

    parfit_ci50 = paruq.ci(50);
    parfit_ci75 = paruq.ci(75);
    parfit_ci95 = paruq.ci(95);

Uncertainty can be propagated to dependent models, e.g. propagating the uncertainty in the fit of ``rmean`` and ``fwhm`` to the resulting Gaussian distance distribution. This can be done via the ``paruq.propagate()`` function: 

.. code-block:: matlab

    % parfit = [rmean fwhm]
    ddmodel = @(parfit)dd_gauss(r,parfit)
    
    %Get the fitted model
    Pfit = ddmodel(parfit);
    
    %Propagate the error in the parameters to the model
    lb = zeros(numel(Pfit),1);
    ub = [];
    Pci = paruq.propagate(ddmodel,lb,ub)

    % Get the 95%-confidence intervals of the fitted distribution
    Pci95 = Pci.ci(95);


Assumptions:
   - The uncertainty in the fitted parameters is described by a Gaussian distribution.
   - The mean of this Gaussian is assumed to be the fitted value and its width by the diagonal elements of the covariance matrix.
   - All parameters are assumed to be unconstrained.


Bootstrap Uncertainty Quantification
------------------------------------------

A more thorough way of assessing parameter uncertainty is bootstrap. In this method, many additional synthetic datasets are generated from the given experimental data and the fitted model and are fitted individually. This yields an ensemble of parameter fits that is analyzed statistically to provide information about the scatter.

Here is an example for a parametric model:

.. code-block:: matlab

    bootuq = bootan(@(V)fitfcn(V,r,K),Vexp,Vfit,1000,'verbose',true);
    
    function parfit = fitfcn(Vin,r,K)
           parfit = fitparamodel(Vin,@dd_gauss,r,K);
    end

The output ``bootuq`` structure can be analogously used to evaluate confidence intervals at different confidence levels, e.g the 50%, 75% and 95% confidence intervals: 

.. code-block:: matlab

    parfit_ci50 = bootuq.ci(50);
    parfit_ci75 = bootuq.ci(75);
    parfit_ci95 = bootuq.ci(95);

The bootstrapped distributions for each parameter can be accessed by using the ``paruq.pardist()`` field, e.g.if the modulation depth is the second fit parameter:

.. code-block:: matlab

    moddepth_dist = bootuq.pardist(2);


Here is an example for a model with a non-parametric distribution:

.. code-block:: matlab

    bootuq = bootan(@(V)fitfcn(V,t,r),Vexp,Vfit,100,'verbose',true);

    function [Pfit, parfit.bg, parfit.ex] = fitfcn(Vin,t,r)
           [~,Pfit,~,parfit] = fitsignal(Vin,t,r,'P',@bg_hom3d,@ex_4pdeer,[],'RegParam',1);
    end

To plot the resulting 95% and 50% confidence interval for the non-parametric distance distribution, use

.. code-block:: matlab
    
    Pci50 = bootuq.ci(50);
    Pci95 = bootuq.ci(95);
    
    plot(r,Pfit,'k')
    fill([r fliplr(r)],[Pci50(:,1); flipud(Pci50(:,2))],'r','FaceAlpha',0.5)
    fill([r fliplr(r)],[Pci95(:,1); flipud(Pci95(:,2))],'r','FaceAlpha',0.2)

Assumptions:
   - ``Vfit`` is a good fit of the experimental data ``Vexp``.

.. _cireference:

Uncertainty Quantification Reference
------------------------------------------
All DeerLab functions which return any kind of uncertainty estimation (covariance-baed or bootstrapped) will return a so-called uncertainty quantification structure. When fitting *N* parameters or e.g. an *N*-element distance distribution, it has the following structure.

``uqstruct`` - Uncertainty quantification structure containing the following fields:

------------------------------------------

    **Confidence intervals**
    
    *   ``.ci(c)`` - Function handle that returns the confidence interval of the fitted parameters for a given coverage or confidence level ``c``

            Inputs:
            
                *   ``c`` - Coverage/Confidence level (scalar, in range [0,100])
            Returns:
            
                *   ``parCI`` - confidence intervals of the *N*-parameters (*Nx2*-matrix, ``parCI(:,1)`` - upper bound, ``parCI(:,2)`` - lower bound)


---------------------------

    **Parameter distributions**

    *   ``.pardist(n)`` - Function handle that returns the distribution of the *n*-th fitted parameter

            Inputs:
            
                *   ``n`` - Index of the fitted parameter (scalar, integer in range [1,N])
            Returns:
            
                *   ``dist`` - Distribution of the n-th fitted parameter (struct)

                        * ``.values`` - evaluated parameter values
                        * ``.pdf`` - probability densities of the parameter values

    *   ``.mean`` - Means of the parameter distributions (*N*-element array)
    *   ``.median`` - Medians of the parameter distributions (*N*-element array)
    *   ``.std`` - Standard deviation of the parameter distributions (*N*-element array)
    *   ``.percentile(p)`` - Function handle that returns the *p*-th percentiles of the parameter distribution

            Inputs:
            
                *   ``p`` - Percentile (scalar, in range [0,100])
            Returns:
            
                *   ``perct`` - Percentiles of the parameter distributions (*N*-element array)

---------------------------

    **Covariance & Propagation**

    *   ``.covmat`` - Covariance matrix for the fit parameters (*NxN* matrix)
    *   ``.propagate(model,lb,ub)`` - Function handle that propagates the uncertainty unto another model (based on quadratic-approximation)


            Inputs:
            
                *   ``@model`` - Function handle of the model/function to propagate the error unto (must accept all *N*-parameters as input)
                *   ``lb`` - Lower bounds of the results returned by ``model``, if empty assumed to be unbounded.
                *   ``ub`` - Upper bounds of the results returned by ``model``, if empty assumed to be unbounded.

            Returns:
            
                *   ``model_cist`` - Confidence interval structure for the results of ``model``.