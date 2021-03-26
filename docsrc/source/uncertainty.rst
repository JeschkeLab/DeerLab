.. _uncertainty:

Uncertainty
=========================================

The presence of any level of noise in the data will result in uncertainty of any fitted quantities based on 
that data. This means that despite the fit returning single values for these quantities, there is a distributions
of values for each quantity which can describe the data. 

These so-called uncertainty distributions are important to determine the accuracy of the results. DeerLab provides a 
complete framework to calculate them for any kind of fitted quantity and to derive related quantities such as confidence
intervals. Confidence intervals (CI) provide a simple mean to report on the uncertainty, by defining a range of values within 
which the true solution might reside with a given probability.

DeerLab provides two means to estimate the uncertainty of fitted quantities: **covariance-matrix uncertainty**; a fast method, but not entirely accurate, and **bootstrapped uncertainty**; a much slower method, but more robust and accurate. 

The ``UQResult`` Object
---------------------------

The uncertainty estimation framework in DeerLab is contained into :ref:`UQResult` (Uncertainty Quantification) objects. 
These objects are variables returned by any function which calculates some kind of uncertainty estimate, and contain different 
fields with all quantities of interested related to uncertainty (see details :ref:`here <UQResult>`). Some of the most basic
attributes is ``UQResult.type``, which identifies whether the uncertainty was estimated by 
covariance-based or bootstrap methods. 

One ``UQResult`` can contain the uncertainty of multiple parameters and information on their correlations. For example, if ``fitmodel`` is used to fit a 
4-pulse DEER signal without background ::

    fit = dl.fitmodel(Vexp,t,r,'P',None,ex_4pdeer)  # Fit a 4-DEER form factor
    Puq = fit.Puncert           # Uncertainty quantification of fitted distance distribution
    lamuq = fit.exparamUncert   # Uncertainty quantification of fitted modulation depth

it will return a fit with a non-parametric distribution of N-points ``fit.P``, and the fitted modulation depth as ``fit.exparam``, then the output ``fit.Puncert`` will be a ``UQResult`` object with the information on all 
N-distribution elements and the output ``fit.exparamUncert`` will contain just the uncertainty of the modulation depth.

Confidence intervals
    As mentioned above, confidence intervals are the most practical quantities to report uncertainty of fit results. The ``UQResult.ci()`` is a method
    that takes the coverage or probability to be covered, and generates the confidence intervals. For the example above, if you want to generate the 95%-confidence 
    intervals you need to call ::

        P_ci = Puq.ci(95)       # 95%-confidence intervals of the distance distribution
        lam_ci = lamuq.ci(95)   # 95%-confidence intervals of the modulation depth parameter

    With this method you can calculate different confidence intervals for the same quantity, for example ::

        P_ci95 = Puq.ci(95) # 95%-confidence intervals of the distance distribution
        P_ci75 = Puq.ci(75) # 75%-confidence intervals of the distance distribution
        P_ci50 = Puq.ci(50) # 50%-confidence intervals of the distance distribution

    The confidence intervals are always returned as a ``Nx2``-array, where each one of the ``N`` parameters has two values, the lower and upper boundaries of the confidence interval. ::

        P_ci[:,0] # lower bound of the 95%-CI of the distance distribution
        P_ci[:,1] # upper bound of the 95%-CI of the distance distribution
        lam_ci[0] # lower bound of the 95%-CI of the modulation depth
        lam_ci[1] # upper bound of the 95%-CI of the modulation depth

Uncertainty distributions 
    A more complete description of the uncertainty are the uncertainty distributions for the fit parameter. These can be requested from the
    the ``UQResult.pardist`` method. Using ``UQResult.pardist(n)`` will return the parameter uncertainty probability density function 
    corresponding to the fit parameter with index ``n`` and its corresponding abscissa values. For example, ::

      P5_dist,P5_vals = Puq.pardit(5)       # Uncertainty distribution of fit.P[5] values
      lam_dist,lam_vals = lamuq.pardist(0)  # Uncertainty distribution of modulation depth values

Uncertainty Propagation
    Sometimes you might fit some parameters and get their corresponding uncertainty analysis, but you might be interested in knowing what the uncertainty
    is of a model that depends on those parameters. Analyzing the effect that the uncertainty of a set of parameters has on a dependent function is called 
    uncertainty or error propagation. 
    
    The ``UQResult.propagate`` method provides a simple interface for propagation uncertainty to dependent models or functions. The method just takes as input the function or model to which you 
    want to propagate the uncertainty. Additionally if the model has some boundaries (e.g. a distance distribution with non-negative values) you can specify the lower and upper bounds as additional inputs. 
    
    For example, if you fitted a Gaussian distribution with ``fitparamodel`` and obtained the uncertainty quantification for its parameters ::

        model  = lambda p: K@dl.dd_gauss(r,p) # Parametric model depending on p=[rmean,sigma] 
        fit = dl.fitparamodel(V,model) # Fit dipolar signal with parametric model 
        Pfit = dl.dd_gauss(r,fit.param) # Fitted distance distribution

        paramuq = fit.paramuq # Uncertainty quantification for fitted parameters

    you can propagate the uncertainty to the model ``dd_gauss`` as follows ::

        lb = np.zeros_like(Pfit)            # Non-negativity constraint of the distance distribution
        Puq = paruq.propagate(ddmodel,lb)   # Uncertainty quantification of fitted distribution

        Pci95 = Puq.ci(95) # 95%-confidence intervals of fitted distance distribution



Covariance Uncertainty Quantification
------------------------------------------

Covariance-baed uncertainty is the fastest method for uncertainty quantification available in DeerLab. Due to its readiness all fit functions 
in DeerLab return covariance-based ``UQResult`` objects for all quantities fitted.

This method estimates the uncertainty based on the curvature of the optimization surface. During optimization, when a minimum is found, the curvature
of the parameter hypersurface at that point is measured. A very sharp minimum means that there is less uncertainty in the found result, whereas a shallow
minimum will indicate a larger uncertainty in the results. 
This curvature is determined via the Jacobian of the optimization objective function and the corresponding covariance matrix. The method then assumes the uncertainty
distributions of all parameters to be normally distributed, to be centered at the fitted values, and their variance to be given by the diagonal elements of the covariance matrix.
In addition, the method does not take into consideration boundaries of the parameters, i.e. they are assumed to be unconstrained. However, in DeerLab confidence intervals and
uncertainty distributions are clipped at the boundaries as it is common practice. 

All these assumptions and approximation can lead to a less accurate estimate of the uncertainty. It is common for covariance-based confidence intervals to be overestimated and broader
than the bootstrapped equivalents.  However, their cheap computation cost makes them ideal for immediate estimations of uncertainty. 

Bootstrap Uncertainty Quantification
------------------------------------------

A more thorough way of assessing parameter uncertainty is bootstrap. In this method, many additional synthetic datasets (called bootstrap samples) are generated from the given experimental data 
and fitted in the same way as the original data. This means that the method samples the uncertainty arising when the same analysis is repeated for multiple noise realizations. 
When all  bootstrap samples have been analyzed, for each fitted quantity a distribution of values are obtained, which are taken as the uncertainty distributions for that quantity. 

Due to the need of repeating the fit for multiple bootstrap samples, this method can take long to estimate the uncertainty. However, since this method does not rely on 
any assumptions, the bootstrap uncertainty estimation are considered some of the most accurate, provided that enough bootstrap samples are taken. While a reduced number of
samples (50-100) can be used when testing workflows or new scripts, for conclusive analysis the minimum standard is to use at least about 1000 bootstrap samples. 

For routine analysis, the function ``fitmodel`` provides the option to switch from covariance based uncertainty quantification to bootstrapped uncertainty quantification by means of the keyword ``uq='bootstrap'``. The function will use 1000 bootstrap samples by default, however this can be changed if wanted by using ``uq=['bootstrap',Nsamples]`` when calling the function::

    fit = dl.fitmodel(Vexp,t,r,'P',dl.bg_hom3d,dl.ex_4pdeer,uq=['bootstrap',Nsamples],verbose=True)
    Puq = fit.Puncert # bootstrapped uncertainty quantification of fitted P(r)

However, in DeerLab you can calculate bootstrap uncertainty estimates of any quantities using the :ref:`bootan` function. The function takes the experimental data, the fit, and the analysis function. This analysis 
function must be a function that takes the experimental data and returns the quantities whose uncertainties are to be calculated. For example, to bootstrap the distance distribution, dipolar signal, and 
regularization parameter obtained from a 4-pulse DEER fit using ``fitregmodel`` you could use the following ::

    fit = dl.fitregmodel(Vexp,K,r,'tikh','aic')
    Vfit = K@fit.P # Fitted signal
    
    # Define the function to be bootstrapped
    def fitfcn(Vexp):
        fit = dl.fitregmodel(Vexp,K,r,'tikh','aic')
        V = K@fit.P
        return fit.P, V, fit.regparam  # bootstrap the fitted distance distribution, dipolar signal and regularization parameter

    bootuq = dl.bootan(fitfcn,Vexp,Vfit,samples=1000,verbose=True) # Bootstrap uncertainty quantification

The output of ``bootuq`` is an ``UQResult`` object equivalent to the ones obtained for covariance-based uncertainty analysis. If the fit procedure is slow or
costly, it is very recommendable to use the ``cores`` option to assign multiple CPU cores to the bootstrapping, in order to run different bootstrap samples in parallel, speeding
up the uncertainty estimation.
