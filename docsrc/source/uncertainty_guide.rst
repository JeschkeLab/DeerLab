.. _uncertainty:

Uncertainty Guide 
=========================================

The mere presence of noise in a dataset results in an uncertain estimation of model parameters. This uncertainty results in a distribution of model parameter values that could describe the dataset. 

Basics
------

Here is a short selection of terms used in DeerLab related to uncertainty: 

Uncertainty distributions 
    The distribution of values for a given variable that quantifies the uncertainty in said variable's estimate. It arises from the noise/errors in the datasets used to estimate that quantity. A very sharp distribution means less uncertainty in the found result, whereas a broad distribution will indicate more significant uncertainty. 
    Many statistical descriptors of the uncertainty distribution, such as the mean, standard deviation, and percentiles, can describe the uncertainty. 
Confidence intervals
    Confidence intervals (CI) provide a simple means to report the uncertainty by defining a range of values within which the ground truth might reside within a given probability. The 95%-confidence intervals are a golden standard for uncertainty quantification and should always be reported along with any fitted quantities.  
Uncertainty propagation 
    Refers to the analysis of the effects of quantity's uncertainty on a dependent function or quantity. It is also commonly referred to as error propagation. Allows the quantification of uncertainty in quantities that might not be directly related to the model fitted to the dataset. 

Uncertainty quantification methods 
----------------------------------

DeerLab provides three different methods for the estimation of uncertainty: 
 - Covariance-based or asymptotic method (least accurate, fastest)
 - Profile-likelihood method   (accurate, slow)
 - Bootstrapping method (most accurate, slowest)

The following section will provide the basic concepts related to these methods: 

Covariance-based or asymptotic method 
*************************************

Covariance-based uncertainty is the fastest method for uncertainty quantification available in DeerLab. Due to its readiness, the ``fit`` function in DeerLab will employ it by default unless another method is requested.

This method estimates the uncertainty based on the curvature of the optimization surface. It assumes the uncertainty distributions of all parameters to be normally distributed, centered at the fitted values (i.e., the maximum-likelihood estimates), and their variance to be given by the curvature of the objective function. This curvature is determined via the Jacobian of the objective optimization function and the corresponding covariance matrix. 

In addition, the method does not consider the boundaries of the parameters, i.e., they are assumed to be unconstrained. However, in DeerLab, confidence intervals and uncertainty distributions are clipped at the boundaries, as is standard practice. 

All these assumptions and approximations can lead to a less accurate estimate of the uncertainty. It is common for covariance-based confidence intervals to be overestimated and broader than those returned by the other methods offered by DeerLab. Due to the assumption of a normal uncertainty distribution, covariance-based confidence intervals will always be symmetrical about the fitted values, which might not always be accurate. However, their cheap computation cost makes them ideal for immediate and preliminary estimations of uncertainty. 


Profile-likelihood method
*************************************

Likelihood-based uncertainty estimated from the profile-likelihood method represents a significant improvement in accuracy with respect to the covariance-based estimates. This method is, however, limited to non-linear model parameters. 

This method estimates the uncertainty based on the profile-likelihood of a non-linear model parameter. The profile of a model parameter is computed by fitting the model to the data repeatedly while keeping that parameter fixed at different values. This results in an isometric profile of the likelihood or objective function. A likelihood uncertainty distribution or confidence intervals can be derived from this profile using a statistical criterion.

The number of values at which the model is repeatedly fitted determines the accuracy of the profile estimation and, therefore, of the uncertainty estimation. For setting up scripts and quick tests, 20 repeats provide rough estimates, but at least 50 repeats are recommended for final results.   

This method does take model parameter boundaries into account, and most importantly, does not assume any shape of the underlying uncertainty distribution. Therefore, it provides better estimates of the uncertainty and allows for asymmetric confidence intervals. However, the repeated optimization nature of the method results in moderately more significant computation times for the uncertainty estimation. It is the restriction to non-linear parameters that makes it also less generally applicable. 

This method is also critical for the identifiability analysis of model parameters.  

Bootstrap method
*************************************

A gold standard to accurately estimate parameter uncertainty is bootstrapping. In this method, many additional synthetic datasets (called bootstrap samples) are generated from the given dataset(s), and the model is fitted against them in the same fashion as the original dataset(s). This means that the method samples the uncertainty when the same analysis is repeated for multiple noise realizations. 
When all bootstrap samples have been analyzed, a distribution of values is obtained for each fitted quantity, which is taken as the uncertainty distribution for that quantity. 

Due to the need to repeat the fit for multiple bootstrap samples, this method can take a long time to estimate the uncertainty. However, since this method does not rely on any assumptions, the bootstrap uncertainty estimation is considered some of the most accurate, provided that enough bootstrap samples are taken. While a reduced number of samples (50-100) can be used when testing workflows or new scripts, the minimum standard is to use at least about 1000 bootstrap samples for conclusive analysis. 

Selecting the uncertainty method 
--------------------------------

The choice of uncertainty quantification method can be easily controlled via optional keyword arguments in the ``fit`` function. The program will then internally employ the appropriate methods to report the requested uncertainty estimates.

Covariance-based method
************************
The covariance-based method is the default method. Therefore, the ``fit`` function can be passed without any additional inputs: ::

    # Fit the model to the data and quantify the uncertainty based on covariance
    fitresult = dl.fit(model, y)

When using bootstrapping, the parameter fits ``fit.<parameter>`` will not correspond to the maximum-likelihood but to the median of the bootstrapped uncertainty quantification (a better and more accurate estimate).   

Bootstrap method
*************************************

The bootstrap method can be requested upon fitting by passing the keyword argument ``bootstrap`` to the ``fit`` function along with the number of bootstrap samples to be taken: ::

    # Fit the model to the data and quantify the uncertainty via bootstrapping
    fitresult = dl.fit(model, y, bootstrap=1000)

When using bootstrapping, the parameter fits ``fit.<parameter>`` will not correspond to the maximum-likelihood but to the median of the bootstrapped uncertainty quantification (a better and more accurate estimate).   

The ``UQResult`` Object
---------------------------

The results uncertainty estimation in DeerLab is contained into :ref:`UQResult` (Uncertainty Quantification result) objects. 
These objects contain all the quantities of interest related to the uncertainty of one or several quantities. 


Confidence intervals
    As mentioned above, confidence intervals are the most practical quantities to report the uncertainty of fit results. They can be computed for arbitrary confidence levels using the  ``ci`` method of the ``UQResult`` object. This method takes the coverage probability (or confidence level) and generates the confidence intervals. For example, to get the 95% confidence intervals of a fitted parameter ::

        # Get the 95% confidence intervals
        ci95 = fitresult.<parameter>Uncert.ci(95)
        # The confidence interval is a list containing the lower/upper bounds
        ci_lower, ci_upper = ci95

    With this method you can calculate different confidence intervals for the same quantity, for example ::

        ci95 = fitresult.<parameter>Uncert.ci(95) # 95%-confidence intervals of the parameter
        ci75 = fitresult.<parameter>Uncert.ci(75) # 75%-confidence intervals of the parameter
        ci50 = fitresult.<parameter>Uncert.ci(50) # 50%-confidence intervals of the parameter

    For vector quantities, confidence intervals are always returned as a ``Nx2``-array, where each of the ``N`` elements of the vector has two values, the lower and upper boundaries of the confidence interval. ::

        # Get the confidence intervals on the model response vector
        response_ci = fitresult.modelUncert.ci(95)

        response_ci[:,0] # lower bound of the 95%-CI of the distance distribution
        response_ci[:,1] # upper bound of the 95%-CI of the distance distribution


Uncertainty distributions 
    A complete description of the uncertainty is the uncertainty distributions for the fit parameter. These can be requested from the ``pardist`` method. Using ``pardist(n)`` will return the uncertainty probability density function and its abscissa values for the corresponding quantity's ``n``-th element. For example, ::

        pardist = fitresult.<parameter>Uncert.pardist(0) # Get the parameter uncertainty distribution
        modeldist5 = fitresult.modelUncert.pardist(4) # Get the uncertainty distribution of the model's response 5th element
