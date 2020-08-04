class FitResult(dict):
    r""" Represents the results of a fit.
 
    Attributes
    ----------
    x : ndarray
        The solution of the optimization.
    success : bool
        Whether or not the optimizer exited successfully.
    cost : float
        Value of the cost function at the solution.
    residuals : ndarray
        Vector of residuals at the solution.
    stats : dict
        Goodness of fit statistical estimators:

        * ``stats['chi2red']`` - Reduced \chi^2 test
        * ``stats['r2']`` - R^2 test
        * ``stats['rmsd']`` - Root-mean squared deviation (RMSD)
        * ``stats['aic']`` - Akaike information criterion
        * ``stats['aicc']`` - Corrected Akaike information criterion
        * ``stats['bic']`` - Bayesian information criterion

    Notes
    -----
    There may be additional attributes not listed above depending of the
    specific fit function. Since this class is essentially a subclass of dict
    with attribute accessors, one can see which attributes are available
    using the `keys()` method. 
    """

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __repr__(self):
        if self.keys():
            m = max(map(len, list(self.keys()))) + 1
            return '\n'.join([k.rjust(m) + ': ' + repr(v)
                              for k, v in sorted(self.items())])
        else:
            return self.__class__.__name__ + "()"

    def __dir__(self):
        return list(self.keys())



























# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import numdifftools as nd
from scipy.stats import norm
from scipy.signal import fftconvolve
import copy

class UncertQuant:
    r""" Represens the uncertainty quantification of fit results.

    Attributes
    ----------
    type : string
        Uncertainty quanification approach:

            * 'covariance' - Covariance-based uncertainty analysis
            * 'bootstrap' - Bootstrapped uncertainty analysis
    
    mean : ndarray
        Mean values of the uncertainty distribution of the parameters.
    median : ndarray
        Median values of the uncertainty distribution of the parameters.
    std : ndarray
        Standard deviations of the uncertainty distribution of the parameters.
    covmat: ndarray
        Covariance matrix
    nparam: int scalar
        Number of parameters in the analysis.

    Methods
    -------

    """


    def __init__(self,uqtype,data,covmat=[],lb=[],ub=[]):

        #Parse inputs schemes
        if uqtype=='covariance':
            
            # Scheme 1: UncertQuant('covariance',parfit,covmat,lb,ub)
            parfit = data
            self.__parfit = parfit
            nParam = len(parfit)
            
        elif uqtype == 'bootstrap':
            # Scheme 2: UncertQuant('bootstrap',samples)
            samples = data
            self.__samples = samples
            nParam = np.shape(samples)[1]
                
        else:
            raise NameError('uqtype not found. Must be: ''covariance'' or ''bootstrap''.')

        if len(lb)==0:
            lb = np.full(nParam, -np.inf)
        
        if len(ub)==0:
            ub = np.full(nParam, np.inf)

        # Create confidence intervals structure
        if uqtype=='covariance':
            self.mean = parfit
            self.median = parfit
            self.std = np.sqrt(np.diag(covmat))
            self.covmat = covmat
                
        # Bootstrap-based CI specific fields
        elif uqtype == 'bootstrap':
            means = np.mean(samples,0)
            covmat = np.squeeze(samples).T@np.squeeze(samples)/np.shape(samples)[0] - means*means.T
            self.mean = means
            self.median = np.squeeze(np.median(samples,0))
            self.std = np.squeeze(np.std(samples,0))
            self.covmat = covmat
        self.type = uqtype

        # Set private variables
        self.__lb = lb
        self.__ub = ub
        self.nparam = nParam


    # Parameter distributions
    #--------------------------------------------------------------------------------
    def pardist(self,n):
        """
        Generate the uncertainty distribution of the n-th parameter

        Parameters
        ----------
        n : int scalar
            Index of the parameter 
        
        Returns
        -------
        ax : ndarray
            Parameter values at which the distribution is evaluated
        pdf : ndarray
            Probability density function of the parameter uncertainty.
        """
        if n > self.nparam or n < 0:
            raise ValueError('The input must be a valid integer number.')
        
        if self.type == 'covariance':
            # Generate Gaussian distribution based on covariance matrix
            sig = np.sqrt(self.covmat[n,n])
            xmean = self.mean[n]
            x = np.linspace(xmean-4*sig,xmean+4*sig,500)
            pdf = 1/sig/np.sqrt(2*np.pi)*np.exp(-((x-xmean)/sig)**2/2)
            
            # Clip the distributions at outside the boundaries
            pdf[x < self.__lb[n]] = 0
            pdf[x > self.__ub[n]] = 0

        if self.type == 'bootstrap':
            # Get bw using silverman's rule (1D only)
            samplen = self.__samples[:, n]
            sigma = np.std(samplen, ddof=1)
            bw = sigma * (len(samplen) * 3 / 4.0) ** (-1 / 5)

            # Make histogram
            maxbin = np.maximum(np.max(samplen),np.mean(samplen)+3*sigma)
            minbin = np.minimum(np.min(samplen),np.mean(samplen)-3*sigma)
            bins = np.linspace(minbin,maxbin, 2**10 + 1)
            count, edges = np.histogram(samplen, bins=bins)

            # Generate kernel
            delta = np.maximum(np.finfo(float).eps,(edges.max() - edges.min()) / (len(edges) - 1))
            kernel_x = np.arange(-4 * bw, 4 * bw + delta, delta)
            kernel = norm(0, bw).pdf(kernel_x)

            # Convolve
            pdf = fftconvolve(count, kernel, mode='same')

            # Set x coordinate of pdf to midpoint of bin
            x = edges[:-1] + delta

        # Ensure normalization of the probability density function
        pdf /= np.trapz(pdf, x)
        
        return x, pdf
    #--------------------------------------------------------------------------------


    # Parameter percentiles
    #--------------------------------------------------------------------------------
    def percentile(self,p):
        """
        Compute the p-th percentiles of the parameters uncertainty distributions

        Parameters
        ----------
        p : float scalar
            Percentile (between 0-100)
        
        Returns
        -------
        prctiles : ndarray
            Percentile values of all parameters
        """
        if p>100 or p<0:
            raise ValueError('The input must be a number between 0 and 100')
 
        x = np.zeros(self.nparam)
        for n in range(self.nparam):
            # Get parameter PDF
            values,pdf = self.pardist(n)
            # Compute corresponding CDF
            cdf = np.cumsum(pdf)
            cdf /= max(cdf)
            # Eliminate duplicates
            cdf, index = np.lib.arraysetops.unique(cdf,return_index=True)
            # Interpolate requested percentile
            x[n] = np.interp(p/100,cdf,values[index])

        return x
    #--------------------------------------------------------------------------------


    # Covariance-based confidence intervals
    #--------------------------------------------------------------------------------
    def ci(self,coverage):
        """
        Compute the confidence intervals for the parameters.

        Parameters
        ----------
        coverage : float scalar
            Coverage (confidence level) of the confidence intervals (between 0-100)
        
        Returns
        -------
        cis : 2D-ndarray
            Confidence intervals for the parameters:

                * ``cis[:,0]`` - Lower confidence intervals
                * ``cis[:,1]`` - Upper confidence intervals
        """
        if coverage>100 or coverage<0:
            raise ValueError('The input must be a number between 0 and 100')
        
        alpha = 1 - coverage/100
        p = 1 - alpha/2 # percentile
        
        x = np.zeros((self.nparam,2))
        if self.type=='covariance':
                # Compute covariance-based confidence intervals
                # Clip at specified box boundaries
                x[:,0] = np.maximum(self.__lb, self.__parfit - norm.ppf(p)*np.sqrt(np.diag(self.covmat)))
                x[:,1] = np.minimum(self.__ub, self.__parfit + norm.ppf(p)*np.sqrt(np.diag(self.covmat)))
                
        elif self.type=='bootstrap':
                # Compute bootstrap-based confidence intervals
                # Clip possible artifacts from the percentile estimation
                x[:,0] = np.minimum(self.percentile(p*100), np.amax(self.__samples))
                x[:,1] = np.maximum(self.percentile((1-p)*100), np.amin(self.__samples))

        return x



    # Error Propagation (covariance-based only)
    #--------------------------------------------------------------------------------
    def propagate(self,model,lbm=[],ubm=[]):
        """
        Uncertainty propagation. This function takes the uncertainty analysis of the 
        parameters and propagates it to another functon depending on those parameters.

        Parameters
        ----------
        model : callable
            Callable model function taking an array of ``nparam`` parameters.
        lbm : ndarray
            Lower bounds of the values returned by ``model``, by default assumed unconstrained.
        ubm : ndarray
            Upper bounds of the values returned by ``model``, by default assumed unconstrained.

        Returns
        -------
        modeluq : :ref:`UncertQuant`
            New uncertainty quantification analysis for the ouputs of ``model``.

        Notes
        -----
        Uncertainty propagation is covariance-based and so will be the resulting uncertainty analysis.
        """
        lbm,ubm = (np.atleast_1d(var) for var in [lbm,ubm])

        parfit = self.mean
        # Evaluate model with fit parameters
        modelfit = model(parfit)
        
        # Validate input boundaries
        if len(lbm)==0:
            lbm = np.full(len(modelfit), -np.inf)

        if len(ubm)==0:
            ubm = np.full(len(modelfit), np.inf)
        
        if len(modelfit)!=len(lbm) or len(modelfit)!=len(ubm):
            raise IndexError ('The 2nd and 3rd input arguments must have the same number of elements as the model output.')
        
        # Get jacobian of model to be propagated with respect to parameters
        J = np.reshape(nd.Jacobian(model)(parfit),(-1,parfit.size))
        
        # Clip at boundaries
        modelfit = np.maximum(modelfit,lbm)
        modelfit = np.minimum(modelfit,ubm)
        
        # Error progation
        modelcovmat = J@self.covmat@J.T
        
        # Construct new CI-structure for the model
        return  UncertQuant('covariance',modelfit,modelcovmat,lbm,ubm)
    #--------------------------------------------------------------------------------

