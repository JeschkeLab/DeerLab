# uqst.py - Uncertainty quantification structure constructor
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.
import numpy as np
from deerlab.utils import jacobianest
from scipy.stats import norm
from scipy.signal import fftconvolve
import copy

class uqst:
    """
    UQST Uncertainty quantification structure constructor

    uqstruct = UQST('covariance',parfit,covmat,lb,ub)
    Constructs the structure for covariance-based uncertainty quantificaton.

    uqstruct = UQST('bootstrap',samples)
    Constructs the structure for bootstrapped uncertainty quantificaton.

    Inputs:
        parfit     N-element array of fitted values used as mean values in covariance-based CI
        covmat     NxN-element covariance matrix
        lb         N-element array of lower bounds used for fitting parfit
        ub         N-element array of upper bounds used for fitting parfit
        samples    MxN-element matrix of bootstrapped samples results
    """


    def __init__(self,uqtype,data,covmat,lb,ub):

        #Parse inputs schemes
        if uqtype=='covariance':
            
            # Scheme 1: uqst('covariance',parfit,covmat,lb,ub)
            parfit = data
            self.__parfit = parfit
            nParam = len(parfit)
            if len(lb)==0:
                lb = np.full(nParam, -np.inf)
            
            if len(ub)==0:
                ub = np.full(nParam, np.inf)
            
        elif uqtype == 'bootstrap':
            # Scheme 2: uqst('bootstrap',samples)
            samples = data
            self.__samples = samples
            nParam = np.shape(samples)[1]
                
        else:
            raise NameError('uqtype not found. Must be: ''covariance'' or ''bootstrap''.')

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
        self.uqtype = uqtype

        # Set private variables
        self.__uqtype = uqtype
        self.__covmat = covmat
        self.__lb = lb
        self.__ub = ub
        self.__nParam = nParam


    # Parameter distributions
    #--------------------------------------------------------------------------------
    def pardist(self,n):
        
        if n > self.__nParam or n < 0:
            raise ValueError('The input must be a valid integer number.')
        
        if self.__uqtype == 'covariance':
            # Generate Gaussian distribution based on covariance matrix
            sig = np.sqrt(self.covmat[n,n])
            xmean = self.mean[n]
            x = np.linspace(xmean-4*sig,xmean+4*sig,500)
            pdf = 1/sig/np.sqrt(2*np.pi)*np.exp(-((x-xmean)/sig)**2/2)
            
            # Clip the distributions at outside the boundaries
            pdf[x < self.__lb[n]] = 0
            pdf[x > self.__ub[n]] = 0

        if self.__uqtype == 'bootstrap':
            # Get bw using silverman's rule (1D only)
            sigma = np.std(self.__samples[:, n], ddof=1)
            bw = sigma * (len(self.__samples) * 3 / 4.0) ** (-1 / 5)

            # Make histogram
            bins = np.linspace(self.__lb[n], self.__ub[n], 2 ** 10 + 1)
            count, edges = np.histogram(self.__samples[:, n], bins=bins)

            # Generate kernel
            delta = (edges.max() - edges.min()) / (len(edges) - 1)
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
        if p>100 or p<0:
            raise ValueError('The input must be a number between 0 and 100')
 
        x = np.zeros(self.__nParam)
        for n in range(self.__nParam):
            # Get parameter PDF
            values,pdf = self.pardist(n)
            # Compute corresponding CDF
            cdf = np.cumsum(pdf)
            # Eliminate duplicates
            cdf, index = np.lib.arraysetops.unique(cdf,return_index=True)
            # Interpolate requested percentile
            x[n] = np.interp(p/100,cdf,values[index])

        return x
    #--------------------------------------------------------------------------------


    # Covariance-based confidence intervals
    #--------------------------------------------------------------------------------
    def ci(self,coverage):
        
        if coverage>100 or coverage<0:
            raise ValueError('The input must be a number between 0 and 100')
        
        alpha = 1 - coverage/100
        p = 1 - alpha/2 # percentile
        
        x = np.zeros((self.__nParam,2))
        if self.__uqtype=='covariance':
                # Compute covariance-based confidence intervals
                # Clip at specified box boundaries
                x[:,0] = np.maximum(self.__lb, self.__parfit - norm.ppf(p)*np.sqrt(np.diag(self.__covmat)))
                x[:,1] = np.minimum(self.__ub, self.__parfit + norm.ppf(p)*np.sqrt(np.diag(self.__covmat)))
                
        elif self.__uqtype=='bootstrap':
                # Compute bootstrap-based confidence intervals
                # Clip possible artifacts from the percentile estimation
                x[:,0] = np.minimum(self.percentile(p*100), np.amax(self.__samples))
                x[:,1] = np.maximum(self.percentile((1-p)*100), np.amin(self.__samples))

        return x



    # Error Propagation (covariance-based only)
    #--------------------------------------------------------------------------------
    def propagate(self,model,lbm=[],ubm=[]):
        
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
        J,_ = jacobianest(model,parfit)
        
        # Clip at boundaries
        modelfit = np.maximum(modelfit,lbm)
        modelfit = np.minimum(modelfit,ubm)
        
        # Error progation
        modelcovmat = J@self.__covmat@J.T
        
        # Construct new CI-structure for the model
        return  uqst('covariance',modelfit,modelcovmat,lbm,ubm)
    #--------------------------------------------------------------------------------

