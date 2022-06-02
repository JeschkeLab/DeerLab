# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2022: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
from deerlab.utils import Jacobian, nearest_psd
from scipy.stats import norm
from scipy.signal import fftconvolve
from scipy.linalg import block_diag
from scipy.optimize import brentq
from scipy.interpolate import interp1d

class FitResult(dict):
# ========================================================================
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

    Methods
    -------
    plot()
        Display the fit results on a Matplotlib window. The script returns a 
        `matplotlib.axes <https://matplotlib.org/api/axes_api.html>`_ object.
        All graphical parameters can be adjusted from this object.

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

    def __str__(self): 
        return self._summary

    def __repr__(self):
        if self.keys():
            m = max(map(len, list(self.keys()))) + 1
            return '\n'.join([k.rjust(m) + ': ' + repr(v)
                              for k, v in sorted(self.items())])
        else:
            return self.__class__.__name__ + "()"

    def __dir__(self):
        return list(self.keys())
# =========================================================================

class UQResult:
# =========================================================================
    r""" Represents the uncertainty quantification of fit results.

    Attributes
    ----------
    type : string
        Uncertainty quantification approach:

            * 'covariance' - Covariance-based uncertainty analysis
            * 'bootstrap' - Bootstrapped uncertainty analysis
    
    mean : ndarray
        Mean values of the uncertainty distribution of the parameters.
    median : ndarray
        Median values of the uncertainty distribution of the parameters.
    std : ndarray
        Standard deviations of the uncertainty distribution of the parameters.
    covmat : ndarray
        Covariance matrix
    nparam : int scalar
        Number of parameters in the analysis.

    Methods
    -------

    """

    def __init__(self,uqtype,data=None,covmat=None,lb=None,ub=None,threshold=None,profiles=None,noiselvl=None):

        #Parse inputs schemes
        if uqtype=='covariance':
            # Scheme 1: UQResult('covariance',parfit,covmat,lb,ub)
            self.type = uqtype
            parfit = data
            nParam = len(parfit)
            
        elif uqtype == 'profile':
            # Scheme 2: UQResult('profile',profiles)
            if not isinstance(profiles,list):
                profiles = [profiles]
            self.type = uqtype
            self.__parfit = data
            self.__noiselvl = noiselvl
            self.profile = profiles
            self.threshold = threshold
            nParam = len(np.atleast_1d(data))

        elif uqtype == 'bootstrap':
            # Scheme 2: UQResult('bootstrap',samples)
            self.type = uqtype
            samples = data
            self.samples = samples
            nParam = np.shape(samples)[1]            

        elif uqtype=='void':
            # Scheme 2: UQResult('void')
            self.type = uqtype
            self.mean, self.median, self.std, self.covmat, self.nparam = ([] for _ in range(5))
            return
        else:
            raise NameError('uqtype not found. Must be: ''covariance'', ''bootstrap'' or ''void''.')

        if lb is None:
            lb = np.full(nParam, -np.inf)
        
        if ub is None:
            ub = np.full(nParam, np.inf)
            
        # Set private variables
        self.__lb = lb
        self.__ub = ub
        self.nparam = nParam

        # Create confidence intervals structure
        if uqtype=='covariance':
            self.mean = parfit
            self.median = parfit
            self.std = np.sqrt(np.diag(covmat))
            self.covmat = covmat

        # Profile-based CI specific fields
        elif uqtype == 'profile':
            xs = [self.pardist(n)[0] for n in range(nParam)]
            pardists = [self.pardist(n)[1] for n in range(nParam)]
            means = [np.trapz(pardist*x,x) for x,pardist in zip(xs,pardists)]
            std = [np.sqrt(np.trapz(pardist*(x-mean)**2,x)) for x,pardist,mean in zip(xs,pardists,means)]
            self.mean = means
            self.median = self.percentile(50)
            self.std = std
            self.covmat = np.diag(np.array(std)**2)
                
        # Bootstrap-based CI specific fields
        elif uqtype == 'bootstrap':
            means = np.mean(samples,0)
            covmat = np.squeeze(samples).T@np.squeeze(samples)/np.shape(samples)[0] - means*means.T
            self.mean = means
            self.median = self.percentile(50)
            self.median = np.array([nth_samples[0] if np.all(nth_samples==nth_samples[0]) else self.median[n] for n,nth_samples in enumerate(samples.T)])
            self.std = np.squeeze(np.std(samples,0))
            self.covmat = covmat


    # Gets called when an attribute is accessed
    #--------------------------------------------------------------------------------
    def __getattribute__(self, attr):
        try:
            # Calling the super class to avoid recursion
            if (attr!='type' and not '__' in attr) and super(UQResult, self).__getattribute__('type') == 'void':
                # Check if the uncertainty quantification has been done, if not report that there is nothing in the object
                raise ValueError(f"The requested attribute/method ('{attr}') is not available. Uncertainty quantification has not been calculated during the fit.")
        except AttributeError:
            # Catch cases where 'type' attribute has still not been defined (e.g. when using copy.deepcopy)
            pass
        # Otherwise return requested attribute
        return super(UQResult, self).__getattribute__(attr)
    #--------------------------------------------------------------------------------


    # Combination of multiple uncertainties
    #--------------------------------------------------------------------------------
    def join(self,*args):
        """
        Combine multiple uncertainty quantification instances.

        Parameters
        ----------
        uq : any number of :ref:`UQResult`
            Uncertainty quantification objects with ``N1,N2,...,Nn`` parameters to be joined 
            to the object calling the method with ``M`` parameters. 
        
        Returns
        -------
        uq_joined : :ref:`UQResult`
            Joined uncertainty quantification object with a total of ``M + N1 + N2 + ... + Nn`` parameters. 
            The parameter vectors are concatenated on the order they are passed. 
        """
        newargs = []
        if self.type=='covariance':
            # Original metadata
            newargs.append(self.mean)
            newargs.append(self.covmat)
            newargs.append(self.__lb)
            newargs.append(self.__ub)
        elif self.type=='bootstrap':
            newargs.append(self.samples)
            
        for uq in args:
            if not isinstance(uq, UQResult):
                raise TypeError('Only UQResult objects can be joined.')
            if uq.type!=self.type:
                raise TypeError(f'A UQResult of type ({uq.type}) and another of type ({self.type}) cannot be joined.')
            if uq.type=='void':
                raise TypeError('Void UQResults cannot be joined.')
            # Concatenate metadata of external UQResult objects
            if self.type=='covariance':
                newargs[0] = np.concatenate([newargs[0], uq.mean])
                newargs[1] = block_diag(newargs[1], uq.covmat)
                newargs[2] = np.concatenate([newargs[2], uq.__lb])
                newargs[3] = np.concatenate([newargs[3], uq.__ub])
            elif self.type=='bootstrap':
                newargs[0] = np.concatenate([newargs[0],uq.samples],axis=1) 
                
        # Return new UQResult object with combined information    
        return UQResult(self.type,*newargs) 
    #--------------------------------------------------------------------------------


    # Parameter distributions
    #--------------------------------------------------------------------------------
    def pardist(self,n=0):
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
        isdelta = False

        if self.type == 'covariance':
            # Generate Gaussian distribution based on covariance matrix
            sig = np.sqrt(self.covmat[n,n])
            xmean = self.mean[n]
            x = np.linspace(xmean-4*sig,xmean+4*sig,500)
            pdf = 1/sig/np.sqrt(2*np.pi)*np.exp(-((x-xmean)/sig)**2/2)
            
        if self.type == 'bootstrap':
            # Get bw using silverman's rule (1D only)
            samplen = self.samples[:, n].real

            # Take limited precision into account to avoid round-off errors
            samplen = np.round(samplen,200)

            if np.all(samplen == samplen[0]):
                # Dirac's delta distribution 
                x = np.array([0.9*samplen[0],samplen[0],1.1*samplen[0]])
                pdf = np.array([0,1,0]).astype(float)
                isdelta = True 
            else:
                sigma = np.std(samplen, ddof=1)
                bw = sigma*(len(samplen)*3/4.0)**(-1/5)

                # Make histogram
                maxbin = np.maximum(np.max(samplen),np.mean(samplen)+3*sigma)
                minbin = np.minimum(np.min(samplen),np.mean(samplen)-3*sigma)
                bins = np.linspace(minbin,maxbin, 2**10 + 1)
                count, edges = np.histogram(samplen, bins=bins)

                # Generate kernel
                delta = np.maximum(np.finfo(float).eps,(edges.max() - edges.min()) / (len(edges) - 1))
                kernel_x = np.arange(-4*bw, 4*bw + delta, delta)
                kernel = norm(0, bw).pdf(kernel_x)

                # Convolve
                pdf = fftconvolve(count, kernel, mode='same')
                
                # Set x coordinate of pdf to midpoint of bin
                x = edges[:-1] + delta

        if self.type=='profile':
            if not isinstance(self.profile,list) and n==0:
                profile = self.profile
            else:
                profile = self.profile[n] 
            σ = self.__noiselvl
            obj2likelihood = lambda f: 1/np.sqrt(σ*2*np.pi)*np.exp(-1/2*(f-np.min(f))/σ**2)
            profileinterp = interp1d(profile['x'], profile['y'], kind='slinear', fill_value=1e6,bounds_error=False)
            x = np.linspace(np.min(profile['x']), np.max(profile['x']), 2**10 + 1)
            pdf = obj2likelihood(profileinterp(x)) 

            # Generate kernel
            sigma = np.sum(x*pdf/np.sum(pdf))
            bw = sigma*(1e12*3/4.0)**(-1/5)
            delta = np.maximum(np.finfo(float).eps,(x.max() - x.min()) / (len(x) - 1))
            kernel_x = np.arange(-5*bw, 5*bw + delta, delta)
            kernel = norm(0, bw).pdf(kernel_x)

            # Convolve
            pdf = fftconvolve(pdf, kernel, mode='same')

        # Clip the distributions outside the boundaries
        pdf[x < self.__lb[n]] = 0
        pdf[x > self.__ub[n]] = 0

        # Enforce non-negativity (takes care of negative round-off errors)
        pdf = np.maximum(pdf,0)

        # Ensure normalization of the probability density function (if not a Dirac delta function)
        if not isdelta:
            if np.trapz(pdf, x)!=0:
                pdf = pdf/np.trapz(pdf, x)
        
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
        ci : 2D-ndarray
            Confidence intervals for the parameters:

                * ``ci[:,0]`` - Lower confidence intervals
                * ``ci[:,1]`` - Upper confidence intervals
        """
        if coverage>100 or coverage<0:
            raise ValueError('The input must be a number between 0 and 100')
        value = self.mean if hasattr(self,'mean') else self.__parfit
        iscomplex = np.iscomplexobj(value)
        
        alpha = 1 - coverage/100
        p = 1 - alpha/2 # percentile
        
        confint = np.zeros((self.nparam,2))
        if iscomplex: confint = confint.astype(complex)

        if self.type=='covariance':
            # Compute covariance-based confidence intervals
            # Clip at specified box boundaries
            standardError = norm.ppf(p)*np.sqrt(np.diag(self.covmat))
            confint[:,0] = np.maximum(self.__lb, self.mean.real - standardError)
            confint[:,1] = np.minimum(self.__ub, self.mean.real + standardError)
            if iscomplex:
                confint[:,0] = confint[:,0] + 1j*np.maximum(self.__lb, self.mean.imag - standardError)
                confint[:,1] = confint[:,1] + 1j*np.minimum(self.__ub, self.mean.imag + standardError)

        elif self.type=='bootstrap':
            # Compute bootstrap-based confidence intervals
            # Clip possible artifacts from the percentile estimation
            confint[:,0] = np.minimum(self.percentile((1-p)*100), np.amax(self.samples))
            confint[:,1] = np.maximum(self.percentile(p*100), np.amin(self.samples))

        elif self.type=='profile':
            # Compute likelihood-profile-based confidence intervals
            for n,profile in enumerate(self.profile):
                
                # Construct interpolator for the profile
                profileinterp = interp1d(profile['x'], profile['y'], kind='slinear', fill_value=1e6,bounds_error=False)

                #-----------------------------------------------------------------
                def getCIbound(boundary,optimum):
                    def getprofile_at(value):
                        return profileinterp(value) - self.threshold(coverage/100)
                    # Evaluate the profile function
                    fbound = getprofile_at(boundary)
                    f0 = getprofile_at(optimum)
                    
                    # Check the signs of the shifted profile
                    if np.sign(fbound)==np.sign(f0):
                        # If both edges have the same sign return one of the edges
                        ci_bound = boundary
                    else:
                        searchrange = [boundary,optimum] if boundary<optimum else [optimum,boundary]
                        ci_bound = brentq(getprofile_at, *searchrange,maxiter=int(1e4))
                    return ci_bound
                #-----------------------------------------------------------------
                # Get the upper and lower bounds of the confidence interval
                confint[n,0] = getCIbound(profile['x'].min(),self.__parfit[n])
                confint[n,1] = getCIbound(profile['x'].max(),self.__parfit[n])

        # Remove singleton dimensions
        confint = np.squeeze(confint)
        return confint



    # Error Propagation (covariance-based only)
    #--------------------------------------------------------------------------------
    def propagate(self,model,lb=None,ub=None,samples=None):
        """
        Uncertainty propagation. This function takes the uncertainty analysis of the 
        parameters and propagates it to another function depending on those parameters.

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
        modeluq : :ref:`UQResult`
            New uncertainty quantification analysis for the outputs of ``model``.
        """
        parfit = self.mean
        # Evaluate model with fit parameters
        modelfit = model(parfit)
        iscomplex = np.iscomplexobj(modelfit)

        # Validate input boundaries
        if lb is None:
            lb = np.full(np.size(modelfit), -np.inf)

        if ub is None:
            ub = np.full(np.size(modelfit), np.inf)

        lb,ub = (np.atleast_1d(var) for var in [lb,ub])

        if np.size(modelfit)!=np.size(lb) or np.size(modelfit)!=np.size(ub):
            raise IndexError ('The 2nd and 3rd input arguments must have the same number of elements as the model output.')

        if samples is None:
            Nsamples = 1000
        else:
            Nsamples = samples

        if self.type=='covariance':

            if iscomplex:
                model_ = model 
                model = lambda p: np.concatenate([model_(p).real,model_(p).imag])

            # Get jacobian of model to be propagated with respect to parameters
            J = Jacobian(model,parfit,self.__lb,self.__ub)

            # Clip at boundaries
            modelfit = np.maximum(modelfit,lb)
            modelfit = np.minimum(modelfit,ub)

            # Error progation
            modelcovmat = nearest_psd(J@self.covmat@J.T)

            if iscomplex:
                N = modelcovmat.shape[0]
                Nreal = np.arange(0,N/2).astype(int)
                Nimag = np.arange(N/2,N).astype(int)
                modelcovmat = modelcovmat[np.ix_(Nreal,Nreal)] + 1j* modelcovmat[np.ix_(Nimag,Nimag)]

            # Construct new uncertainty object
            return  UQResult('covariance',modelfit,modelcovmat,lb,ub)

        elif self.type=='bootstrap':

            sampled_parameters = [[] for _ in range(self.nparam)]
            for n in range(self.nparam):
                # Get the parameter uncertainty distribution
                values,pdf = self.pardist(n)
                # Random sampling form the uncertainty distribution
                np.random.seed(0)
                sampled_parameters[n] =  [np.random.choice(values, p=pdf/sum(pdf)) for _ in range(Nsamples)]
            # Convert to matrix
            sampled_parameters = np.atleast_2d(sampled_parameters)

            # Bootstrap sampling of the model response
            sampled_model = [model(sampled_parameters[:,n]) for n in range(Nsamples)]
            sampled_model = [np.maximum(sample,lb) for sample in sampled_model]
            sampled_model = [np.minimum(sample,ub) for sample in sampled_model]

            # Convert to matrix
            sampled_model = np.atleast_2d(sampled_model)
            
            # Construct new uncertainty object
            return  UQResult('bootstrap',data=sampled_model,lb=lb,ub=ub)
 
    #--------------------------------------------------------------------------------

# =========================================================================
