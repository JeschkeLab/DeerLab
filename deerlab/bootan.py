# bootan.py - Bootstrap analysis for uncertainty estimation
# --------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
from joblib import Parallel, delayed
import types
from deerlab.classes import UncertQuant

def bootan(fcn,Vexp,Vfit, samples=1000, resampling='gaussian', verbose=False, cores=1):
    r""" Bootstrap analysis for uncertainty quantification

    Parameters
    ----------
    fcn : callable
        Function to be analyzed. Must be a callable function accepting a signal
        array as input and returning a tuple with all variables to be analyzed.
        All variables must be numerical arrays (no strings or booleans) and 
        must preserve shape between calls.

    Vexp : array_like or list of array_like
        Experimental dataset(s).

    Vfit : array or list of array_like
        Fit of the dataset(s).

    samples : scalar
        Number of bootstrap samples to analyze. The quality of bootstrapping 
        results improve with the number of boostrap samples evaluated, the 
        default is 1000.
        
    Returns
    -------
    bootuq : :ref:`UncertQuant` or list of :ref:`UncertQuant`
        Bootstrap uncertainty quantification for each variable returned by ``fcn``. 

    Other Parameters
    ----------------
    resampling : string
        Specifies the method employed for re-sampling new bootstrap samples.

        * ``'gaussian'`` - Sample noise from a Gaussian distribution.
        * ``'residual'`` - Sample noise from the fit residuals.
        The default is ``'gaussian'``.

    cores : scalar
        Number of CPU cores/processes for parallel computing. If ``cores=1`` no parallel 
        computing is used. If ``cores=-1`` all available CPUs are used. The default is 
        one core (no parallelization).

    verbose : boolean
        Specifies whether to print the progress of the bootstrap analysis on the 
        command window, the default is false.

    Examples
    --------
    To analyze several variables during the same run, the function must return tuple 
    containing all of them::

        def myfcn(V):
            Pfit1 = fitparamodel(V,dd_gauss,r,K)
            Pfit2 = fitparamodel(V,dd_randcoil,r,K)
            return Pfit1,Pfit2

        bootuq = bootan(@(V)myfcn(p,varargin),V,Vfit)
        Pfit1_uq = bootuq[0]
        Pfit2_uq = bootuq[1]


    """

    if samples<2:
        raise ValueError('At least two bootstrap samples are required.')
    nSamples = int(samples)

    # Validation of input signals
    Vexp,Vfit = ( [V] if type(V) is np.ndarray else V for V in [Vexp,Vfit])
    if len(Vexp)!=len(Vfit):
        raise KeyError('The same number of signals V and fits Vfit must be provided.')
    nSignals = len(Vexp)
    for i in range(len(Vfit)):
        if len(Vexp[i])!=len(Vfit[i]):
            raise KeyError('The V and Vfit #{} must have the same number of elements.'.format(i))
    if type(fcn) is not types.FunctionType:
        raise KeyError('The 1st argument must be a callable function accepting dataset(s) as input.')

    # Get residuals and estimate standard deviation
    residuals,sigma = ([],[])
    for i in range(nSignals):
        residuals.append(Vfit[i] - Vexp[i])
        sigma.append(np.std(residuals[i]))

    # Use original data on the first run
    varargout = fcn(Vexp)
    if type(varargout) is not tuple:
        varargout = (varargout,)
    nout = len(varargout)

    # Assert that all outputs are strictly numerical
    for var in varargout:
        if not all((isinstance(x,(float,int)) for x in var)):
            raise ValueError('Non-numeric output arguments by the analyzed function are not accepted.')
    
    shapeout = []
    evals = []
    for j in range(nout):
        out = np.atleast_1d(varargout[j])
        shapeout.append(np.shape(out))
        evals.append(out[np.newaxis,:])

    # Generate all bootstrap samples
    Vsamples = []
    for iSample in range(nSamples-1):
        sample = [0]*nSignals
        # Inform of progress if requested
        if verbose:
            print('Bootstrapping: #{}/#{} samples finished'.format(iSample+1,nSamples), end='\r', flush=True)
        
        for i in range(nSignals):
            #Determine requested re-sampling method
            if resampling == 'gaussian':
                # Resample from a Gaussian distribution with variance estimated from the residuals
                sample[i] = Vfit[i] + np.random.normal(0, sigma[i], len(Vfit[i]))

            elif resampling == 'residual':
                # Resample from the residual directly
                sample[i] =  Vfit[i] + residuals[i][np.random.permutation(len(Vfit[i]))]
            else:
                raise KeyError("Resampling method not found. Must be either 'gaussian' or 'residual'.")
        Vsamples.append(sample)
        
    # Run analysis function on each sample in series (cores==1) or in parallel (cores>1)    
    out = Parallel(n_jobs=cores)(delayed(fcn)(Vsample) for Vsample in Vsamples)

    # Post-process the outputs
    for varargout in out:
        if type(varargout) is not tuple:
            varargout = (varargout,)
        nout = len(varargout)

        # Loop over the different outputs of the function
        for j in range(nout):
            out = np.atleast_1d(varargout[j])
            if np.shape(out)!=shapeout[j]: 
                raise TypeError('Inconsistent output variable size. One of the outputs of the analyzed function is changing its size in between runs. To solve this, fix the axis of the output and interpolate the result.')
            # Store outputs in an N-dimensional array
            evals[j] = np.concatenate((evals[j],out[np.newaxis,:]),0)
            
        
    # Compile statistics for all parameters from bootstrap samples
    #-------------------------------------------------------------------------------
    stats = []
    for bootsamples in evals:
        stats.append(UncertQuant('bootstrap',bootsamples))
    if len(stats)==1:
        stats = stats[0]
    return stats