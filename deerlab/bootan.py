# bootan.py - Bootstrap analysis for uncertainty estimation
# --------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import types
from deerlab import uqst

def bootan(fcn,Vexp,Vfit, samples=1000, resampling='gaussian', verbose = False):
    """
    Bootstrap analysis for uncertainty quantification
    =================================================
    
    Bootstrap analysis of for uncertainty quantification.

    Usage:
    ------
    bootuq = bootan(fcn,V,Vfit,samples)

    Arguments:
    ----------
    fcn (callable)
        Function to be analyzed.
    V (array or list of arrays)
        Experimental dataset(s).
    Vfit (array or list of arrays)
        Fit of the dataset(s).
    samples (scalar, default=1000)
        Number of bootstrap samples to analyze.

    Returns:
    --------
    bootuq (obj)
        Bootstrap uncertainty quantification object.

    Additional keyword arguments:
    ----------------------------- 
    resampling (string, default='gaussian')
        Specifies the method employed for re-sampling new bootstrap samples.
            'gaussian' - Sample noise from a Gaussian distribution
            'residual' - Sample noise from the fit residuals
    verbose (boolean, default=False)
        Specifies whether to print the progress of the bootstrap analysis on the command window.
    """

    nSamples = samples

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


    # Bootsrap analysis
    #-------------------------------------------------------------------------------
    Vresample = [0]*nSignals
    for iSample in range(nSamples-1):
        
        # Change the random seed each iteration, but keeping reproducibility between runs
        np.random.seed(iSample)

        # Inform of progress if requested
        if verbose:
            print('\r Bootstrapping: #{}/#{} samples finished'.format(iSample,nSamples))
            
        for i in range(nSignals):

            # Get a bootstrap sample
            if iSample>0:
                
                #Determine requested re-sampling method
                if resampling is 'gaussian':
                    # Resample from a Gaussian distribution with variance estimated from the residuals
                    Vresample[i] = Vfit[i] + np.random.normal(0, sigma[i], len(Vfit[i]))

                elif resampling is 'residual':
                    # Resample from the residual directly
                    Vresample[i] =  Vfit[i] + residuals[i][np.random.permutation(len(Vfit[i]))]
                else:
                    raise KeyError("Resampling method not found. Must be either 'gaussian' or 'residual'.")
            else:
                # Use original data on the first run
                Vresample[i] = Vexp[i]

        # Run the model function with bootstrap sample
        varargout = fcn(Vresample)
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
        stats.append(uqst('bootstrap',bootsamples))
    if len(stats)==1:
        stats = stats[0]
    return stats