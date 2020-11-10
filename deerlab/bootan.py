# bootan.py - Bootstrap analysis for uncertainty estimation
# --------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import types
from tqdm.auto import tqdm
from joblib import Parallel, delayed
from deerlab.classes import UncertQuant

def bootan(fcn,Vexp,Vfit, samples=1000, resampling='gaussian', verbose = False, cores=1):
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

    samples : scalar, optional
        Number of bootstrap samples to analyze. The quality of bootstrapping 
        results improve with the number of boostrap samples evaluated, the 
        default is 1000.

    resampling : string, optional
        Specifies the method employed for re-sampling new bootstrap samples.

        * ``'gaussian'`` - Sample noise from a Gaussian distribution.
        * ``'residual'`` - Sample noise from the fit residuals.
        The default is ``'gaussian'``.

    cores : scalar, optional
        Number of CPU cores/processes for parallel computing. If ``cores=1`` no parallel 
        computing is used. If ``cores=-1`` all available CPUs are used. The default is 
        one core (no parallelization).

    verbose : boolean, optional
        Specifies whether to print the progress of the bootstrap analysis on the 
        command window, the default is false.


    Returns
    -------
    bootuq : :ref:`UncertQuant` or list of :ref:`UncertQuant`
        Bootstrap uncertainty quantification for each variable returned by ``fcn``. 


    Examples
    --------
    To analyze several variables during the same run, the function must return tuple 
    containing all of them::

        def myfcn(V):
            Pfit1 = dl.fitparamodel(V,dl.dd_gauss,r,K)
            Pfit2 = dl.fitparamodel(V,dl.dd_randcoil,r,K)
            return Pfit1,Pfit2

        bootuq = bootan(myfcn,V,Vfit)
        Pfit1_uq = bootuq[0]
        Pfit2_uq = bootuq[1]


    """

    # Catch invalid sample sizes
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

    # Prepare bootstrap sampler (reduced I/O-stream when parallelized) 
    def sample():
        Vsample = [0]*nSignals
        for i in range(nSignals):
            #Determine requested re-sampling method
            if resampling == 'gaussian':
                # Resample from a Gaussian distribution with variance estimated from the residuals
                Vsample[i] = Vfit[i] + np.random.normal(0, sigma[i], len(Vfit[i]))
            elif resampling == 'residual':
                # Resample from the residual directly
                Vsample[i] =  Vfit[i] + residuals[i][np.random.permutation(len(Vfit[i]))]
            else:
                raise KeyError("Resampling method not found. Must be either 'gaussian' or 'residual'.")
        return Vsample       

    # Use original data on the first run
    varargout = fcn(Vexp)
    if type(varargout) is not tuple:
        varargout = (varargout,)
    nout = len(varargout)

    # Assert that all outputs are strictly numerical
    for var in varargout:
        if not all((isinstance(x,(float,int)) for x in var)):
            raise ValueError('Non-numeric output arguments by the analyzed function are not accepted.')
    
    # Get ndarray shapes of all outputs
    shapeout = []
    evals = []
    for j in range(nout):
        out = np.atleast_1d(varargout[j])
        shapeout.append(np.shape(out))
        evals.append(out[np.newaxis,:])

    # Main bootstrap function
    def bootsample():
        # Get bootstrap sample
        Vsample = sample()
        # Run the analysis function
        out = fcn(Vsample)
        return out

    # Run bootsample() for all samples in series (cores=1) or in parallel (cores>1)
    if verbose : print('Bootstrap analysis with {0} cores:'.format(cores))          
    out = _ProgressParallel(n_jobs=cores,total=nSamples,use_tqdm=verbose)(delayed(bootsample)() for _ in range(nSamples))
    
    # Post-process the outputs
    for varargout in out:
        if type(varargout) is not tuple:
            varargout = (varargout,)
        nout = len(varargout)

        # Loop over the different outputs of the function
        for j in range(nout):
            out = np.atleast_1d(varargout[j])
            if np.shape(out)!=shapeout[j]: 
                raise TypeError('Inconsistent output variable shape. One of the outputs of the analyzed function is changing its size/shape in between runs.')
            # Store outputs in an N-dimensional array
            evals[j] = np.concatenate((evals[j],out[np.newaxis,:]),0)
            
        
    # Compile statistics for all parameters from bootstrap samples
    #-------------------------------------------------------------------------------
    stats = []
    for variables in evals:
        stats.append(UncertQuant('bootstrap',variables))
    if len(stats)==1:
        stats = stats[0]
    return stats


#-------------------------------------------------------------------------------
class _ProgressParallel(Parallel):
    """
    Patch for joblib.Parallel

    Overrides the print_progress() method to enable the synchronous use of the TQDM bar
    even for parallel processing.  
    """
    def __init__(self, use_tqdm=True, total=None, *args, **kwargs):
        self._use_tqdm = use_tqdm
        self._total = total
        super().__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        with tqdm(disable=not self._use_tqdm, total=self._total) as self._pbar:
            return Parallel.__call__(self, *args, **kwargs)

    def print_progress(self):
        if self._total is None:
            self._pbar.total = self.n_dispatched_tasks
        self._pbar.n = self.n_completed_tasks
        self._pbar.refresh()
#-------------------------------------------------------------------------------