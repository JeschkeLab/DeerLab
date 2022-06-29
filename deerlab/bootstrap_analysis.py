# bootstrap_analysis.py - Bootstrap analysis for uncertainty estimation
# --------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md).
# Copyright(c) 2019-2022: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import types
from tqdm.auto import tqdm
from joblib import Parallel, delayed
from deerlab.classes import UQResult
from deerlab.utils import isnumeric
from deerlab.noiselevel import noiselevel

def bootstrap_analysis(fcn,Vexp,Vfit, samples=1000, noiselvl=None, resampling='gaussian', verbose = False, cores=1, memorylimit=8):
    r""" 
    Bootstrap analysis for uncertainty quantification

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

    noiselvl : scalar or list thereof
        Noise level of the input dataset(s), specified as standard deviations.
        If not specified, these are automatically estimated from the experimental
        dataset(s). 

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

    memorylimit : 
        Memory limit to be allocated for the full boostrap analysis. If the requested analysis exceeds this limit, the 
        execution will be stopped. The default is 12GB.  

    Returns
    -------
    bootuq : :ref:`UQResult` or list of :ref:`UQResult`
        Bootstrap uncertainty quantification for each variable returned by ``fcn``. 


    Examples
    --------
    To analyze several variables during the same run, the function must return tuple 
    containing all of them::

        def myfcn(V):
            Pfit1 = dl.nlls(V,dl.dd_gauss,r,K)
            Pfit2 = dl.nlls(V,dl.dd_randcoil,r,K)
            return Pfit1,Pfit2

        bootuq = bootstrap_analysis(myfcn,V,Vfit)
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
            raise KeyError(f'The Vexp and Vfit #{i} must have the same number of elements.')
    if type(fcn) is not types.FunctionType:
        raise KeyError('The 1st argument must be a callable function accepting dataset(s) as input.')

    # Get residuals and estimate standard deviation
    for i in range(nSignals):
        residuals = [Vfit[i] - Vexp[i] for i  in range(nSignals)]
    if noiselvl is None:
        noiselvl = [noiselevel(Vexp[i]) for i in range(nSignals)]
    noiselvl = np.atleast_1d(noiselvl)

    # Prepare bootstrap sampler (reduced I/O-stream when parallelized) 
    def sample():
        Vsample = [0]*nSignals
        for i in range(nSignals):
            #Determine requested re-sampling method
            if resampling == 'gaussian':
                # Resample from a Gaussian distribution with variance estimated from the residuals
                Vsample[i] = Vfit[i] + np.random.normal(0, noiselvl[i], len(Vfit[i]))
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
        if not all(isnumeric(x) for x in var):
            raise ValueError('Non-numeric output arguments by the analyzed function are not accepted.')
    
    # Check that the full bootstrap analysis will not exceed the memory limits
    memory_requirements = np.sum([nSamples*np.atleast_1d(var).size*np.atleast_1d(var).itemsize])/1e9 # GB
    if memory_requirements > memorylimit: 
        raise MemoryError(f'The requested bootstrap analysis requires {memory_requirements:.2f}GB, exceeding the current memory limit {memorylimit:.2f}GB.')

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
    if verbose : print(f'Bootstrap analysis with {cores} cores:')          
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
        stats.append(UQResult('bootstrap',variables))
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