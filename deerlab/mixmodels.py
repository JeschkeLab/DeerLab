# mixmodels.py - Parametric model mixer
# ---------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2020: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
import types

def mixmodels(*models):
    r""" 
    Combine DeerLab parametric models into a mixed multi-component model.

    Parameters
    ----------
    models : callables
        Parametric DeerLab models to be combined.

    Returns
    --------
    newmodel : callable
        Mixed model.   

    Examples
    --------
    If one mixes a single Gaussian model (2 parameters) with a WLC model (2 parameters) into a single model::

        newmodel = mixmodels(dd_gauss,dd_wormchain)


    the resulting model newmodel will contain 6 parameters in the following order: the 2 single-Gaussian parameters,
    1 amplitude parameter for the Gaussian model, the 2 WLC parameters, and 1 amplitude parameter for the WLC model.
    """
    if len(models)==0:
        raise KeyError('At least one model must be provided.')

    if np.any([type(model) is not types.FunctionType for model in models]):
        raise TypeError('Input arguments must all be function handles.')

    # Detemine number of models to be mixed
    nModels = len(models)

    # Combine info structures from all models
    idx = 0
    Info = dict(Parameters=[],Units=[],Start=[],Lower=[],Upper=[])
    pidx = []
    pidx_amp = []
    for i in range(nModels):
        info = models[i]()
        nparam = len(info['Start'])
        pidx.append(idx + np.arange(0,nparam))
        idx = idx + nparam
        for j in range(nparam):
            Info['Parameters'].append('Model {}: {}'.format(i+1,info['Parameters'][j]))
            Info['Units'].append(info['Units'][j])
            Info['Lower'].append(info['Lower'][j])
            Info['Upper'].append(info['Upper'][j])
            Info['Start'].append(info['Start'][j])

        # Add amplitudes for each model
        Info['Parameters'].append('Model {}: Amplitude'.format(i+1))
        Info['Units'].append('')
        Info['Lower'].append(0)
        Info['Upper'].append(1)
        Info['Start'].append(1/nModels)
        pidx_amp.append(len(Info['Start'])-1)
        idx = idx + 1

    # Convert the numerical fields to numpy arrays
    Info['Lower'] = np.asarray(Info['Lower'])
    Info['Upper'] = np.asarray(Info['Upper'])
    Info['Start'] = np.asarray(Info['Start'])

    def mixedFunction(*args):
    # =======================================================================
        """
        Mixed model function handle
        ---------------------------
        Function to allow request of information structure or model values
        """
        if not args:    
            return Info
        
        if len(args)<2:
            raise KeyError('At least two input arguments are required.')
        elif len(args)>3:
            raise KeyError('Only two input arguments are allowed.')
        
        ax, params = args        
        params = np.atleast_1d(params)
        evaled = 0
        for k in range(nModels):
            evaled = evaled + params[pidx_amp[k]]*models[k](ax,params[pidx[k]])
        
        #Normalize the distribution if it is a distance distribution model
        isddmodel = any(['dd' in model.__name__ for model in models])
        if isddmodel and not np.all(evaled==0):
            evaled = evaled/np.trapz(evaled,ax)
        return evaled
    # =======================================================================

    return mixedFunction