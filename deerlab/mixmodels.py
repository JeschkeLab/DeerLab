# mixmodels.py - Parametric model mixer
# ---------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

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
    If one mixes a single Gaussian model (2 parameters) with a worm-like chain (WLC) model (2 parameters) into a single model::

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
    for i,model in enumerate(models):
        nparam = len(model.start)
        pidx.append(idx + np.arange(0,nparam))
        idx = idx + nparam
        for j in range(nparam):
            Info['Parameters'].append(f'Model {i+1}: {model.parameters[j]}')
            Info['Units'].append(model.units[j])
            Info['Lower'].append(model.lower[j])
            Info['Upper'].append(model.upper[j])
            Info['Start'].append(model.start[j])

        # Add amplitudes for each model
        Info['Parameters'].append(f'Model {i+1}: Amplitude')
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


    # =================================================================
    def setmetadata(parameters,units,start,lower,upper):
        """
        Decorator: Set model metadata as function attributes 
        """
        def _setmetadata(func):
            func.parameters = parameters
            func.units = units
            func.start = start
            func.lower = lower
            func.upper = upper
            return func
        return _setmetadata
    # =================================================================


    # =================================================================
    @setmetadata(
    parameters = Info['Parameters'],
    units = Info['Units'],
    start = Info['Start'],
    lower = Info['Lower'],
    upper = Info['Upper'])
    def mixedFunction(ax, params):
        """
        Mixed model function handle
        ---------------------------
        Function to allow request of information structure or model values
        """

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