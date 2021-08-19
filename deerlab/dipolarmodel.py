# dipolarmodel.py - DeerLab's dipolar EPR model generator
# ---------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
from deerlab.dipolarkernel import dipolarkernel  
from deerlab.model import Model
from deerlab import bg_hom3d

#===============================================================================
def dipolarmodel(t,r,Pmodel=None,Bmodel=bg_hom3d,npathways=1,harmonics=None):

    # Input parsing and validation
    if not isinstance(Pmodel,Model) and Pmodel is not None:
        raise TypeError('The argument Pmodel must be a valid Model object')
    if not isinstance(Bmodel,Model) and Bmodel is not None:
        raise TypeError('The argument Bmodel must be a valid Model object or None.')
    if not isinstance(npathways,int) or npathways<=0:
        raise ValueError('The number of pathway must be an integer number larger than zero.')
    if isinstance(harmonics,float) or isinstance(harmonics,int): 
        harmonics = [int(harmonics)]
    if not isinstance(harmonics,list) and harmonics is not None:
        raise TypeError('The harmonics must be specified as a list of integer values.')    

    #------------------------------------------------------------------------
    def _importparameter(parameter):
        """ Private function for importing a parameter's metadata """
        return {
            'lb' : parameter.lb,
            'ub' : parameter.ub,
            'par0' : parameter.par0,
            'description' : parameter.description,
            'units' : parameter.units,
            'linear' : parameter.linear
        }
    #------------------------------------------------------------------------

    # Parse the harmonics of the dipolar pathways
    if harmonics is None: 
        harmonics = np.ones(npathways)
    if len(harmonics)!=npathways: 
        raise ValueError('The number of harmonics must match the number of dipolar pathways.')
    #------------------------------------------------------------------------
    def dipolarpathways(*param):
        """ Parametric constructor of the dipolar pathways definition """
        param = np.atleast_1d(param)
        if npathways==1:
            # Single-pathway model, use modulation depth notation
            lam,reftime = param
            pathways = [[1-lam],[lam,reftime,harmonics[0]]]
        else:
            # Otherwise, use general notation
            Lam0 = param[0]
            lams = param[np.arange(1,len(param),2)]
            reftimes = param[np.arange(2,len(param),2)]
            # Unmodulated pathways ccontribution
            pathways = [[Lam0]]
            # Modulated pathways
            for n in range(npathways):
                pathways.append([lams[n], reftimes[n], harmonics[n]])    
        return  pathways
    #------------------------------------------------------------------------

    # Construct the signature of the dipolarpathways() function
    if npathways == 1:
        variables = ['mod','reftime']
    else:
        variables = ['Lam0']
        for n in range(npathways):
            variables.append(f'lam{n+1}')
            variables.append(f'reftime{n+1}')

    # Create the dipolar pathways model object
    PathsModel = Model(dipolarpathways,signature=variables)

    if Pmodel is None:
        def nonparametric(r):
            return np.eye(len(r))
        # Create model
        dd_nonparametric = Model(nonparametric,constants='r')
        dd_nonparametric.description = 'Non-parametric distribution model'
        # Parameters
        dd_nonparametric.addlinear('P',vec=len(r),lb=0,par0=0,description='Non-parametric distance distribution')
        Pmodel = dd_nonparametric

    # Populate the basic information on the dipolar pathways parameters
    if npathways==1:
        # Special case: use modulation depth notation instead of general pathway amplitude
        getattr(PathsModel,f'mod').set(lb=0,ub=1,par0=0.5,description=f'Modulation depth',units='')
        getattr(PathsModel,f'reftime').set(par0=0,description=f'Refocusing time',units='μs')
    else:
        # General case: use pathway ampltiudes and refocusing times
        getattr(PathsModel,f'Lam0').set(lb=0,ub=1,par0=0.5,description=f'Amplitude of pathway #{n+1}',units='')
        for n in range(npathways):
            getattr(PathsModel,f'lam{n+1}').set(lb=0,ub=1,par0=0.5,description=f'Amplitude of pathway #{n+1}',units='')
            getattr(PathsModel,f'reftime{n+1}').set(par0=0,description=f'Refocusing time of pathway #{n+1}',units='μs')

    # Construct the signature of the dipolar signal model function
    signature = []
    parameters,linearparam,vecparam = [],[],[]
    for model in [PathsModel,Bmodel,Pmodel]:
        if model is not None:
            for param in model._parameter_list(order='vector'):
                if np.any(getattr(model,param).linear):
                    parameters.append(getattr(model,param))
                    linearparam.append({'name':param,'vec':len(np.atleast_1d(getattr(model,param).idx))})
                elif not (model==Bmodel and param=='lam'):
                    signature.append(param)
                    parameters.append(getattr(model,param))

    # Initialize lists of indices to access subsets of nonlinear parameters 
    Psubset = np.zeros(Pmodel.Nnonlin,dtype=int)
    PathsSubset = np.zeros(PathsModel.Nnonlin,dtype=int)
    if Bmodel is None:
        Bsubset = []
    else:
        Bsubset = np.zeros(Bmodel.Nnonlin-('lam' in Bmodel._parameter_list()),dtype=int)

    # Determine subset indices based on main function signature
    idx = 0
    for model,subset in zip([Pmodel,Bmodel,PathsModel],[Psubset,Bsubset,PathsSubset]): 
        if model is not None:
            for idx,param in enumerate(signature):
                if param in model._parameter_list(order='vector'):
                    subset[getattr(model,param).idx] = idx

    #------------------------------------------------------------------------
    def Vnonlinear_fcn(*nonlin):
        """ Non-linear part of the dipolar signal function """
        # Make input arguments as array to access subsets easily
        nonlin = np.atleast_1d(nonlin)
        # Construct the basis function of the intermolecular contribution
        if Bmodel is None: 
            Bfcn = np.ones_like(t)
        else:
            Bfcn = lambda t,lam: Bmodel.nonlinmodel(t,*np.concatenate([nonlin[Bsubset],[lam]]))
        # Construct the definition of the dipolar pathways
        pathways = PathsModel.nonlinmodel(*nonlin[PathsSubset])
        # Construct the dipolar kernel
        Kdipolar = dipolarkernel(t,r,pathways=pathways, bg=Bfcn)
        # Compute the non-linear part of the distance distribution
        if Pmodel is None:        
            Pnonlin = np.eye(r)
        else: 
            Pnonlin = Pmodel.nonlinmodel(r,*nonlin[Psubset])
        # Forward calculation of the non-linear part of the dipolar signal
        Vnonlin = Kdipolar@Pnonlin
        return Vnonlin
    #------------------------------------------------------------------------

    # Create the dipolar model object
    DipolarSignal =Model(Vnonlinear_fcn,signature=signature)

    # Add the linear parameters from the subset models            
    for lparam in linearparam:
        DipolarSignal.addlinear(lparam['name'],vec=lparam['vec'])

    if Pmodel is None:
        DipolarSignal.addlinear()

    # If there are no linear paramters, add a linear scaling parameter 
    if DipolarSignal.Nlin==0:
        DipolarSignal.addlinear('scale',lb=0,par0=1)

    # Import all parameter information from the subset models
    for name,param in zip(DipolarSignal._parameter_list(order='vector'),parameters):
        getattr(DipolarSignal,name).set(**_importparameter(param))

    # Set other dipolar model specific attributes
    DipolarSignal.Pmodel = Pmodel
    DipolarSignal.Bmodel = Pmodel
    DipolarSignal.Npathways = npathways

    return DipolarSignal
#===============================================================================
