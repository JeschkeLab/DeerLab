# dipolarmodel.py - DeerLab's dipolar EPR model generator
# ---------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
from deerlab.dipolarkernel import dipolarkernel  
from deerlab.dd_models import freedist
from deerlab.model import Model
from deerlab import bg_hom3d

#===============================================================================
def dipolarmodel(t,r,Pmodel=None,Bmodel=bg_hom3d,npathways=1,harmonics=None):
    """
    Construct a dipolar EPR signal model. 

    Parameters
    ----------
    t : array_like 
        Vector of dipolar time increments, in microseconds.
    r : array_like 
        Vector of intraspin distances, in nanometers.
    Pmodel : :ref:`Model`, optional 
        Model for the distance distribution. If not speficied, a non-parametric
        distance distribution is assumed. 
    Bmodel : :ref:`Model`, optional 
        Model for the intermolecular (background) contribution. If not specified, 
        a background arising from a homogenous 3D distribution of spins is assumed. 
    npathways : integer scalar
        Number of dipolar pathways. If not specified, a single dipolar pathway is assumed. 
    harmonics: list of integers 
        Harmonics of the dipolar pathways. Must be a list with `npathways` harmonics for each
        defined dipolar pathway. 

    Returns
    -------
    Vmodel : :ref:`Model`
        Dipolar signal model object.
    """


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
            lams = param[np.arange(0,len(param),2)]
            reftimes = param[np.arange(1,len(param),2)]
            Lam0 = np.maximum(0,1 - np.sum(lams)) 
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
        variables = []
        for n in range(npathways):
            variables.append(f'lam{n+1}')
            variables.append(f'reftime{n+1}')

    # Create the dipolar pathways model object
    PathsModel = Model(dipolarpathways,signature=variables)

    if Pmodel is None:
        Pmodel = freedist(r)
    Nconstants = len(Pmodel._constantsInfo)

    # Populate the basic information on the dipolar pathways parameters
    if npathways==1:
        # Special case: use modulation depth notation instead of general pathway amplitude
        getattr(PathsModel,f'mod').set(lb=0,ub=1,par0=0.2,description=f'Modulation depth',units='')
        getattr(PathsModel,f'reftime').set(par0=0,description=f'Refocusing time',units='μs')
    else:
        # General case: use pathway ampltiudes and refocusing times
        for n in range(npathways):
            getattr(PathsModel,f'lam{n+1}').set(lb=0,ub=1,par0=0.2,description=f'Amplitude of pathway #{n+1}',units='')
            getattr(PathsModel,f'reftime{n+1}').set(par0=0,lb=-20,ub=20,description=f'Refocusing time of pathway #{n+1}',units='μs')

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
        elif hasattr(Bmodel,'lam'):
            Bfcn = lambda t,lam: Bmodel.nonlinmodel(t,*np.concatenate([nonlin[Bsubset],[lam]]))
        else:
            
            Bfcn = lambda t,_: Bmodel.nonlinmodel(t,*nonlin[Bsubset])
        # Construct the definition of the dipolar pathways
        pathways = PathsModel.nonlinmodel(*nonlin[PathsSubset])
        # Construct the dipolar kernel
        Kdipolar = dipolarkernel(t,r,pathways=pathways, bg=Bfcn)
        # Compute the non-linear part of the distance distribution
        Pnonlin = Pmodel.nonlinmodel(*[r]*Nconstants,*nonlin[Psubset])
        # Forward calculation of the non-linear part of the dipolar signal
        Vnonlin = Kdipolar@Pnonlin
        return Vnonlin
    #------------------------------------------------------------------------

    # Create the dipolar model object
    DipolarSignal = Model(Vnonlinear_fcn,signature=signature)

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


# -----------------------------------------------------------------------------------------
def _dipolarmodel_with_prior_information(t,r,reftimes,lams_par0,**kargs):
    # Generate the dipolar model
    Vmodel = dipolarmodel(t,r,**kargs)

    # Set prior knowledge on the parameters
    npathways = kargs['npathways']
    if npathways>1:
        for n in range(npathways):
            getattr(Vmodel,f'reftime{n+1}').set(par0=reftimes[n], lb=reftimes[n]-0.1, ub=reftimes[n]+0.1)
            getattr(Vmodel,f'lam{n+1}').set(par0=lams_par0[n])   
    else:
        getattr(Vmodel,f'reftime').set(par0=reftimes[0], lb=reftimes[0]-0.1, ub=reftimes[0]+0.1)
        getattr(Vmodel,f'mod').set(par0=lams_par0[0])   
    return Vmodel     
# -----------------------------------------------------------------------------------------


#===============================================================================
def model3pdeer(t,r,tau,npathways=1,**kargs):
    r"""
    Generate a 3-pulse DEER dipolar model. 
    
    The theoretically predicted refocusing times of its dipolar pathways are
    automatically computed from the pulse sequence delays and set into the model.  

    Parameters 
    ----------
    t : array_like 
        Vector of dipolar time increments, in microseconds.

    r : array_like 
        Vector of intraspin distances, in nanometers.

    tau : float scalar
        Static interpulse delay. 

    npathways : integer scalar, optional
        Number of dipolar pathways. If not specified, a single dipolar pathway is assumed. 

    Pmodel : :ref:`Model`, optional 
        Model for the distance distribution. If not speficied, a non-parametric
        distance distribution is assumed. 

    Bmodel : :ref:`Model`, optional 
        Model for the intermolecular (background) contribution. If not specified, 
        a background arising from a homogenous 3D distribution of spins is assumed. 

    Returns
    -------
    Vmodel : :ref:`Model`
        Dipolar signal model object.

    """
    # Check number of pathways does not exceed reality
    if npathways>2: 
        raise ValueError('A 3-pulse DEER signal can have up to two dipolar pathways.')
    # Theoretical refocusing pathways
    reftimes = [ tau, 0]
    # Initial guesses for the pathway amplitudes
    lams_par0 = [ 0.3, 0.05]


    Vmodel = _dipolarmodel_with_prior_information(t,r,reftimes,lams_par0,kargs)
    Vmodel.description = f'3-pulse DEER dipolar model ({npathways} dipolar pathways)'
    return Vmodel
#===============================================================================

#===============================================================================
def model4pdeer(t,r,tau1,tau2,npathways=1,**kargs):
    r"""
    Generate a 4-pulse DEER dipolar model. 
    
    The theoretically predicted refocusing times of its dipolar pathways are
    automatically computed from the pulse sequence delays and set into the model.  

    Parameters 
    ----------
    t : array_like 
        Vector of dipolar time increments, in microseconds.
   
    r : array_like 
        Vector of intraspin distances, in nanometers.
   
    tau1 : float scalar
        1st static interpulse delay. 
   
    tau2 : float scalar
        2nd static interpulse delay. 
   
    npathways : integer scalar, optional
        Number of dipolar pathways. If not specified, a single dipolar pathway is assumed. 
   
    Pmodel : :ref:`Model`, optional 
        Model for the distance distribution. If not speficied, a non-parametric
        distance distribution is assumed. 
   
    Bmodel : :ref:`Model`, optional 
        Model for the intermolecular (background) contribution. If not specified, 
        a background arising from a homogenous 3D distribution of spins is assumed. 

    Returns
    -------
    Vmodel : :ref:`Model`
        Dipolar signal model object.

    """

    # Check number of pathways does not exceed reality
    if npathways>4: 
        raise ValueError('A 4-pulse DEER signal can have up to four dipolar pathways.')
    # Theoretical refocusing pathways
    reftimes = [ tau1, tau1+tau2, 0, tau2 ]
    # Initial guesses for the pathway amplitudes
    lams_par0 = [ 0.3, 0.05, 0.05, 0.05]

    Vmodel = _dipolarmodel_with_prior_information(t,r,reftimes,lams_par0,npathways,**kargs)
    Vmodel.description = f'4-pulse DEER dipolar model ({npathways} dipolar pathways)'
    return Vmodel
#===============================================================================

#===============================================================================
def model5pdeer(t,r,tau1,tau2,tau3,npathways=1,**kargs):
    r"""
    Generate a 5-pulse DEER dipolar model. 
    
    The theoretically predicted refocusing times of its dipolar pathways are
    automatically computed from the pulse sequence delays and set into the model.  

    Parameters 
    ----------
  
    t : array_like 
        Vector of dipolar time increments, in microseconds.
  
    r : array_like 
        Vector of intraspin distances, in nanometers.
  
    tau1 : float scalar
        1st static interpulse delay. 
  
    tau2 : float scalar
        2nd static interpulse delay. 
  
    tau3 : float scalar
        3rd static interpulse delay. 
  
    npathways : integer scalar, optional
        Number of dipolar pathways. If not specified, a single dipolar pathway is assumed. 

    Pmodel : :ref:`Model`, optional 
        Model for the distance distribution. If not speficied, a non-parametric
        distance distribution is assumed. 
  
    Bmodel : :ref:`Model`, optional 
        Model for the intermolecular (background) contribution. If not specified, 
        a background arising from a homogenous 3D distribution of spins is assumed. 
  
    Returns
    -------
    Vmodel : :ref:`Model`
        Dipolar signal model object.

    """
    # Check number of pathways does not exceed reality
    if npathways>8: 
        raise ValueError('A 5-pulse DEER signal can have up to eight dipolar pathways.')
    # Theoretical refocusing pathways
    reftimes = [ tau3, tau2, tau2-tau3, tau1+tau3, 0, tau1+tau2, tau1+tau2-tau3, tau1]
    # Initial guesses for the pathway amplitudes
    lams_par0 = [ 0.3, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05]

    Vmodel = _dipolarmodel_with_prior_information(t,r,reftimes,lams_par0,**kargs)
    Vmodel.description = f'5-pulse DEER dipolar model ({npathways} dipolar pathways)'
    return Vmodel
#===============================================================================


