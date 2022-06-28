# dipolarmodel.py - DeerLab's dipolar EPR model generator
# ---------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2022: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
from deerlab.dipolarkernel import dipolarkernel  
from deerlab.regoperator import regoperator
from deerlab.dd_models import freedist
from deerlab.model import Model,Penalty
from deerlab import bg_hom3d
from deerlab.constants import *

#===============================================================================
def dipolarmodel(t, r, Pmodel=None, Bmodel=bg_hom3d, npathways=1, harmonics=None, experiment=None,
                    excbandwidth=np.inf, orisel=None, g=[ge,ge]):
    """
    Generate a dipolar EPR signal model.

    Parameters
    ----------
    t : array_like
        Vector of dipolar evolution times, in microseconds.
    r : array_like
        Vector of spin-spin distances, in nanometers.
    Pmodel : :ref:`Model`, optional
        Model for the distance distribution. If not specified, a non-parametric
        distance distribution defined over ``r`` is used.
    Bmodel : :ref:`Model`, optional
        Model for the intermolecular (background) contribution. If not specified,
        a background arising from a homogenous 3D distribution of spins is used.
    npathways : integer scalar, optional
        Number of dipolar pathways. If not specified, a single dipolar pathway is used.
    experiment : :ref:`ExperimentInfo`, optional
        Experimental information obtained from experiment models (``ex_``). If specified, the 
        boundaries and start values of the dipolar pathways' refocusing times and amplitudes 
        are refined based on the specific experiment's delays.
    harmonics : list of integers, optional
        Harmonics of the dipolar pathways. Must be a list with ``npathways`` harmonics, one
        for each pathway.
    orisel : callable  or ``None``, optional
        Probability distribution of possible orientations of the interspin vector to account
        for orientation selection. Must be a function taking a value of the angle θ ∈ [0,π/2]
        between the interspin vector and the external magnetic field and returning the corresponding
        probability density. If specified as ``None`` (default), a uniform distribution is used.
    excbandwidth : scalar, optional
        Excitation bandwidth of the pulses in MHz to account for limited excitation bandwidth.
        If not specified, an infinite excitation bandwidth is used.
    g : scalar, 2-element array, optional
        Electron g-values of the two spins, ``[g1, g2]``. If a single g is specified, ``[g, g]`` is used.
        If not specified, g = 2.002319... is used for both spins.

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
        raise ValueError('The number of pathways must be a positive integer.')
    if isinstance(harmonics,float) or isinstance(harmonics,int): 
        harmonics = [int(harmonics)]
    if not isinstance(harmonics,list) and harmonics is not None:
        raise TypeError('The harmonics must be specified as a list of integer values.')
    if experiment is not None:
        if not isinstance(experiment,ExperimentInfo): 
            raise TypeError('The experiment must be a valid deerlab.ExperimentInfo object.')
        # Check that the number of requested pathways does not exceed the theoretical limit of the experiment
        npathways = experiment.npathways
        maxpathways = len(experiment.reftimes)
        if npathways>maxpathways:
            raise ValueError(f'The {experiment.name} experiment can only have up to {maxpathways} dipolar pathways.')

    #------------------------------------------------------------------------
    def _importparameter(parameter):
        """ Private function for importing a parameter's metadata """
        return {
            'lb' : parameter.lb,
            'ub' : parameter.ub,
            'par0' : parameter.par0,
            'description' : parameter.description,
            'unit' : parameter.unit,
            'linear' : parameter.linear
        }
    #------------------------------------------------------------------------

    # Parse the harmonics of the dipolar pathways
    if harmonics is None: 
        if experiment is not None:
            harmonics = np.array(experiment.harmonics)[:npathways]
        else:
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

    Pnonparametric = Pmodel is None
    if Pnonparametric:
        Pmodel = freedist(r)
    Nconstants = len(Pmodel._constantsInfo)

    # Populate the basic information on the dipolar pathways parameters
    if npathways==1:
        # Special case: use modulation depth notation instead of general pathway amplitude
        getattr(PathsModel,f'mod').set(lb=0,ub=1,par0=0.01,description=f'Modulation depth',unit='')
        getattr(PathsModel,f'reftime').set(par0=0,description=f'Refocusing time',unit='μs')
    else:
        # General case: use pathway ampltiudes and refocusing times
        for n in range(npathways):
            getattr(PathsModel,f'lam{n+1}').set(lb=0,ub=1,par0=0.01,description=f'Amplitude of pathway #{n+1}',unit='')
            getattr(PathsModel,f'reftime{n+1}').set(par0=0,lb=-20,ub=20,description=f'Refocusing time of pathway #{n+1}',unit='μs')

    # Construct the signature of the dipolar signal model function
    signature = []
    parameters,linearparam = [],[]
    for model in [PathsModel,Bmodel,Pmodel]:
        if model is not None:
            for param in model._parameter_list(order='vector'):
                if np.any(getattr(model,param).linear):
                    parameters.append(getattr(model,param))
                    linearparam.append({'name':param,'vec':len(np.atleast_1d(getattr(model,param).idx)),'normalization':getattr(model,param).normalization})
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

    kernelmethod = 'fresnel' if orisel is None and np.isinf(excbandwidth) else 'grid' 

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
        Kdipolar = dipolarkernel(t,r,pathways=pathways, bg=Bfcn,
                                 excbandwidth=excbandwidth, orisel=orisel,
                                  g=g, method=kernelmethod)
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
        DipolarSignal.addlinear(lparam['name'],vec=lparam['vec'],normalization=lparam['normalization'])

    if Pmodel is None:
        DipolarSignal.addlinear()

    # If there are no linear paramters, add a linear scaling parameter 
    if DipolarSignal.Nlin==0:
        DipolarSignal.addlinear('scale',lb=0,par0=1)

    # Import all parameter information from the subset models
    for name,param in zip(DipolarSignal._parameter_list(order='vector'),parameters):
        getattr(DipolarSignal,name).set(**_importparameter(param))

    # Set prior knowledge on the parameters if experiment is specified
    if experiment is not None:

        # Compile the parameter names to change in the model
        if npathways>1:
            reftime_names = [f'reftime{n+1}' for n in range(npathways)]
        else: 
            reftime_names = ['reftime']

        # Specify start values and boundaries according to experimental timings
        for n in range(npathways):
            getattr(DipolarSignal,reftime_names[n]).set(
                par0 = experiment.reftimes[n],
                lb = experiment.reftimes[n] - 0.05,
                ub = experiment.reftimes[n] + 0.05
                )

    # Set other dipolar model specific attributes
    DipolarSignal.description = 'Dipolar signal model'

    return DipolarSignal
#===============================================================================

#===============================================================================
def dipolarpenalty(Pmodel, r, type, selection=None):
    r"""
    Construct penalties based on the distance distribution.

    Parameters
    ----------
    Pmodel : ``Model`` object 
        The Pmodel of the distance distribution, where ``None`` represents a non-parametric distance distribution. 
    r : array_like 
        Distance axis vector, in nanometers. 
    type : string
        Type of property to be imposed by the penalty. 
        
        - ``'smoothness'`` : Smoothness of the distance distribution
        - ``'compactness'`` : Compactness of the distance distribution   
        
    selection : string, optional
        Selection functional for the outer optimization of the penalty weight.

            - ``'aic'`` - Akaike information criterion
            - ``'bic'`` - Bayesian information criterion
            - ``'aicc'`` - COrrected Akaike information criterion
            - ``'icc'`` - Informational complexity criterion 
   
    Returns
    -------
    penalty : ``Penalty`` object 
        Penalty object to be passed to the ``fit`` function.
    """
    if Pmodel is None: 
        Pmodel = freedist(r)
    Nconstants = len(Pmodel._constantsInfo)

    # If include compactness penalty
    if type=='compactness':

        if selection is None: 
            selection = 'icc'

        # Define the compactness penalty function
        def compactness_penalty(*args): 
            P = Pmodel(*[r]*Nconstants,*args)
            if not np.all(P==0):
                P = P/np.trapz(P,r)
            return np.sqrt(P*(r - np.trapz(P*r,r))**2*np.mean(np.diff(r)))
        # Add the penalty to the Pmodel
        penalty = Penalty(compactness_penalty,selection,
                    signature = Pmodel._parameter_list(),
                    description = 'Distance distribution compactness penalty.')
        penalty.weight.set(lb=1e-6, ub=1e1)

    # If include smoothness penalty
    elif type=='smoothness':
        if selection is None: 
            selection = 'aic'

        # Define the smoothness penalty function
        L = regoperator(r,2)
        def smoothness_penalty(*args): 
            return L@Pmodel(*[r]*Nconstants,*args)
        # Add the penalty to the Pmodel
        penalty = Penalty(smoothness_penalty,selection,
                    signature = Pmodel._parameter_list(),
                    description = 'Distance distribution smoothness penalty.')
        penalty.weight.set(lb=1e-9, ub=1e3)

    else:
        raise KeyError(f"The requested {type} is not a valid penalty. Must be 'compactness' or 'smoothness'.")

    return penalty
#===============================================================================


#===============================================================================
def _checkpathways(pathways,Nmax):    
    # Check that pathways are correctly specified
    if len(pathways)>Nmax: 
        raise ValueError(f"The number of pathways cannot exceed {Nmax}.")
    if np.any(np.array(pathways)<1) or not np.all([not p%1 for p in pathways]): 
        raise ValueError(f"The pathways must be specified by integer numbers in the range 1-{Nmax}.")
#===============================================================================



#===============================================================================
class ExperimentInfo():
    r"""
    Represents information about a dipolar EPR experiment"""

    def __init__(self,name,reftimes,harmonics):
        self.npathways = len(reftimes)
        self.reftimes = reftimes
        self.harmonics = harmonics
        self.name = name
#===============================================================================

#===============================================================================
def ex_3pdeer(tau, pathways=None):
    r"""
    Generate a 3-pulse DEER dipolar experiment model. 


    .. image:: ../images/sequence_3pdeer.svg
        :width: 450px


    This experiment model has the following modulated dipolar pathways. The 
    theoretical refocusing times of the individual pathways are calculated 
    from the input pulse sequence delays:    

    ========= ================== ===========
     Pathway    Refocusing time    Harmonic
    ========= ================== ===========
        1             0              1
        2             𝜏              1
    ========= ================== ===========

    Parameters 
    ----------

    tau : float scalar
        Static interpulse delay. 

    pathways :  array_like, optional 
        Pathways to include in the model. The pathways are specified based to the pathways numbers. 
        By default, both pathways are included in the order given in the table above.

    Returns
    -------
    experiment : ``ExperimentInfo`` object
        Experiment object. Can be passed to ``dipolarmodel`` to introduce better
        constraints into the model.

    """

    # Theoretical refocusing times
    reftimes = [ tau, 0]
    # Theoretical dipolar harmonics
    harmonics = [ 1, 1]

    # Sort according to pathways order
    if pathways is not None:
        _checkpathways(pathways,Nmax=len(reftimes))
        reftimes = [reftimes[pathway-1] for pathway in pathways]
        harmonics = [harmonics[pathway-1] for pathway in pathways] 

    return ExperimentInfo('3-pulse DEER',reftimes,harmonics)
#===============================================================================

#===============================================================================
def ex_4pdeer(tau1, tau2, pathways=None):
    r"""
    Generate a 4-pulse DEER dipolar experiment model. 
    
    .. image:: ../images/sequence_4pdeer.svg
        :width: 450px


    This experiment model has the following modulated dipolar pathways. The 
    theoretical refocusing times of the individual pathways are calculated 
    from the input pulse sequence delays:

    ========= ================== ===========
     Pathway    Refocusing time    Harmonic
    ========= ================== ===========
        1            𝜏1              1
        2          𝜏1 + 𝜏2           1
        3             0              1
        4            𝜏2              1
    ========= ================== ===========

    Parameters 
    ----------

    tau1 : float scalar
        1st static interpulse delay. 
   
    tau2 : float scalar
        2nd static interpulse delay. 

    pathways :  array_like, optional 
        Pathways to include in the model. The pathways are specified based to the pathways numbers. 
        By default, all 4 pathways are included in the order given in the table above.

    Returns
    -------
    experiment : ``ExperimentInfo`` object
        Experiment object. Can be passed to ``dipolarmodel`` to introduce better
        constraints into the model.

    """
    # Theoretical refocusing pathways
    reftimes = [tau1, tau1+tau2, 0, tau2]
    # Theoretical dipolar harmonics
    harmonics = [1, 1, 1, 1]

    # Sort according to pathways order
    if pathways is not None:
        _checkpathways(pathways,Nmax=len(reftimes))
        reftimes = [reftimes[pathway-1] for pathway in pathways]
        harmonics = [harmonics[pathway-1] for pathway in pathways] 

    return ExperimentInfo('4-pulse DEER', reftimes, harmonics)
#===============================================================================

#===============================================================================
def ex_rev5pdeer(tau1, tau2, tau3, pathways=None):
    r"""
    Generate a reverse 5-pulse DEER dipolar experiment model. 
    
    .. image:: ../images/sequence_5pdeer_reverse.svg
        :width: 450px

    This experiment model has the following modulated dipolar pathways. The 
    theoretical refocusing times of the individual pathways are calculated 
    from the input pulse sequence delays:

    ========= ================== ===========
     Pathway    Refocusing time    Harmonic
    ========= ================== ===========
        1            𝜏3              1
        2            𝜏2              1
        3           𝜏2-𝜏3            1
        4           𝜏1+𝜏3            1
        5          𝜏1+𝜏2-𝜏3          1
        6             0              1
        7           𝜏1+𝜏2            1
        8            𝜏1              1
    ========= ================== ===========

    Parameters 
    ----------

    tau1 : float scalar
        1st static interpulse delay. 
   
    tau2 : float scalar
        2nd static interpulse delay. 

    tau3 : float scalar
        3rd static interpulse delay. 

    pathways :  array_like, optional 
        Pathways to include in the model. The pathways are specified based to the pathways numbers. 
        By default, all 8 pathways are included in the order given in the table above.

    Returns
    -------
    experiment : ``ExperimentInfo`` object
        Experiment object. Can be passed to ``dipolarmodel`` to introduce better
        constraints into the model.

    """
    # Theoretical refocusing pathways
    reftimes = [ tau3, tau2, tau2-tau3, tau1+tau3, tau1+tau2-tau3, 0, tau1+tau2, tau1]
    # Theoretical dipolar harmonics
    harmonics = [ 1, 1, 1, 1, 1, 1, 1, 1]
    
    # Sort according to pathways order
    if pathways is not None:
        _checkpathways(pathways,Nmax=len(reftimes))
        reftimes = [reftimes[pathway-1] for pathway in pathways]
        harmonics = [harmonics[pathway-1] for pathway in pathways]
    
    return ExperimentInfo('Reverse 5-pulse DEER',reftimes,harmonics)
#===============================================================================


#===============================================================================
def ex_fwd5pdeer(tau1, tau2, tau3, pathways=None):
    r"""
    Generate a forward 5-pulse DEER dipolar experiment model. 
    
    .. image:: ../images/sequence_5pdeer_forward.svg
        :width: 450px

    This experiment model has the following modulated dipolar pathways. The 
    theoretical refocusing times of the individual pathways are calculated 
    from the input pulse sequence delays:    

    ========= ================== ===========
     Pathway    Refocusing time    Harmonic
    ========= ================== ===========
        1            𝜏3              1
        2            𝜏1              1
        3           𝜏1-𝜏3            1
        4           𝜏2+𝜏3            1
        5          𝜏1+𝜏2-𝜏3          1
        6             0              1
        7           𝜏1+𝜏2            1
        8            𝜏2              1
    ========= ================== ===========

    Parameters 
    ----------

    tau1 : float scalar
        1st static interpulse delay. 
   
    tau2 : float scalar
        2nd static interpulse delay. 

    tau3 : float scalar
        3rd static interpulse delay. 

    pathways :  array_like, optional 
        Pathways to include in the model. The pathways are specified based to the pathways numbers. 
        By default, all pathways are included in the order given in the table above.

    Returns
    -------
    experiment : ``ExperimentInfo`` object
        Experiment object. Can be passed to ``dipolarmodel`` to introduce better
        constraints into the model.

    """
    # Theoretical refocusing pathways
    reftimes = [tau3, tau1, tau1-tau3, tau2+tau3, tau1+tau2-tau3, 0, tau1+tau2, tau2]
    # Theoretical dipolar harmonics
    harmonics = [1, 1, 1, 1, 1, 1, 1, 1]

    # Sort according to pathways order
    if pathways is not None:
        _checkpathways(pathways,Nmax=len(reftimes))
        reftimes = [reftimes[pathway-1] for pathway in pathways]
        harmonics = [harmonics[pathway-1] for pathway in pathways] 

    return ExperimentInfo('Forward 5-pulse DEER', reftimes, harmonics)
#===============================================================================

#===============================================================================
def ex_sifter(tau1, tau2, pathways=None):
    r"""
    Generate a SIFTER dipolar experiment model. 

    .. image:: ../images/sequence_sifter.svg
        :width: 450px

    This experiment model has the following modulated dipolar pathways. The 
    theoretical refocusing times of the individual pathways are calculated 
    from the input pulse sequence delays:    

    ========= ================== ===========
     Pathway    Refocusing time    Harmonic
    ========= ================== ===========
        1            𝜏2-𝜏1            1
        2             2𝜏2            1/2
        3            -2𝜏1            1/2
    ========= ================== ===========

    Parameters 
    ----------

    tau1 : float scalar
        1st static interpulse delay. 

    tau2 : float scalar
        2nd static interpulse delay. 

    pathways :  array_like, optional 
        Pathways to include in the model. The pathways are specified based to the pathways numbers. 
        By default, all 3 pathways are included and ordered as given in the table above.

    Returns
    -------
    experiment : ``ExperimentInfo`` object
        Experiment object. Can be passed to ``dipolarmodel`` to introduce better
        constraints into the model.

    """

    # Theoretical refocusing times
    reftimes = [ tau2-tau1, 2*tau2, -2*tau1]
    # Theoretical dipolar harmonics
    harmonics = [ 1, 1/2, 1/2]

    # Sort according to pathways order
    if pathways is not None:
        _checkpathways(pathways,Nmax=len(reftimes))
        reftimes = [reftimes[pathway-1] for pathway in pathways]
        harmonics = [harmonics[pathway-1] for pathway in pathways]

    return ExperimentInfo('SIFTER',reftimes,harmonics)
#===============================================================================


#===============================================================================
def ex_ridme(tau1, tau2, pathways=None):
    r"""
    Generate a RIDME dipolar experiment model. 

    .. image:: ../images/sequence_ridme.svg
        :width: 450px

    This experiment model has the following modulated dipolar pathways. The 
    theoretical refocusing times of the individual pathways are calculated 
    from the input pulse sequence delays:

    ========= ================== ===========
     Pathway    Refocusing time    Harmonic
    ========= ================== ===========
        1             0               1
        2             𝜏2              1
        3            -𝜏1              1
        4            𝜏2-𝜏1            1
    ========= ================== ===========

    Parameters 
    ----------

    tau1 : float scalar
        1st static interpulse delay. 

    tau2 : float scalar
        2nd static interpulse delay.  

    pathways :  array_like, optional 
        Pathways to include in the model. The pathways are specified based to the pathways numbers. 
        By default, all 4 pathways are included in the order given in the table above.

    Returns
    -------
    experiment : ``ExperimentInfo`` object
        Experiment object. Can be passed to ``dipolarmodel`` to introduce better
        constraints into the model.

    """

    # Theoretical refocusing times
    reftimes = [ 0, tau2, -tau1, tau2-tau1]
    # Theoretical dipolar harmonics
    harmonics = [ 1, 1, 1, 1]

    # Sort according to pathways order
    if pathways is not None:
        _checkpathways(pathways,Nmax=len(reftimes))
        reftimes = [reftimes[pathway-1] for pathway in pathways]
        harmonics = [harmonics[pathway-1] for pathway in pathways] 

    return ExperimentInfo('RIDME',reftimes,harmonics)
#===============================================================================
