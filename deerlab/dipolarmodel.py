# dipolarmodel.py - DeerLab's dipolar EPR model generator
# ---------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2022: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
from deerlab.dipolarkernel import dipolarkernel, dipolarbackground
from deerlab.regoperator import regoperator
from deerlab.dd_models import freedist
from deerlab.model import Model,Penalty, link
from deerlab import bg_hom3d
from deerlab.utils import choleskycovmat
from deerlab.constants import *
from math import factorial
from itertools import permutations
from scipy.stats import multivariate_normal

#===============================================================================
def dipolarmodel(t, r=None, Pmodel=None, Bmodel=bg_hom3d, npathways=1,  spins=2, harmonics=None, experiment=None,
                    excbandwidth=np.inf, orisel=None, g=[ge,ge], gridsize=1000, minamp=1e-3, specperm=True, triangles=None):
    """
    Generate a dipolar EPR signal model.

    Parameters
    ----------
    t : array_like
        Vector of dipolar evolution times, in microseconds.
    r : array_like, optional
        Vector of spin-spin distances, in nanometers. If not specified, it defaults to the range 
        1.5nm-10nm with a resolution of 0.05nm. Has no effect on models with ``spins>2``.
    Pmodel : :ref:`Model`, optional
        Model for the distance distribution. If not specified, a non-parametric
        distance distribution defined over ``r`` is used. For models with ``spins>2`` a multivariate 
        Gaussian distance distribution is assumed.
    Bmodel : :ref:`Model`, optional
        Model for the intermolecular (background) contribution. If not specified,
        a background arising from a homogenous 3D distribution of spins is used.
    npathways : integer scalar, optional
        Number of dipolar pathways. If not specified, a single dipolar pathway is used.
    harmonics : list of integers, optional
        Harmonic prefactors of the dipolar pathways. Must be a list with ``npathways`` harmonics, one
        for each pathway.
    experiment : :ref:`ExperimentInfo`, optional
        Experimental information obtained from experiment models (``ex_``). If specified, the 
        boundaries and start values of the dipolar pathways' refocusing times and amplitudes 
        are refined based on the specific experiment's delays.
    spins : integer scalar, optional 
        Number of spins in system. If not specified it defaults to two-spin systems. For multi-spins
        systems indicated by ``spins>2`` three-spin interactions are accounted for by the model automatically. 
    triangles : list of lists, optional 
        Only required when ``spins>3``. List of ``(Nspins)!/(3!*(Nspins-3)!)`` triangles formed by the different interspin distances. 
        Each triangle is specified by a list of indices of the interspin 
        distances forming the triangle, e.g. for a three-spin system ``triangles=[[1,2,3]]`` where the first, 
        second, and third interspin distances connect to form a triangle.  
    specperm : boolean, optional 
        Only when ``spins>3``. Enables the assumption of spectral permutability, i.e. the assumption that all spins in the system
        have the same spectral distribution. If enabled, all pairwise pathways of the different dipolar interactions are assumed to 
        have the same amplitude. If disabled, all the individual pathway amplitudes are parametrized separately. Enabled by default.   
    orisel : callable  or ``None``, optional
        Probability distribution of possible orientations of the interspin vector to account
        for orientation selection. Must be a function taking a value of the angle θ ∈ [0,π/2]
        between the interspin vector and the external magnetic field and returning the corresponding
        probability density. If specified as ``None`` (default), a uniform distribution is used.
    gridsize : scalar, optional
        Number of points on the grid of powder orientations to be used in the ``'grid'`` kernel 
        calculation method when evaluating models with orientation selection or multi-spin effects.
        By default set to 1000 points.
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

    # If distance axis not specified default to a long range one
    if r is None: 
        r = np.arange(1.5,10,0.05)

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

    # Determine number of dipolar interactions in the spins system
    Q = int(spins*(spins-1)/2)

    if spins>2:
        triangles_required = factorial(spins)/(factorial(3)*factorial(spins-3))
        if spins==3 and triangles is None: 
            triangles = [[1,2,3]]
        elif spins>3 and triangles is None: 
            raise SyntaxError(f'For systems with {spins} spins, a list of {triangles_required:n} triangles is required.')        
        elif spins>3 and len(triangles)!=triangles_required:
            raise SyntaxError(f'For systems with {spins} spins, a list of {triangles_required:n} triangles is required.') 
        if np.min(triangles)!=0: 
            triangles = [np.array(triangle)-1 for triangle in triangles]

    #------------------------------------------------------------------------
    def dipolarpathways(*param):
        """ 
        Generator of dipolar pathways based on a parameter set. 

        Constructs all possible two-spin dipolar pathways based on a parameter set, 
        and then constructs all possible multi-spin dipolar pathways describing three-spin
        interactions in systems with more than two spins based on the pair-wise dipolar pathways.  

        The true signature of this function is defined dynamically after its definition below.
        """
        param = np.atleast_1d(param)

        if spins==2:
            # Construct pair-wise dipolar pathways given by the parameter set
            λs = param[np.arange(0,len(param),2)]
            trefs = param[np.arange(1,len(param),2)]

            # Unmodulated pathway
            Λ0 = np.maximum(0,1 - np.sum(λs)) 
            pairpathways = [{'amp':Λ0}]

            # Modulated pathways
            for λ,tref,δ in zip(λs,trefs,harmonics):
                pairpathways.append({'amp': λ, 'reftime': tref, 'harmonic': δ})    

            # For two-spin systems these are all the required pathways
            return pairpathways

        # Construct multispin dipolar pathways given by the parameter set
        λs = [param[np.arange(q,len(param),Q)] for q in range(Q)]
        trefs = param[np.arange(Q,len(param),Q)]

        # Amplitude of the unmodulated pair-wise dipolar pathway
        λu = param[-1]
        # Initialize containers
        pathways,threespin_pathways = [],[]
        Λ0 = 1
        # Construct all possible multi-spin dipolar pathways
        for idx in range(npathways):
            # Construct all the permutations of the two-spin interaction pathways (without repetitions)
            for perm in set(set(permutations([0]+[None]*(Q-1)))):
                q = int(np.where(np.array(perm)==0)[0])
                # Compute amplitude of the two-spin interaction pathway 
                λp = λs[q][idx]
                λ2k = factorial(Q-1)*λp*λu**(Q-1)
                tref = [trefs[idx] if n==0 else None for n in perm]
                δs = [harmonics[idx] if n==0 else 0 for n in perm]
                # Add two-spin interaction dipolar pathway to the list
                pathways.append({'amp':λ2k, 'reftime': tuple(tref), 'harmonic': tuple(δs)})
                # Update the unmodulated amplitude
                Λ0 -= λ2k
                # Consider now three-spin interactions                
                for idx2 in np.arange(idx + 1,npathways):
                    # If the pathway has negligible amplitude, ignore it (speed-up) 
                    if λ3k>minamp:
                        # Construct all the permutations of the multi-spin pathway (without repetitions) 
                        for perm in set(set(permutations([0,1]+[None]*(Q-2)))):
                            # Check if the permutation corresponds to a valid interaction triangle
                            q1,q2 = [np.where(perm==n)[0] for n in [0,1]]
                            for triangle in triangles: 
                                if q1 not in triangle and q2 not in triangle:
                                    continue # If not valid, move on to the next one
                            # Compute amplitude of the three-spin interaction pathway 
                            λm = λs[q2][idx]
                            λ3k = 2*factorial(Q-2)*λp*λm*λu**(Q-2)
                            tref = [trefs[idx] if n==0 else trefs[idx2] if n==1 else None for n in perm]
                            δs = [harmonics[idx] if n==0 else harmonics[idx2] if n==1 else 0 for n in perm]
                            # Add three-spin interaction dipolar pathway to the list
                            threespin_pathways.append({'amp':λ3k, 'reftime': tuple(tref), 'harmonic': tuple(δs)})
                            # Update the unmodulated amplitude
                            Λ0 -= λ3k
        # Add pathway with unmodulated contribution
        pathways.append({'amp':Λ0})

        return pathways,threespin_pathways
    #------------------------------------------------------------------------

    # Construct the signature of the dipolarpathways() function dynamically
    if npathways == 1 and spins==2:
        # Use single-pathway notation, with amplitude as modulation depth
        variables = ['mod','reftime']
    else:
        # Use general notation
        variables = []
        for n in range(npathways):
            if spins==2:
                variables.append(f'lam{n+1}')
            else:
                for q in range(Q):
                    variables.append(f'lam{n+1}_{q+1}')
            variables.append(f'reftime{n+1}')
    # Add unmodulated pairwise pathway amplitude for multi-spin systems
    if spins>2:
        variables.append('lamu')

    # Create the dipolar pathways model object
    PathsModel = Model(dipolarpathways,signature=variables)

    if spins>2 and specperm:
        links = {f'lam{n+1}': [f'lam{n+1}_{q+1}' for q in range(Q)] for n in range(npathways)}
        PathsModel = link(PathsModel, **links)

    #-----------------------------------------------------------------------
    def _Pmultivar(param):
        rmeans,cholesky_factors = param[:Q], param[Q:]
        np.random.seed(seed=1)
        Σ = choleskycovmat(Q,cholesky_factors)
        Pmultivar = multivariate_normal(rmeans, cov=Σ)
        return Pmultivar
    #-----------------------------------------------------------------------

    # Determine if distance distributions is non-parametric
    Pnonparametric = Pmodel is None
    if Pnonparametric:
        Pmodel = freedist(r)
    # For multi-spin systems, there is no longer freedom of choice for distance distributions
    if spins>2:
        # Prepare the normal multivariate distance distribution model
        variables = []
        variables += [f'rmean{q+1}' for q in range(Q)]
        variables += [f'chol{q+1}{q+1}' for q in range(Q)]
        for q in range(Q):
            variables += [f'chol{j+1}{q+1}' for j in np.arange(q+1,Q)]
        Pmodel = Model(_Pmultivar,signature=variables)
        [getattr(Pmodel,f'rmean{q+1}').set(lb=0,ub=10,par0=3,unit='nm',description=f'Average inter-spin distance  #{q+1}') for q in range(Q)]
        [getattr(Pmodel,f'chol{q+1}{q+1}').set(lb=0,ub=6,par0=0.4,unit='nm',description=f'Cholesky factor ℓ{q+1}{q+1}') for q in range(Q)]
        for q in range(Q):
            [getattr(Pmodel,f'chol{j+1}{q+1}').set(lb=-1,ub=1,par0=0.0,unit='nm',description=f'Cholesky factor ℓ{j+1}{q+1}') for j in np.arange(q+1,Q)]
    Nconstants = len(Pmodel._constantsInfo)

    # Populate the basic information on the dipolar pathways parameters
    if spins==2 and npathways==1:
        # Special case: use modulation depth notation instead of general pathway amplitude
        getattr(PathsModel,f'mod').set(lb=0,ub=1,par0=0.01,description=f'Modulation depth',unit='')
        getattr(PathsModel,f'reftime').set(par0=0,description=f'Refocusing time',unit='μs')
    else:
        # General case: use pathway ampltiudes and refocusing times
        for n in range(npathways):
            pairwise = ''
            if spins>2: 
                pairwise = ' pairwise'   
            if spins==2 or specperm:
                getattr(PathsModel,f'lam{n+1}').set(lb=0,ub=1,par0=0.01,description=f'Amplitude of{pairwise} pathway #{n+1}',unit='')
            else:
                for q in range(Q):
                    getattr(PathsModel,f'lam{n+1}_{q+1}').set(lb=0,ub=1,par0=0.01,description=f'Amplitude of{pairwise} pathway #{n+1} on interaction #{q+1}',unit='')
            getattr(PathsModel,f'reftime{n+1}').set(par0=0,lb=-20,ub=20,description=f'Refocusing time of{pairwise} pathway #{n+1}',unit='μs')
    if spins>2: 
        getattr(PathsModel,f'lamu').set(lb=0,ub=1,par0=0.5,description='Amplitude of unmodulated pairwise pathway',unit='')
        
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
    def Vtwospin_nonlinear_fcn(*nonlin):
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
                                  g=g, method=kernelmethod, gridsize=gridsize)
        # Compute the non-linear part of the distance distribution
        Pnonlin = Pmodel.nonlinmodel(*[r]*Nconstants,*nonlin[Psubset])
        # Forward calculation of the non-linear part of the dipolar signal
        Vnonlin = Kdipolar@Pnonlin
        return Vnonlin
    #------------------------------------------------------------------------

    #------------------------------------------------------------------------
    def Vmultispin_nonlinear_fcn(*nonlin):
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

        # Sample from the multi-variate distance distribution
        Nsamples = 1000000
        np.random.seed(seed=1)
        rsamples = Pmodel.nonlinmodel(nonlin[Psubset]).rvs(Nsamples)
        rsamples = np.maximum(rsamples,1e-16) # Avoid values exactly at zero 

        # Enforce the generalized triangle inequiality
        for triangle in triangles: # Loop over all triangle combinations
            idx = triangle
            triangle_condition = np.full_like(rsamples,False)
            for q in range(Q): # Evaluate all triangle inequalities
                idx = np.roll(idx,-1,axis=0)
                triangle_condition[:,q] = np.sum(rsamples[:,idx[:-1]],axis=1) > rsamples[:,idx[-1]]
            triangle_condition = np.all(triangle_condition,axis=1)
            # Discard samples that do not satisfy triangle inequalities
            rsamples = rsamples[triangle_condition,:]

        # Generate all multi-spin dipolar pathways
        twospin_pathways, threespin_pathways = PathsModel.nonlinmodel(*nonlin[PathsSubset])

        # Start by accounting for the unmodulated contribution
        Λ0 = twospin_pathways.pop(-1)['amp']
        Vintra = Λ0

        # Account for all two-spin intramolecular contributions
        for q in range(Q):

            # Get the marginal univariate distance distribution
            Punimarg,bins = np.histogram(rsamples[:,q], bins=100)
            rq = (bins[:-1] + bins[1:])/2
            rs = [rq]*Q 
            if np.sum(Punimarg)>0:
                Punimarg = Punimarg/np.sum(Punimarg)

            # Get the two-spin pathways modulated by the q-th dipolar interaction
            pathways_q = [pathway for pathway in twospin_pathways if pathway['harmonic'][q]!=0]

            # Construct the corresponding dipolar kernel
            Ktwospin = dipolarkernel(t,rs, pathways=pathways_q, integralop=False,
                                        excbandwidth=excbandwidth, orisel=orisel,
                                        g=g, method=kernelmethod, gridsize=gridsize)

            # Two-spin intramolecular contribution from all pathways modulated by the q-th interaction
            Vintra += Ktwospin@Punimarg 

        # Account for all three-spin intramolecular contributions
        if threespin_pathways:
            for triangle in triangles:
                
                # Get the marginal trivariate distance distribution
                Ptrimarg,bins = np.histogramdd(rsamples[:,triangle],bins=10, density=True)
                r1,r2,r3 = [ (bins_q[:-1] + bins_q[1:])/2 for bins_q in bins ]

                # Get the subset of the most distance combinations with most distribution mass (speed-up)
                subset = Ptrimarg > np.max(Ptrimarg)/100
                Ptrimarg = Ptrimarg[subset]
                if np.sum(Ptrimarg)>0:
                    Ptrimarg = Ptrimarg/np.sum(Ptrimarg)
                r1,r2,r3 = [rq[subset] for rq in np.meshgrid(r1,r2,r3, indexing='ij')]

                # Get the three-spin pathways modulated by the triangle interaction
                pathways_triangle = [pathway for pathway in threespin_pathways if np.sum([pathway['harmonic'][q]!=0 for q in triangle])==2]

                # Construct the corresponding dipolar kernel
                Kthreespin = dipolarkernel(t,[r1,r2,r3], pathways=pathways_triangle, integralop=False, 
                                            excbandwidth=excbandwidth, orisel=orisel,
                                            g=g, method='grid', gridsize=gridsize)

                # Three-spin intramolecular contribution from all pathways modulated by the triangle interaction
                Vintra += Kthreespin@Ptrimarg 

        # Compute the multi-spin dipolar intermolecular contribution
        Vinter = dipolarbackground(t, twospin_pathways+threespin_pathways, Bfcn)

        # Construct the full dipolar signal 
        Vnonlin = Vintra*Vinter

        return Vnonlin
    #------------------------------------------------------------------------

    # Create the dipolar model object
    if spins==2:
        DipolarSignal = Model(Vtwospin_nonlinear_fcn,signature=signature)
    else:
        DipolarSignal = Model(Vmultispin_nonlinear_fcn,signature=signature)

    # Add the linear parameters from the subset models            
    for lparam in linearparam:
        DipolarSignal.addlinear(lparam['name'],vec=lparam['vec'],normalization=lparam['normalization'])

    if Pmodel is None:
        DipolarSignal.addlinear()

    # If there are no linear paramters, add a linear scaling parameter 
    if DipolarSignal.Nlin==0:
        DipolarSignal.addlinear('scale',lb=0,par0=1,description='Overall echo amplitude/scale')

    # Import all parameter information from the subset models
    for name,param in zip(DipolarSignal._parameter_list(order='vector'),parameters):
        getattr(DipolarSignal,name).set(**_importparameter(param))

    # Set prior knowledge on the parameters if experiment is specified
    if experiment is not None:

        # Compile the parameter names to change in the model
        if spins>2 or npathways>1:
            reftime_names = [f'reftime{n+1}' for n in range(npathways)]
        else: 
            reftime_names = ['reftime']

        # Specify start values and boundaries according to experimental timings
        for n in range(npathways):
            getattr(DipolarSignal,reftime_names[n]).set(
                par0 = experiment.reftimes[n],
                lb = experiment.reftimes[n] - 3*experiment.pulselength,
                ub = experiment.reftimes[n] + 3*experiment.pulselength
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

    def __init__(self,name,reftimes,harmonics,pulselength):
        self.npathways = len(reftimes)
        self.reftimes = reftimes
        self.harmonics = harmonics
        self.pulselength = pulselength
        self.name = name
#===============================================================================

#===============================================================================
def ex_3pdeer(tau, pathways=None, pulselength=0.016):
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

    pulselength : float scalar, optional 
        Length of the longest microwave pulse in the sequence in microseconds. Used to determine the uncertainty in the 
        boundaries of the pathay refocusing times.

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

    return ExperimentInfo('3-pulse DEER',reftimes,harmonics,pulselength)
#===============================================================================

#===============================================================================
def ex_4pdeer(tau1, tau2, pathways=None, pulselength=0.016):
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

    pulselength : float scalar, optional 
        Length of the longest microwave pulse in the sequence in microseconds. Used to determine the uncertainty in the 
        boundaries of the pathay refocusing times.

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

    return ExperimentInfo('4-pulse DEER', reftimes, harmonics, pulselength)
#===============================================================================

#===============================================================================
def ex_rev5pdeer(tau1, tau2, tau3, pathways=None, pulselength=0.016):
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

    pulselength : float scalar, optional 
        Length of the longest microwave pulse in the sequence in microseconds. Used to determine the uncertainty in the 
        boundaries of the pathay refocusing times.

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
    
    return ExperimentInfo('Reverse 5-pulse DEER',reftimes,harmonics,pulselength)
#===============================================================================


#===============================================================================
def ex_fwd5pdeer(tau1, tau2, tau3, pathways=None, pulselength=0.016):
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

    pulselength : float scalar, optional 
        Length of the longest microwave pulse in the sequence in microseconds. Used to determine the uncertainty in the 
        boundaries of the pathay refocusing times.

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

    return ExperimentInfo('Forward 5-pulse DEER',reftimes,harmonics,pulselength)
#===============================================================================

#===============================================================================
def ex_sifter(tau1, tau2, pathways=None, pulselength=0.016):
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

    pulselength : float scalar, optional 
        Length of the longest microwave pulse in the sequence in microseconds. Used to determine the uncertainty in the 
        boundaries of the pathay refocusing times.

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

    return ExperimentInfo('SIFTER',reftimes,harmonics,pulselength)
#===============================================================================


#===============================================================================
def ex_ridme(tau1, tau2, pathways=None, pulselength=0.016):
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

    pulselength : float scalar, optional 
        Length of the longest microwave pulse in the sequence in microseconds. Used to determine the uncertainty in the 
        boundaries of the pathay refocusing times.

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

    return ExperimentInfo('RIDME',reftimes,harmonics,pulselength)
#===============================================================================
