# dipolarmodel.py - DeerLab's dipolar EPR model generator
# ---------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2022: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
from deerlab.dipolarkernel import dipolarkernel
from deerlab.dipolarbackground import dipolarbackground, dipolarbackgroundmodel, _multispin_pathways
from deerlab.regoperator import regoperator
from deerlab.dd_models import freedist
from deerlab.model import Model,Penalty, link
from deerlab import bg_hom3d
from deerlab.utils import choleskycovmat
from deerlab.constants import *
from math import factorial
import inspect
from scipy.stats import multivariate_normal, norm, qmc


#===============================================================================
def _normalsobol(dim, m, seed=0):
    """
    Fixed quasi-Monte Carlo point set of the standard normal distribution in ``dim`` dimensions: ``2**m`` points
    of a scrambled Sobol sequence (with a fixed seed) mapped by the inverse normal CDF, and their (equal) weights.
    Compared to a Cartesian grid with the same number of points, it does not alias rapidly oscillating integrands.
    """
    try:
        sobol = qmc.Sobol(d=dim, scramble=True, rng=seed)
    except TypeError: # SciPy<1.15
        sobol = qmc.Sobol(d=dim, scramble=True, seed=seed)
    z = norm.ppf(sobol.random_base2(m))
    return z, np.full(len(z),1/len(z))
#===============================================================================

#===============================================================================
def _normalgrid(dim, spacing, radius):
    """
    Fixed quadrature grid of the standard normal distribution in ``dim`` dimensions: the nodes of a
    Cartesian grid with the given ``spacing`` lying within a sphere of the given ``radius``, and their
    normalized standard normal weights.
    """
    nodes1d = np.arange(-np.floor(radius/spacing),np.floor(radius/spacing)+1)*spacing
    z = np.stack(np.meshgrid(*[nodes1d]*dim,indexing='ij'),axis=-1).reshape(-1,dim)
    z = z[np.sum(z**2,axis=1)<=radius**2]
    w = np.exp(-np.sum(z**2,axis=1)/2)
    return z, w/np.sum(w)
#===============================================================================

#===============================================================================
def _triangleweight(r1, r2, r3, width=0.05):
    """
    C1-smooth indicator of the triangle inequality for the distances ``r1``, ``r2``, ``r3`` (in nm). It is exactly
    zero outside the allowed region (which also excludes non-positive distances), one if all inequalities are
    fulfilled by a margin larger than ``width``, and rises smoothly in between, so that a multivariate
    distribution truncated by it remains a smooth function of its parameters.
    """
    def smoothstep(margin):
        x = np.clip(margin/width,0,1)
        return x*x*(3 - 2*x)
    return smoothstep(r2 + r3 - r1)*smoothstep(r1 + r3 - r2)*smoothstep(r1 + r2 - r3)
#===============================================================================

#===============================================================================
def _cholesky2x2(Σ):
    """ Lower Cholesky factor of a 2x2 covariance matrix, robust to (numerically) singular matrices """
    a = np.sqrt(max(Σ[0,0],0))
    b = Σ[1,0]/a if a>0 else 0
    c = np.sqrt(max(Σ[1,1] - b**2,0))
    return np.array([[a,0],[b,c]])
#===============================================================================

# Fixed quadrature point sets of the standard normal distribution used to evaluate the multivariate normal distance
# distribution of multi-spin models. At each evaluation they are mapped onto the distance distribution via
# r = rmean + L@z, such that the model is deterministic and a smooth function of all its parameters.
_Z1 = np.linspace(-4.5,4.5,121)            # Univariate marginals (two-spin contributions)
_W1 = np.exp(-_Z1**2/2)/np.sum(np.exp(-_Z1**2/2))
_Z2,_W2 = _normalgrid(2,0.5,3.5)            # Conditional bivariate distributions (triangle inequality of univariate marginals)
_Z3,_W3 = _normalsobol(3,9)                 # Trivariate marginals (three-spin contributions)

#===============================================================================
def _modulated(pathway, q):
    """ Whether a multi-spin dipolar pathway is modulated by the q-th dipolar interaction """
    return 'harmonic' in pathway and np.any(np.atleast_1d(pathway['harmonic'][q])!=0)
#===============================================================================

#===============================================================================
def _univariate_marginal(rmeans, Σ, q, triangles):
    """
    Univariate marginal distribution of the q-th distance of a multivariate normal distance distribution truncated by
    the triangle inequality of all triangles containing the q-th distance (exact for three-spin systems, treating
    the triangles as independent otherwise). Returns the distances and the normalized probabilities at the nodes with
    non-zero probability.
    """
    Σqq = max(Σ[q,q],1e-12)
    rq = rmeans[q] + np.sqrt(Σqq)*_Z1
    P = _W1.copy()
    for triangle in triangles:
        if q not in triangle:
            continue
        # Probability of fulfilling the triangle inequality conditioned on the q-th distance, from the conditional
        # normal distribution of the two other distances in the triangle
        o = [qq for qq in triangle if qq!=q]
        Σoq = Σ[o,q]
        μcond = rmeans[o] + np.outer(rq - rmeans[q], Σoq/Σqq)
        Lcond = _cholesky2x2(Σ[np.ix_(o,o)] - np.outer(Σoq,Σoq)/Σqq)
        ro = μcond[:,np.newaxis,:] + (_Z2@Lcond.T)[np.newaxis,:,:]
        P = P*np.sum(_W2*_triangleweight(rq[:,np.newaxis], ro[...,0], ro[...,1]),axis=1)
    keep = P>0
    rq,P = rq[keep],P[keep]
    return rq, P/np.sum(P) if np.sum(P)>0 else P
#===============================================================================

#===============================================================================
def _trivariate_marginal(rmeans, Σ, L, triangle):
    """
    Trivariate marginal distribution of the distances forming a triangle, of a multivariate normal distance
    distribution truncated by the triangle inequality. Returns the distances (as an array of shape (N,3)) and
    the normalized probabilities at the nodes with non-zero probability.
    """
    if len(rmeans)==3:
        # The triangle contains all distances, the Cholesky factor of the distribution can be used directly
        Ltri = L[np.ix_(triangle,triangle)]
    else:
        Ltri = np.linalg.cholesky(Σ[np.ix_(triangle,triangle)] + 1e-12*np.eye(3))
    r = rmeans[triangle] + _Z3@Ltri.T
    P = _W3*_triangleweight(*r.T)
    keep = P>0
    r,P = r[keep],P[keep]
    return r, P/np.sum(P) if np.sum(P)>0 else P
#===============================================================================


#===============================================================================
def _populate_dipolar_pathways_parameters(model,npathways,spins=2,samespins=True,experiment=None,parametrization='reftimes',Q=None):
    """
    A private function to populate the dipolar pathways parameters in the model.

    Parameters
    ----------
    model : :ref:`Model`
        Model object to be populated with the dipolar pathways parameters.
    npathways : int
        Number of dipolar pathways.
    spins : int
        Number of spins in the system.
    samespins : bool    
        Enables the assumption of spectral permutability.
    experiment : :ref:`ExperimentInfo`
        Experimental information obtained from experiment models (``ex_``).
    parametrization : str
        Parametrization strategy of the dipolar pathway refocusing times.
    Q : int
        Number of interactions in the spins system.

    Returns
    -------
    model : :ref:`Model`
        Model object populated with the dipolar pathways parameters.
    """

    # Populate the basic information on the dipolar pathways parameters
    if experiment is not None:
        labels = experiment.labels
        pulsedelay_names = inspect.signature(experiment.reftimes).parameters.keys()
    else:
        # Otherwise, just label them in order
        labels = np.arange(1,npathways+1)
    
    if spins>2:
        if Q is None:
            raise ValueError('If spins>2, the number of interactions Q must be specified.')
        
    
    for n in range(npathways):
        if spins==2 and npathways==1:
            # Special case: use modulation depth notation instead of general pathway amplitude
            getattr(model,f'mod').set(lb=0,ub=1,par0=0.01,description=f'Modulation depth',unit='')
        else:
            pairwise = ''
            if spins>2: 
                pairwise = ' pairwise'   
            if spins==2 or samespins:
                getattr(model,f'lam{labels[n]}').set(lb=0,ub=1,par0=0.01,description=f'Amplitude of{pairwise} pathway #{labels[n]}',unit='')
            else:
                for q in range(Q):
                    getattr(model,f'lam{labels[n]}_{q+1}').set(lb=0,ub=1,par0=0.01,description=f'Amplitude of{pairwise} pathway #{labels[n]} on interaction #{q+1}',unit='')
        if parametrization=='reftimes':
            if experiment is None:
                theoretical_reftime = 0
                reftime_variability = 20
            else:
                theoretical_reftime = experiment.reftimes(*experiment.delays)[n]
                reftime_variability = 3*experiment.pulselength
            if spins==2 and npathways==1:
                # Special case: use reftime notation
                getattr(model,f'reftime').set(description=f'Refocusing time',unit='μs',
                    par0 = theoretical_reftime,  
                    lb   = theoretical_reftime - reftime_variability,
                    ub   = theoretical_reftime + reftime_variability)
            else: 
                # Special case: use reftime notation
                getattr(model,f'reftime{labels[n]}').set(description=f'Refocusing time of{pairwise} pathway #{labels[n]}',unit='μs',
                    par0 = theoretical_reftime,
                    lb   = theoretical_reftime - reftime_variability,
                    ub   = theoretical_reftime + reftime_variability)
    if parametrization=='delays': 
        experimental_delays = experiment.delays
        delay_variability = 5*experiment.pulselength
        for n,delay in enumerate(pulsedelay_names):
            getattr(model,delay).set(description=f'Pulse delay {delay} of the experiment',unit='μs',
                par0 = experimental_delays[n],
                lb   = experimental_delays[n] - delay_variability,
                ub   = experimental_delays[n] + delay_variability)
    elif parametrization=='shift':
        variability = 5*experiment.pulselength
        getattr(model,'tshift').set(par0=0,lb=-variability,ub=variability,description=f'Variability of experimental pulse delays',unit='μs')
    if spins>2: 
        getattr(model,f'lamu').set(lb=0,ub=1,par0=0.5,description='Amplitude of unmodulated pairwise pathway',unit='')
    
    return model


def dipolarmodel(t, r=None, Pmodel=None, Bmodel=bg_hom3d, experiment=None, parametrization='reftimes', npathways=1,  spins=2, harmonics=None,
                    excbandwidth=np.inf, orisel=None, g=[ge,ge], gridsize=1000, minamp=1e-3, samespins=True, triangles=None, interp=True,):
    """
    Generate a dipolar EPR signal model.
    
    This function generates a model for the dipolar EPR signal arising from the dipolar interactions
    between the electron spins of the system. The model is defined in terms of the evolution time vector 
    ``t``, the spin-spin distances ``r``, the distance distribution model ``Pmodel``, the intermolecular contribution
    ``Bmodel``, and an ``experiment`` object containing information about the specific experiment. 
    The number of spins in the system ``spins`` can be specified, as well as the parametrization strategy 
    ``parametrization`` for the dipolar pathway refocusing times. 
    
    For multi-spin systems (``spins > 2``), the function allows for the inclusion of three-spin interactions by
    specifying the ``triangles`` list of interspin triangles and enabling the ``samespins`` assumption of spectral
    permutability. The ``minamp`` parameter can be used to threshold the amplitude of the three-spin interactions. 


    Additional effects such as orientation selection or limited excitation bandwidth can be included via the ``orisel`` 
    and ``excbandwidth`` arguments, respectively. 

    Parameters
    ----------
    t : array_like
        Vector of dipolar evolution times, in microseconds.
    
    r : array_like, optional
        Vector of spin-spin distances, in nanometers. If not specified, it defaults to the range 
        1.5nm-8nm with a resolution of 0.05nm. Has no effect on models with ``spins>2``.
    
    Pmodel : :ref:`Model`, optional
        Model for the distance distribution. If not specified, a non-parametric
        distance distribution defined over ``r`` is used. For models with ``spins>2`` a multivariate
        Gaussian distance distribution is assumed, truncated by the triangle inequality. Its univariate and
        trivariate marginals are evaluated by deterministic quadrature on fixed point sets that scale with its
        means and covariance.
    
    Bmodel : :ref:`Model`, optional
        Model for the intermolecular (background) contribution. If not specified,
        a background arising from a homogenous 3D distribution of spins is used.
    
    experiment : :ref:`ExperimentInfo`, optional
        Experimental information obtained from experiment models (``ex_``). If specified, the 
        boundaries and start values of the dipolar pathways' refocusing times and amplitudes 
        are refined based on the specific experiment's delays.
    
    parametrization : string, optional
        Parametrization strategy of the dipolar pathway refocusing times. Must be one of the following:
            
            * ``'reftimes'`` - Each refocusing time is represented individually as a parameter.
            * ``'delays'`` - The pulse delays are introduced as parameters from which the refocusing times are computed. Requires ``experiment`` to be specified.
            * ``'shift'`` - A time shift is introduced as a parameter to represent the variability of the refocusing times from their theoretical values. Requires ``experiment`` to be specified.

        The default is ``'reftimes'``. 
    
    npathways : integer scalar, optional
        Number of dipolar pathways. If not specified, a single dipolar pathway is used.
    
    harmonics : list of integers, optional
        Harmonic prefactors of the dipolar pathways. Must be a list with ``npathways`` harmonics, one
        for each pathway.

    spins : integer scalar, optional 
        Number of spins in system. If not specified it defaults to two-spin systems. For multi-spins
        systems indicated by ``spins>2`` three-spin interactions are accounted for by the model automatically. 
    
    triangles : list of lists, optional 
        Only required when ``spins>3``. List of ``(Nspins)!/(3!*(Nspins-3)!)`` triangles formed by the different interspin distances. 
        Each triangle is specified by a list of indices of the interspin 
        distances forming the triangle, e.g. for a three-spin system ``triangles=[[1,2,3]]`` where the first, 
        second, and third interspin distances connect to form a triangle.  
    
    samespins : boolean, optional
        Only when ``spins>2``. Enables the assumption of spectral permutability, i.e. the assumption that all spins in the system
        have the same spectral distribution. If enabled, all pairwise pathways of the different dipolar interactions are assumed to 
        have the same amplitude. If disabled, all the individual pathway amplitudes are parametrized separately. Enabled by default.   
    
    minamp : float scalar, optional
        Only when ``spins>2``. Threshold amplitude for neglecting three-spin interaction pathways. To enhance performance, all dipolar
        pathways arising from three-spin interactions with an amplitude ``lam<minamp`` are approximated as having zero amplitude.
        Between ``minamp`` and ``2*minamp`` their amplitude is switched on linearly, to keep the model continuous in its parameters.
        The default threshold is ``1e-3``.
    
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
    
    interp : boolean, optional 
        Enable dipolar kernel interpolation for computation time reduction. By default enabled. 
    
    Returns
    -------
    Vmodel : :ref:`Model`
        Dipolar signal model object. With these additional attributes:
            * ``t`` - The time axis of the model.
            * ``r`` - The distance axis of the model, only for two-spin systems.
            * ``Bmodel`` - The multi-pathway background model object used in the construction of the dipolar signal model.
            * ``Pmodel`` - The distance distribution model object used in the construction of the dipolar signal model.
            

    """

    # If distance axis not specified default to a long range one
    if r is None: 
        r = np.arange(1.5,8,0.05)

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
        maxpathways = len(experiment.harmonics)
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
            'linear' : parameter.linear,
            'frozen' : parameter.frozen,
            'value' : parameter.value
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

    if parametrization!='reftimes': 
        if experiment is None:
            raise SyntaxError(f"The '{parametrization}' parametrization requires 'experiment' to be specified.")
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
        if np.ndim(triangles)!=2 or np.shape(triangles)[1]!=3:
            raise SyntaxError('The triangles must be specified as a list of lists of three interaction indices each, e.g. triangles=[[1,2,3]].')
        if np.min(triangles)!=0:
            triangles = [np.array(triangle)-1 for triangle in triangles]
        triangles = [np.array(triangle,dtype=int) for triangle in triangles]
        if np.min(triangles)<0 or np.max(triangles)>=Q:
            raise SyntaxError(f'The triangles must contain interaction indices between 1 and {Q}.')

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
            if parametrization=='reftimes':
                λs = param[np.arange(0,len(param),2)]
                trefs = param[np.arange(1,len(param),2)]
            elif parametrization=='delays':
                delays = param[0:len(experiment.delays)]
                λs = param[len(delays):]
                trefs = experiment.reftimes(*delays)
            elif parametrization=='shift':
                shift = param[0]
                λs = param[1:]
                delays = np.array(experiment.delays) + shift
                trefs = experiment.reftimes(*delays)

            # Unmodulated pathway
            Λ0 = np.maximum(0,1 - np.sum(λs)) 
            pairpathways = [{'amp':Λ0}]

            # Modulated pathways
            for λ,tref,δ in zip(λs,trefs,harmonics):
                pairpathways.append({'amp': λ, 'reftime': tref, 'harmonic': δ})    

            # For two-spin systems these are all the required pathways
            return pairpathways

        # Amplitude of the unmodulated pair-wise dipolar pathway
        λu = param[-1]
        # Parse the pair-wise pathway amplitudes λs[q][p] and refocusing times
        if parametrization=='reftimes':
            # Parameters ordered as (lam{p}_1,...,lam{p}_Q, reftime{p}) for each pathway p
            blocks = np.reshape(param[:-1],(npathways,Q+1))
            λs, trefs = blocks[:,:Q].T, blocks[:,Q]
        elif parametrization=='delays':
            delays = param[0:len(experiment.delays)]
            λs = np.reshape(param[len(delays):-1],(npathways,Q)).T
            trefs = experiment.reftimes(*delays)
        elif parametrization=='shift':
            delays = np.array(experiment.delays) + param[0]
            λs = np.reshape(param[1:-1],(npathways,Q)).T
            trefs = experiment.reftimes(*delays)

        # Construct all two-spin and three-spin dipolar pathways
        return _multispin_pathways(λs, trefs, harmonics, λu, Q, triangles=triangles, minamp=minamp)
    #------------------------------------------------------------------------

    # Construct the signature of the dipolarpathways() function dynamically
    if experiment is not None:
        labels = experiment.labels
    else:
        # Otherwise, just label them in order
        labels = np.arange(1,npathways+1)

    variables = []
    if parametrization=='delays': 
        pulsedelay_names = inspect.signature(experiment.reftimes).parameters.keys()
        variables = list(pulsedelay_names)
    if parametrization=='shift': 
        variables = ['tshift']
    for n in range(npathways):
        if npathways == 1 and spins==2:
            # Use single-pathway notation, with amplitude as modulation depth
            variables.append('mod')
        else:
            if spins==2:
                variables.append(f'lam{labels[n]}')
            else:
                for q in range(Q):
                    variables.append(f'lam{labels[n]}_{q+1}')
        if parametrization=='reftimes':
            if npathways == 1 and spins==2:
                variables.append(f'reftime')
            else: 
                variables.append(f'reftime{labels[n]}')
    # Add unmodulated pairwise pathway amplitude for multi-spin systems
    if spins>2:
        variables.append('lamu')

    # Create the dipolar pathways model object
    PathsModel = Model(dipolarpathways,signature=variables)

    if spins>2 and samespins:
        links = {f'lam{labels[n]}': [f'lam{labels[n]}_{q+1}' for q in range(Q)] for n in range(npathways)}
        PathsModel = link(PathsModel, **links)

    #-----------------------------------------------------------------------
    def _Pmultivar(param):
        rmeans,cholesky_factors = param[:Q], param[Q:]
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

    PathsModel = _populate_dipolar_pathways_parameters(
        PathsModel,npathways,spins=spins,samespins=samespins,experiment=experiment,parametrization=parametrization,Q=Q)
            
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

    dt = np.mean(np.diff(t))
    tinterp = None

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
        Kdipolar = dipolarkernel(t, r, pathways=pathways, bg=Bfcn, tinterp=tinterp,
                                 excbandwidth=excbandwidth, orisel=orisel,
                                  g=g, method=kernelmethod, gridsize=gridsize)

        # Compute the non-linear part of the distance distribution
        Pnonlin = Pmodel.nonlinmodel(*[r]*Nconstants,*nonlin[Psubset])

        # Forward calculation of the non-linear part of the dipolar signal
        Vnonlin = Kdipolar@Pnonlin
        return Vnonlin
    #------------------------------------------------------------------------

    #------------------------------------------------------------------------
    def multispin_contributions(*nonlin):
        """ 
        Two-spin (including the unmodulated contribution) and three-spin intramolecular contributions,
        and intermolecular contribution of the non-linear part of the multi-spin dipolar signal 
        """
        # Make input arguments as array to access subsets easily
        nonlin = np.atleast_1d(nonlin)
        # Construct the basis function of the intermolecular contribution
        if Bmodel is None: 
            Bfcn = np.ones_like(t)
        elif hasattr(Bmodel,'lam'):
            Bfcn = lambda t,lam: Bmodel.nonlinmodel(t,*np.concatenate([nonlin[Bsubset],[lam]]))
        else:
            Bfcn = lambda t,_: Bmodel.nonlinmodel(t,*nonlin[Bsubset])

        # Parameters of the multivariate normal distance distribution: means and lower Cholesky factor
        Pparam = nonlin[Psubset]
        rmeans = Pparam[:Q]
        L = np.zeros((Q,Q))
        L[np.diag_indices(Q)] = Pparam[Q:2*Q]
        cols,rows = np.triu_indices(Q,k=1)
        L[rows,cols] = Pparam[2*Q:]
        Σ = L@L.T

        # Generate all multi-spin dipolar pathways
        twospin_pathways, threespin_pathways = PathsModel.nonlinmodel(*nonlin[PathsSubset])

        # Start by accounting for the unmodulated contribution
        Λ0 = twospin_pathways.pop(-1)['amp']
        Vtwospin = np.full(len(t),Λ0,dtype=float)
        Vthreespin = np.zeros(len(t))

        # Account for all two-spin intramolecular contributions
        for q in range(Q):

            # Get the two-spin pathways modulated by the q-th dipolar interaction
            pathways_q = [pathway for pathway in twospin_pathways if _modulated(pathway,q)]
            if not pathways_q:
                continue

            # Get the marginal univariate distance distribution
            rq,Punimarg = _univariate_marginal(rmeans,Σ,q,triangles)

            # Construct the corresponding dipolar kernel
            Ktwospin = dipolarkernel(t,[rq]*Q, pathways=pathways_q, integralop=False, tinterp=tinterp,
                                        excbandwidth=excbandwidth, orisel=orisel,
                                        g=g, method=kernelmethod, gridsize=gridsize)

            # Two-spin intramolecular contribution from all pathways modulated by the q-th interaction
            Vtwospin += Ktwospin@Punimarg

        # Account for all three-spin intramolecular contributions
        for triangle in triangles:

            # Get the three-spin pathways modulated by the triangle interaction, restricted to the triangle's
            # three interactions and with identical pathways merged (speed-up)
            pathways_triangle = {}
            for pathway in threespin_pathways:
                if np.sum([_modulated(pathway,q) for q in triangle])!=2:
                    continue
                tref = tuple(pathway['reftime'][q] for q in triangle)
                δs = tuple(pathway['harmonic'][q] for q in triangle)
                key = repr((tref,δs))
                if key in pathways_triangle:
                    pathways_triangle[key]['amp'] += pathway['amp']
                else:
                    pathways_triangle[key] = {'amp':pathway['amp'], 'reftime':tref, 'harmonic':δs}
            if not pathways_triangle:
                continue

            # Get the marginal trivariate distance distribution
            r3,Ptrimarg = _trivariate_marginal(rmeans,Σ,L,triangle)

            # Construct the corresponding dipolar kernel
            Kthreespin = dipolarkernel(t,list(r3.T), pathways=list(pathways_triangle.values()), integralop=False,
                                        excbandwidth=excbandwidth, orisel=orisel,
                                        g=g, method='grid', gridsize=gridsize)

            # Three-spin intramolecular contribution from all pathways modulated by the triangle interaction
            Vthreespin += Kthreespin@Ptrimarg

        # Compute the multi-spin dipolar intermolecular contribution
        if Bmodel is None:
            Vinter = 1
        else:
            Vinter = dipolarbackground(t, twospin_pathways+threespin_pathways, Bfcn)

        return Vtwospin, Vthreespin, Vinter
    #------------------------------------------------------------------------

    #------------------------------------------------------------------------
    def Vmultispin_nonlinear_fcn(*nonlin):
        """ Non-linear part of the dipolar signal function """
        Vtwospin, Vthreespin, Vinter = multispin_contributions(*nonlin)

        # Construct the full dipolar signal 
        Vnonlin = (Vtwospin + Vthreespin)*Vinter

        return Vnonlin
    #------------------------------------------------------------------------

    #------------------------------------------------------------------------
    def contributions(**params):
        """
        Two-spin and three-spin contributions to the multi-spin dipolar signal, such that ``V = Vtwospin + Vthreespin``.
        Takes the same keyword arguments as the model call. The two-spin contribution includes the unmodulated
        contribution, and both contributions include the intermolecular contribution and the overall scale.
        """
        missing = [param for param in signature+['scale'] if param not in params]
        if missing:
            raise KeyError(f'Missing model parameters: {missing}')
        Vtwospin, Vthreespin, Vinter = multispin_contributions(*[params[param] for param in signature])
        return {'twospin': params['scale']*Vtwospin*Vinter, 'threespin': params['scale']*Vthreespin*Vinter}
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

    DipolarSignal.description = 'Dipolar signal model'
    # Set prior knowledge on the parameters if experiment is specified
    if experiment is not None:

        # Specify experiment in model description
        DipolarSignal.description = f'{experiment.name} dipolar signal model'

        # Specify start values and boundaries according to experimental timings
        reftimes = experiment.reftimes(*experiment.delays)
        if interp:
            tinterp = np.arange( 
                    min(t)-max(reftimes) - 3*experiment.pulselength - dt,
                    max(t)+max(reftimes) + 3*experiment.pulselength + dt,
                    dt)

    # 
    if Bmodel is not None:
        DipolarSignal.Bmodel = dipolarbackgroundmodel(experiment=experiment, basis=Bmodel, parametrization=parametrization, spins=spins, samespins=samespins, npathways=npathways,
                                                      triangles=triangles if spins>2 else None, minamp=minamp if spins>2 else 0)
    
    if Pmodel is not None:
        DipolarSignal.Pmodel = Pmodel
    DipolarSignal.t = t
    if spins==2:
        DipolarSignal.r = r
    else:
        DipolarSignal.contributions = contributions

    return DipolarSignal
#===============================================================================

#===============================================================================
def dipolarpenalty(Pmodel, r, type, selection=None):
    r"""
    Construct penalties based on the distance distribution.


    This function constructs penalties for a given distance distribution that can be used to impose
    certain desired properties on the distribution. The type of penalty and the function to optimize 
    it are determined by the ``type`` parameter. The ``selection`` parameter allows the user to choose 
    a method for optimizing the penalty weight.

    Parameters
    ----------
    Pmodel : `Model` object 
        The model of the distance distribution. If ``Pmodel`` is ``None``, a non-parametric distance distribution is assumed. 
    r : array_like 
        Distance axis vector, in nanometers. This vector specifies the distances at which the distance distribution is defined.
    type : string
        Type of property to be imposed by the penalty. This parameter determines the type of penalty 
        that will be constructed, as well as the function that will be optimized to determine the 
        penalty weight. The following values are supported for this parameter:
        
        - ``'smoothness'`` : Imposes a smoothness constraint on the distance distribution. This penalty 
          is optimized by minimizing the L2 norm of the second derivative of the distribution.
        - ``'compactness'`` : Imposes a compactness constraint on the distance distribution. This penalty 
          is optimized by minimizing the square root of the product of the distance distribution and 
          the square of the difference between the distances and their mean.
        
    selection : string, optional
        Selection functional for the outer optimization of the penalty weight. Possible values:

            - ``'aic'`` - Akaike information criterion
            - ``'bic'`` - Bayesian information criterion
            - ``'aicc'`` - COrrected Akaike information criterion
            - ``'icc'`` - Informational complexity criterion 
   
       If no value is provided, ``'icc'`` is used for ``'compactness'`` and ``'aic'`` for ``'smoothness'``.

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
                P = P/np.trapezoid(P,r)
            return np.sqrt(P*(r - np.trapezoid(P*r,r))**2*np.mean(np.diff(r)))
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
    if np.any(np.array(pathways)<1) or not np.all([not p%1 for p in pathways]) or np.any(np.array(pathways)>Nmax): 
        raise ValueError(f"The pathways must be specified by integer numbers in the range 1-{Nmax}.")
    if len(np.unique(pathways))!=len(pathways):
        raise ValueError("Some pathways are duplicated.")
#===============================================================================



#===============================================================================
class ExperimentInfo():
    r"""
    Represents information about a dipolar EPR experiment
    
    Attributes    
    ----------
    name : string
        Name of the experiment.
    reftimes : list of floats
        List of refocusing times, in microseconds.
    harmonics : list
        List of harmonics of the refocusing times.
    pulselength : float
        Length of the longest microwave pulse in the sequence in microseconds.
    pathwaylabels : list
        List of pathway labels.
    delays : list
        List of pulse delays.
            
    """

    def __init__(self,name,reftimes,harmonics,pulselength,pathwaylabels,delays):
        """ Populates the attributes of the ExperimentInfo object.

        Parameters
        ----------
        name : string
            Name of the experiment.
        reftimes : list of floats
            List of refocusing times, in microseconds.
        harmonics : list
            List of harmonics of the refocusing times.  
        pulselength : float
            Length of the longest microwave pulse in the sequence in microseconds.
        pathwaylabels : list
            List of pathway labels.
        delays : list
            List of pulse delays.
        """
        self.npathways = len(harmonics)
        self.reftimes = reftimes
        self.harmonics = harmonics
        self.pulselength = pulselength
        self.labels = pathwaylabels
        self.delays = delays
        self.name = name
#===============================================================================

#===============================================================================
def ex_3pdeer(tau, pathways=[1,2], pulselength=0.016):
    r"""
    Generate a 3-pulse DEER dipolar experiment model. 
    
    The figure below shows the dipolar pathways in 3-pulse DEER. The observer (black) and pump (grey) pulses and their interpulse delays are shown on the top.
    The middle table summarizes all detectable modulated dipolar pathways `p` along their dipolar phase accumulation factors `\mathbf{s}_p`,
    harmonics `\delta_p` and refocusing times `t_{\mathrm{ref},p}`. The most commonly encountered pathways are highlighted in color. 
    The bottom panel shows a decomposition of the dipolar signal into the individual intramolecular contributions (shown as colored lines).

    .. figure:: ../images/ex_3pdeer_pathways.png
        :width: 350px

        Source: L. Fábregas-Ibáñez, Advanced data analysis and modeling in dipolar EPR spectroscopy, Doctoral dissertation, 2022


    Parameters 
    ----------

    tau : float scalar
        Static interpulse delay `\tau`. 

    pathways :  array_like, optional 
        Pathways to include in the model. The pathways are specified as a list of pathways labels `p`. 
        By default, both pathways are included as shown in the table above.

    pulselength : float scalar, optional 
        Length of the longest microwave pulse in the sequence in microseconds. Used to determine the uncertainty in the 
        boundaries of the pathway refocusing times.

    Returns
    -------
    experiment : ``ExperimentInfo`` object
        Experiment object. Can be passed to ``dipolarmodel`` to introduce better
        constraints into the model.

    """

    # Theoretical refocusing times
    def reftimes(tau):
        tref = [ tau, 0]
        # Sort according to pathways order
        if pathways is not None:
            tref = [tref[pathway-1] for pathway in pathways]
        return tref
    # Pulse delays 
    delays = [tau]
    # Theoretical dipolar harmonics
    harmonics = [ 1, 1]
    # Pathway labels
    pathwaylabels = np.arange(1,len(harmonics)+1)

    # Sort according to pathways order
    if pathways is not None:
        _checkpathways(pathways,Nmax=len(harmonics))
        harmonics = [harmonics[pathway-1] for pathway in pathways] 
        pathwaylabels = [pathwaylabels[pathway-1] for pathway in pathways]

    return ExperimentInfo('3-pulse DEER', reftimes, harmonics, pulselength, pathwaylabels,delays)
#===============================================================================

#===============================================================================
def ex_4pdeer(tau1, tau2, pathways=[1,2,3,4], pulselength=0.016):
    r"""
    Generate a 4-pulse DEER dipolar experiment model. 
    
    The figure below shows the dipolar pathways in 4-pulse DEER. The observer (black) and pump (grey) pulses and their interpulse delays are shown on the top.
    The middle table summarizes all detectable modulated dipolar pathways `p` along their dipolar phase accumulation factors `\mathbf{s}_p`,
    harmonics `\delta_p` and refocusing times `t_{\mathrm{ref},p}`. The most commonly encountered pathways are highlighted in color. 
    The bottom panel shows a decomposition of the dipolar signal into the individual intramolecular contributions (shown as colored lines).
    
    The dipolar time axis is defined such that `t=0` right after the second observer pulse.

    .. figure:: ../images/ex_4pdeer_pathways.png
        :width: 350px

        Source: L. Fábregas-Ibáñez, Advanced data analysis and modeling in dipolar EPR spectroscopy, Doctoral dissertation, 2022


    Parameters 
    ----------

    tau1 : float scalar
        1st static interpulse delay `\tau_1`. 
   
    tau2 : float scalar
        2nd static interpulse delay `\tau_2`. 

    pathways :  array_like, optional 
        Pathways to include in the model. The pathways are specified as a list of pathways labels `p`. 
        By default, pathways 1-4 are included as shown in the table above.

    pulselength : float scalar, optional 
        Length of the longest microwave pulse in the sequence in microseconds. Used to determine the uncertainty in the 
        boundaries of the pathway refocusing times.

    Returns
    -------
    experiment : ``ExperimentInfo`` object
        Experiment object. Can be passed to ``dipolarmodel`` to introduce better constraints into the model.

    """
    # Theoretical refocusing times
    def reftimes(tau1,tau2):
        tref = [tau1, 
            tau1+tau2,
                    0,
                 tau2]
        # Sort according to pathways order
        if pathways is not None:
            tref = [tref[pathway-1] for pathway in pathways]
        return tref
    # Pulse delays 
    delays = [tau1,tau2]
    # Theoretical dipolar harmonics
    harmonics = [1, 1, 1, 1]
    # Pathway labels
    pathwaylabels = np.arange(1,len(harmonics)+1)

    # Sort according to pathways order
    if pathways is not None:
        _checkpathways(pathways,Nmax=len(harmonics))
        harmonics = [harmonics[pathway-1] for pathway in pathways] 
        pathwaylabels = [pathwaylabels[pathway-1] for pathway in pathways]

    return ExperimentInfo('4-pulse DEER', reftimes, harmonics, pulselength, pathwaylabels, delays)
#===============================================================================

#===============================================================================
def ex_rev5pdeer(tau1, tau2, tau3, pathways=[1,2,3,4,5], pulselength=0.016):
    r"""
    Generate a reverse 5-pulse DEER dipolar experiment model. 
        
    The figure below shows the dipolar pathways in reverse 5-pulse DEER. The observer (black) and pump (grey) pulses and their interpulse delays are shown on the top.
    The middle table summarizes all detectable modulated dipolar pathways `p` along their dipolar phase accumulation factors `\mathbf{s}_p`,
    harmonics `\delta_p` and refocusing times `t_{\mathrm{ref},p}`. The most commonly encountered pathways are highlighted in color. 
    The bottom panel shows a decomposition of the dipolar signal into the individual intramolecular contributions (shown as colored lines).

    .. figure:: ../images/ex_rev5pdeer_pathways.png
        :width: 350px

        Source: L. Fábregas-Ibáñez, Advanced data analysis and modeling in dipolar EPR spectroscopy, Doctoral dissertation, 2022


    Parameters 
    ----------

    tau1 : float scalar
        1st static interpulse delay `\tau_1`. 
   
    tau2 : float scalar
        2nd static interpulse delay `\tau_2`. 

    tau3 : float scalar
        3rd static interpulse delay `\tau_3`. 

    pathways :  array_like, optional 
        Pathways to include in the model. The pathways are specified as a list of pathways labels `p`. 
        By default, pathways 1-5 are included as shown in the table above.

    pulselength : float scalar, optional 
        Length of the longest microwave pulse in the sequence in microseconds. Used to determine the uncertainty in the 
        boundaries of the pathway refocusing times.

    Returns
    -------
    experiment : ``ExperimentInfo`` object
        Experiment object. Can be passed to ``dipolarmodel`` to introduce better
        constraints into the model.

    """

    # Theoretical refocusing times
    def reftimes(tau1,tau2,tau3):
        tref = [  tau3,
             tau2-tau3,
             tau1+tau3,
        tau1+tau2-tau3,
                  tau2,
             tau1+tau2,
                     0,
                  tau1]
        # Sort according to pathways order
        if pathways is not None:
            tref = [tref[pathway-1] for pathway in pathways]
        return tref
    # Pulse delays 
    delays = [tau1,tau2,tau3]
    # Theoretical dipolar harmonics
    harmonics = [ 1, 1, 1, 1, 1, 1, 1, 1]
    # Pathway labels
    pathwaylabels = np.arange(1,len(harmonics)+1)

    # Sort according to pathways order
    if pathways is not None:
        _checkpathways(pathways,Nmax=len(harmonics))
        harmonics = [harmonics[pathway-1] for pathway in pathways]
        pathwaylabels = [pathwaylabels[pathway-1] for pathway in pathways]

    return ExperimentInfo('Reverse 5-pulse DEER', reftimes, harmonics, pulselength, pathwaylabels, delays)
#===============================================================================


#===============================================================================
def ex_fwd5pdeer(tau1, tau2, tau3, pathways=[1,2,3,4,5], pulselength=0.016):
    r"""
    Generate a forward 5-pulse DEER dipolar experiment model. 
        
    The figure below shows the dipolar pathways in forward 5-pulse DEER. The observer (black) and pump (grey) pulses and their interpulse delays are shown on the top.
    The middle table summarizes all detectable modulated dipolar pathways `p` along their dipolar phase accumulation factors `\mathbf{s}_p`,
    harmonics `\delta_p` and refocusing times `t_{\mathrm{ref},p}`. The most commonly encountered pathways are highlighted in color. 
    The bottom panel shows a decomposition of the dipolar signal into the individual intramolecular contributions (shown as colored lines).

    .. figure:: ../images/ex_fwd5pdeer_pathways.png
        :width: 350px

        Source: L. Fábregas-Ibáñez, Advanced data analysis and modeling in dipolar EPR spectroscopy, Doctoral dissertation, 2022


    Parameters 
    ----------

    tau1 : float scalar
        1st static interpulse delay `\tau_1`. 
   
    tau2 : float scalar
        2nd static interpulse delay `\tau_2`. 

    tau3 : float scalar
        3rd static interpulse delay `\tau_3`. 

    pathways :  array_like, optional 
        Pathways to include in the model. The pathways are specified as a list of pathways labels `p`. 
        By default, pathways 1-5 are included as shown in the table above.

    pulselength : float scalar, optional 
        Length of the longest microwave pulse in the sequence in microseconds. Used to determine the uncertainty in the 
        boundaries of the pathway refocusing times.

    Returns
    -------
    experiment : ``ExperimentInfo`` object
        Experiment object. Can be passed to ``dipolarmodel`` to introduce better constraints into the model.


    """
    # Theoretical refocusing times
    def reftimes(tau1,tau2,tau3):
        tref = [  tau3, 
             tau1-tau3,
             tau2+tau3,
        tau1+tau2-tau3, 
                  tau1,
             tau1+tau2,
                     0,
                  tau2]
        # Sort according to pathways order
        if pathways is not None:
            tref = [tref[pathway-1] for pathway in pathways]
        return tref
    # Pulse delays 
    delays = [tau1,tau2,tau3]
    # Theoretical dipolar harmonics
    harmonics = [1, 1, 1, 1, 1, 1, 1, 1]
    # Pathway labels
    pathwaylabels = np.arange(1,len(harmonics)+1)

    # Sort according to pathways order
    if pathways is not None:
        _checkpathways(pathways,Nmax=len(harmonics))
        harmonics = [harmonics[pathway-1] for pathway in pathways] 
        pathwaylabels = [pathwaylabels[pathway-1] for pathway in pathways]

    return ExperimentInfo('Forward 5-pulse DEER', reftimes, harmonics, pulselength, pathwaylabels, delays)
#===============================================================================

#===============================================================================
def ex_sifter(tau1, tau2, pathways=[1,2,3], pulselength=0.016):
    r"""
    Generate a 4-pulse SIFTER dipolar experiment model. 
        
    The figure below shows the dipolar pathways in 4-pulse SIFTER. The pulses and the interpulse delays are shown on the top.
    The middle table summarizes all detectable modulated dipolar pathways `p` along their dipolar phase accumulation factors `\mathbf{s}_p`,
    harmonics `\delta_p` and refocusing times `t_{\mathrm{ref},p}`. The most commonly encountered pathways are highlighted in color. 
    The bottom panel shows a decomposition of the dipolar signal into the individual intramolecular contributions (shown as colored lines).

    .. figure:: ../images/ex_sifter_pathways.png
        :width: 350px

        Source: L. Fábregas-Ibáñez, Advanced data analysis and modeling in dipolar EPR spectroscopy, Doctoral dissertation, 2022


    Parameters 
    ----------

    tau1 : float scalar
        1st static interpulse delay `\tau_1`. 

    tau2 : float scalar
        2nd static interpulse delay `\tau_2`. 

    pathways :  array_like, optional 
        Pathways to include in the model. The pathways are specified as a list of pathways labels `p`. 
        By default, all 3 pathways are included as shown in the table above.

    pulselength : float scalar, optional 
        Length of the longest microwave pulse in the sequence in microseconds. Used to determine the uncertainty in the 
        boundaries of the pathway refocusing times.

    Returns
    -------
    experiment : ``ExperimentInfo`` object
        Experiment object. Can be passed to ``dipolarmodel`` to introduce better constraints into the model.

    """
    # Theoretical refocusing times
    def reftimes(tau1,tau2):
        tref = [ tau2-tau1, 
                    2*tau2,
                   -2*tau1]
        # Sort according to pathways order
        if pathways is not None:
            tref = [tref[pathway-1] for pathway in pathways]
        return tref
    # Pulse delays 
    delays = [tau1,tau2]
    # Theoretical dipolar harmonics
    harmonics = [ 1, 1/2, 1/2]
    # Pathway labels
    pathwaylabels = np.arange(1,len(harmonics)+1)

    # Sort according to pathways order
    if pathways is not None:
        _checkpathways(pathways,Nmax=len(harmonics))
        harmonics = [harmonics[pathway-1] for pathway in pathways]
        pathwaylabels = [pathwaylabels[pathway-1] for pathway in pathways]

    return ExperimentInfo('4-pulse SIFTER', reftimes, harmonics, pulselength, pathwaylabels, delays)
#===============================================================================


#===============================================================================
def ex_ridme(tau1, tau2, pathways=[1,2,3,4], pulselength=0.016):
    r"""
    Generate a 5-pulse RIDME dipolar experiment model. 
        
    The figure below shows the dipolar pathways in 5-pulse RIDME. The pulses and their interpulse delays are shown on the top.
    The middle table summarizes all detectable modulated dipolar pathways `p` along their dipolar phase accumulation factors `\mathbf{s}_p`,
    harmonics `\delta_p` and refocusing times `t_{\mathrm{ref},p}`. The most commonly encountered pathways are highlighted in color. 
    The bottom panel shows a decomposition of the dipolar signal into the individual intramolecular contributions (shown as colored lines).
    
    The dipolar time axis is defined such that `t=0` right after the second pulse.

    (Note that the model does not account for any relaxation-induced effects)

    .. figure:: ../images/ex_ridme_pathways.png
        :width: 350px

    Parameters 
    ----------

    tau1 : float scalar
        1st static interpulse delay `\tau_1`.

    tau2 : float scalar
        2nd static interpulse delay `\tau_2`.

    pathways :  array_like, optional 
        Pathways to include in the model. The pathways are specified as a list of pathways labels `p`. 
        By default, pathways 1-4 are included as shown in the table above.

    pulselength : float scalar, optional
        Length of the longest microwave pulse in the sequence, in microseconds. Used to determine
        the uncertainty in the boundaries of the pathway refocusing times.

    Returns
    -------
    experiment : ``ExperimentInfo`` object
        Experiment object. Can be passed to ``dipolarmodel`` to introduce better
        constraints into the model.

    """
    # Theoretical refocusing times
    def reftimes(tau1, tau2):
        tref = [tau1,
                tau1+tau2,
                0,
                tau2]
        # Sort according to pathways order
        if pathways is not None:
            tref = [tref[pathway-1] for pathway in pathways]
        return tref
    # Pulse delays 
    delays = [tau1, tau2]
    # Theoretical dipolar harmonics
    harmonics = [1, 1, 1, 1]
    # Pathway labels
    pathwaylabels = np.arange(1, len(harmonics)+1)

    # Sort according to pathways order
    if pathways is not None:
        _checkpathways(pathways, Nmax=len(harmonics))
        harmonics = [harmonics[pathway-1] for pathway in pathways] 
        pathwaylabels = [pathwaylabels[pathway-1] for pathway in pathways]

    return ExperimentInfo('5-pulse RIDME', reftimes, harmonics, pulselength, pathwaylabels, delays)
#===============================================================================


#===============================================================================
def ex_dqc(tau1, tau2, tau3, pathways=[1,2,3], pulselength=0.016):
    r"""
    Generate a 6-pulse DQC dipolar experiment model. 
        
    The figure below shows the dipolar pathways in 6-pulse DQC. The pulses and their interpulse delays are shown on the top.
    The middle table summarizes all detectable modulated dipolar pathways `p` along their dipolar phase accumulation factors `\mathbf{s}_p`,
    harmonics `\delta_p` and refocusing times `t_{\mathrm{ref},p}`. The most commonly encountered pathways are highlighted in color. 
    The bottom panel shows a decomposition of the dipolar signal into the individual intramolecular contributions (shown as colored lines).

    .. figure:: ../images/ex_dqc_pathways.png
        :width: 700px

        Source: L. Fábregas-Ibáñez, Advanced data analysis and modeling in dipolar EPR spectroscopy, Doctoral dissertation, 2022


    Parameters 
    ----------

    tau1 : float scalar
        1st static interpulse delay `\tau_1`. 

    tau2 : float scalar
        2nd static interpulse delay `\tau_2`. 

    tau3 : float scalar
        3rd static interpulse delay `\tau_3`. 

    pathways :  array_like, optional 
        Pathways to include in the model. The pathways are specified as a list of pathways labels `p`. 
        By default, pathways 1, 2, and 3 are included as shown in the table above.

    pulselength : float scalar, optional 
        Length of the longest microwave pulse in the sequence in microseconds. Used to determine the uncertainty in the 
        boundaries of the pathway refocusing times.

    Returns
    -------
    experiment : ``ExperimentInfo`` object
        Experiment object. Can be passed to ``dipolarmodel`` to introduce better constraints into the model.

    """
    # Theoretical refocusing times
    def reftimes(tau1,tau2,tau3):
        tref = [ tau2-tau1,
                    2*tau2, 
                   -2*tau1,
            tau2-tau1-tau3,
            -2*(tau1-tau3),
            -2*(tau1+tau3),
             2*(tau2+tau3),
             2*(tau1-tau3)]
        # Sort according to pathways order
        if pathways is not None:
            tref = [tref[pathway-1] for pathway in pathways]
        return tref
    # Pulse delays 
    delays = [tau1,tau2,tau3]
    # Theoretical dipolar harmonics
    harmonics = [ 1, 
                1/2,
                1/2,
                  1,
                1/2,
                1/2,
                1/2,
                1/2]
    # Pathway labels
    pathwaylabels = np.arange(1,len(harmonics)+1)

    # Sort according to pathways order
    if pathways is not None:
        _checkpathways(pathways,Nmax=len(harmonics))
        harmonics = [harmonics[pathway-1] for pathway in pathways]
        pathwaylabels = [pathwaylabels[pathway-1] for pathway in pathways]

    return ExperimentInfo('6-pulse DQC', reftimes, harmonics, pulselength, pathwaylabels, delays)
#===============================================================================