# model.py - DeerLab's modelling interface
# ---------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
from numpy.core.shape_base import atleast_1d, block
from scipy.sparse.construct import block_diag
from deerlab.solvers import snlls
from deerlab.classes import FitResult, UQResult
from deerlab.bootan import bootan
import inspect 
from copy import deepcopy
import difflib

#===================================================================================
class Parameter(): 
    r""" Represents a model parameter or a single parameter vector. 

    Attributes
    ----------
    name : string 
        Name/description of the parameter

    units : string 
        Physical units of the parameter

    par0 : float or array_like 
        Value at which to initialize the parameter at the start of a fit routine. 
        Must be specified in the model or latest upon fitting.

    lb : float or array_like 
        Lower bound of the parameter. If not specified it is assumed to unbounded.

    ub : float or array_like 
        Upper bound of the parameter. If not specified it is assumed to unbounded.

    linear : boolean 
        Describes whether the model behaves linearly with respect to the parameter.

    frozen : boolean 
        Describes whether the parameter will be frozen at a specific value during a fit.

    value : float
        Value at which the parameter will be frozen during a fit.


    Methods
    -------

    """

    #=======================================================================================
    #                                         Constructor
    #=======================================================================================

    #---------------------------------------------------------------------------------------
    def __init__(self, parent=None, idx=None, description=None, par0=None, frozen=False, lb=-np.inf, ub=np.inf,value=None, units=None, linear=False): 
        # Attributes
        self._parent = parent # Parent 
        self.idx = idx
        self.description = description # Description
        self.units = units    # Units
        self.par0 = par0      # Start values
        self.lb = lb          # Lower bounds
        self.ub = ub          # Upper bounds
        self.value = value
        self.frozen = frozen   # Frozen
        self.linear = linear  # Linearity
    #---------------------------------------------------------------------------------------

    #=======================================================================================
    #                                         Methods
    #=======================================================================================

    #---------------------------------------------------------------------------------------
    def set(self,**attributes):
        """
        Set multiple attributes for a parameter. 

        Parameters
        ----------
        attributes : keyword/values pairs
            Pairs of keywords defining the parameter attribute and the value to be assignes.
            For example ``parameter.set(lb=0,ub=1)`` to set both the ``lb`` and ``ub`` attributes.
        """
        for key in attributes:
            if not hasattr(self,key):
                raise AttributeError(f"'{key}' is not a valid parameter attribute.")
            setattr(self,key,attributes[key])
        return 
    #---------------------------------------------------------------------------------------
    
    #---------------------------------------------------------------------------------------
    def freeze(self,value):
        """
        Freeze a parameter during a fit/optimization.

        Parameters
        ----------
        value : float or array_like
            Value at which to freeze the parameter during optimization.
        """
        N = len(np.atleast_1d(self.frozen))
        if N>1:
            self.frozen = np.full(N,True)
            self.value = np.full(N,value)        
        else:
            self.frozen = True
            self.value = value        
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def unfreeze(self):
        """
        Release a frozen parameter during a fit/optimization.
        """
        N = len(np.atleast_1d(self.frozen))
        if N>1:
            self.frozen = np.full(N,False)
            self.value = np.full(N,None)        
        else:
            self.frozen = False
            self.value = None
    #---------------------------------------------------------------------------------------
#===================================================================================



class Model():
#===================================================================================


    #=======================================================================================
    #                                         Constructor
    #=======================================================================================

    #---------------------------------------------------------------------------------------
    def __init__(self,Amodel,constants=None,signature=None): 
        """
        Model object constructor
        """
        if not callable(Amodel):
            Amatrix = Amodel.copy()
            Amodel = lambda *_: Amatrix 
        self.nonlinmodel = Amodel
        self.description = None
        self._constantsInfo = []
        self.parents = None
        if signature is None:
            # Get list of parameter names from the function signature
            signature = inspect.getfullargspec(Amodel).args

        parameters = signature.copy()
        # Check if one of the arguments is an axis argument     
        if constants is not None:
            if not isinstance(constants,list):
                constants = [constants]
            for argname in constants:
                for n,par in enumerate(parameters): 
                    if par==argname: 
                        self._constantsInfo.append({"argkey":argname,'argidx':n})
            for argname in constants:
                for n,par in enumerate(parameters): 
                    if par==argname: 
                        parameters.remove(argname)
        Nconstants = len(self._constantsInfo)


        # Use a wrapper function to facilitate internal arguments manipulation        
        #-----------------------------------
        def model_with_constants(*inputargs):
            constants = inputargs[:Nconstants]
            θ = inputargs[Nconstants:]
            args = list(θ)
            if self._constantsInfo is not None:
                for info,constant in zip(self._constantsInfo,constants):
                    args.insert(info['argidx'],constant)
            return Amodel(*args)
        #-----------------------------------    
        self.nonlinmodel = model_with_constants

        # Update the number of parameters in the model
        self.Nparam = len(parameters)
        self.Nnonlin = len(parameters)
        self.Nlin = 0

        for n,param in enumerate(parameters): 
            newparam = Parameter(parent=self, idx=n)
            setattr(self,param,newparam)
        self.signature = signature
    #---------------------------------------------------------------------------------------

    # Gets called when an attribute is accessed
    #--------------------------------------------------------------------------------
    def __getattribute__(self, attr):
        try:
            return super(Model, self).__getattribute__(attr)
        except AttributeError:
            errstr = f'The model has no attribute {attr}.'
            attributes = [key for key in self.__dict__]
            proposal = difflib.get_close_matches(attr, attributes)
            if len(proposal)>0:
                errstr += f' \n Did you mean: {proposal} ?'
            raise AttributeError(errstr)
    #--------------------------------------------------------------------------------



    #=======================================================================================
    #                                    Private methods
    #=======================================================================================

    #---------------------------------------------------------------------------------------
    def _getNparents(self):
        if self.parents is None: 
            return 1
        else: 
            return len(self.parents)
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def _core_model(self,Amodel,θnonlin,θlin): 
        """ 
        Calculates the core model response ``y`` based on the mathematical expression ``y = A(θnonlin)@θlin``.
        """
        # Calculate the design matrix
        if len(θnonlin)>0:
            A = Amodel(*θnonlin)
        else:
            A = Amodel()

        # Ensure the proper matrix properties
        A = np.atleast_2d(A)
        θlin = np.atleast_1d(np.squeeze(θlin))

        # If there are no linear parameters defined
        if len(θlin)==0: 
            θlin = np.array([1])

        if A.shape[1]!=len(θlin): 
            A = A.T

        # Full model calculation 
        y = A@θlin
        
        return y
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def _parameter_list(self, order='alphabetical'):
        "Get the list of parameters defined in the model sorted alphabetically or by vector definition"
        if order=='alphabetical':
            keylist = [param for param in dir(self) if isinstance(getattr(self,param),Parameter)]
        elif order=='vector':
            keylist = [param for param in dir(self) if isinstance(getattr(self,param),Parameter)]
            # If there are any parameters in vector form...
            for n,key in enumerate(keylist): 
                if isinstance(getattr(self,key).idx,np.ndarray):
                    # ...insert the key string for the same number of linear parameters in that vector 
                    keylist = np.insert(keylist,n*np.ones(len(getattr(self,key).idx)-1,dtype=int),key)  
            keylist = self._vecsort(keylist)
        # Remove any duplicates
        keylist = list(dict.fromkeys(keylist))
        return keylist
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def _vecsort(self,list):
        "Sort vectorized parameters attributes from alphabetical ordering to vector indexing"
        list = np.squeeze(np.atleast_1d(list))
        indices = np.concatenate([np.atleast_1d(getattr(self,param).idx) for param in dir(self) if isinstance(getattr(self,param),Parameter)])
        orderedlist = np.atleast_1d(list.copy())
        orderedlist[indices] = list

        return orderedlist
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def _split_linear(self,variable):
        "Split a vector in non-linear and linear parameter subset vectors"
        variable = np.atleast_1d(variable)
        linear = np.concatenate([np.atleast_1d(getattr(self,param).linear) for param in dir(self) if isinstance(getattr(self,param),Parameter)])
        linear = self._vecsort(linear)
        variable_nonlin = variable[~linear]
        variable_lin = variable[linear]
        return variable_lin, variable_nonlin
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def _getvector(self,attribute):
        "Get the list of parameters attributes defined in the model sorted alphabetically"
        return np.concatenate([np.atleast_1d(getattr(getattr(self,param),attribute)) for param in dir(self) if isinstance(getattr(self,param),Parameter)])
    #---------------------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    def _getparamuq(self,uq_full,paramidx):
        "Get the uncertainty quantification of an individual parameter"
        subset_model = lambda x: x[paramidx]
        param_lb =  self._vecsort(self._getvector('lb'))[paramidx]
        param_ub =  self._vecsort(self._getvector('ub'))[paramidx]
        frozen = self._vecsort(self._getvector('frozen'))[paramidx]
        if np.all(frozen): 
            param_uq = UQResult('void')
        else:
            param_uq = uq_full.propagate(subset_model,lbm=param_lb, ubm=param_ub)

        return param_uq
    #-----------------------------------------------------------------------------

    #=======================================================================================
    #                                         Methods
    #=======================================================================================

    #---------------------------------------------------------------------------------------
    def addlinear(self, key, vec=1, lb=-np.inf, ub=np.inf, par0=None, name=None, units=None, description=None):
        """
        Add a new linear :ref:`Parameter` object. 

        Parameters
        ----------
        key : string
            identifier of the parameter.

        vec : int scalar, optional
            Number of elements in the parameter. If ``vec>1`` then the parameter will represent a 
            vector of linear parameters of length ``vec``. 

        lb : float or array_like, optional
            Lower bound of the parameter. For vectorized parameters, must be a vector with ``vec`` elements. 

        ub : float or array_like, optional
            Lower bound of the parameter. For vectorized parameters, must be a vector with ``vec`` elements. 

        name : string, optional 
            Name/descriptrion of the parameter. 

        units : string, optional
            Physical units of the parameter.
        """
        if vec>1: 
            idx = np.arange(self.Nparam,self.Nparam+vec) 
            self.Nparam += vec        
            self.Nlin += vec
            newparam = Parameter(linear=np.full(vec,True), parent=self, idx=idx, par0=np.full(vec,par0), lb=np.full(vec,lb), ub=np.full(vec,ub), value=np.full(vec,None),frozen=np.full(vec,False), units=units, description=description)
        else:
            idx = self.Nparam
            self.Nparam += 1
            self.Nlin += 1
            newparam = Parameter(linear=True, parent=self, idx=idx, par0=par0, lb=lb, ub=ub, units=units, description=description)
        setattr(self,key,newparam)
        self.signature.append(key)
    #---------------------------------------------------------------------------------------

    def __call__(self,*args,**kargs):
    #---------------------------------------------------------------------------------------
        if kargs and args: 
            raise SyntaxError('The model must be called either with positional or keyword arguments. Not both.')
        
        Nrequired = len(self._parameter_list())
        Nrequired += len(self._constantsInfo)

        if len(args)!=Nrequired and len(kargs)!=Nrequired:
            raise SyntaxError(f'The model requires {Nrequired} arguments, but {len(args)+len(kargs)} have been specified.')

        if args:         
            constants= [np.atleast_1d(args[info['argidx']]) for info in self._constantsInfo]
            args_list = [np.atleast_1d(arg) for arg in args]
            args_list = [np.atleast_1d(arg) for arg in args_list  if not np.any([np.all(arg==constant) for constant in constants]) ]
            # Values are already ordered
            θ = np.concatenate(args_list) 
        elif kargs:
            constants = [np.atleast_1d(kargs[info['argkey']]) for info in self._constantsInfo]
            args_list = [np.atleast_1d(kargs[param]) for param in self._parameter_list(order='vector')]
            args_list = [ arg for arg in args_list  if not np.any([np.all(arg==constant) for constant in constants]) ]
            # Values must be ordered
            θ = np.concatenate(args_list)

        if len(θ)!=self.Nparam:
            raise SyntaxError(f'The model requires {self.Nparam} parameters, but {len(θ)} were specified.')   
            

        # Determine which parameters are linear and which nonlinear
        θlin, θnonlin = self._split_linear(θ)

        # Evaluate the core model
        y = self._core_model(lambda *θ: self.nonlinmodel(*constants,*θ),θnonlin,θlin)

        # Evaluate whether the response has 
        if hasattr(self,'_posteval_fcn'):
            y = self._posteval_fcn(y,*constants,*θ)
        return y
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def getmetadata(self):
        """
        Utility function to quickly request all metadata attributes of the model in vector form. 
        All elements are sorted according to the model function signature.
        """
        return {
            'names': self._parameter_list(order='vector'),
            'ub' : self._vecsort(self._getvector('ub')),
            'lb' : self._vecsort(self._getvector('lb')),
            'par0' : self._vecsort(self._getvector('par0')),
            'frozen' : self._vecsort(self._getvector('frozen')),
            'linear' : self._vecsort(self._getvector('linear')),
            'values' : self._vecsort(self._getvector('value')),
            'units' : self._vecsort(self._getvector('units')),
            }
    #---------------------------------------------------------------------------------------
#===================================================================================




#==============================================================================================
def fit(model,y,*constants,par0=None,bootstrap=0,**kwargs):
    r"""
    Fit the model to the data ``y`` via one of the three following approaches: 
    
    - Non-linear least-squares 
    - Regularized linear-least-squares 
    - Separable non-linear least-squares 

    The most appropiate solver is chosen automatically based on the model structure. 

    Parameters
    ----------
    y : array_like 
        Data to be fitted. 
    par0 : array_like, optional 
        Value at which to initialize the parameter at the start of a fit routine. 
        Must be specified if not defined in the model. Otherwise, it overrides the definition in the model. 

    Any additional solver-specific keyword arguments can be specified. See :ref:`nlls`, :ref:`rlls` and :ref:`snlls` for a full reference. 

    Returns
    -------
    :ref:`FitResult` with the following fields defined:
    nonlin : ndarray
        Fitted non-linear parameters.
    lin : ndarray
        Fitted linear parameters.
    model : ndarray
        Fitted model.
    nonlinUncert : :ref:`UQResult`
        Uncertainty quantification of the non-linear parameter set.
    linUncert : :ref:`UQResult`
        Uncertainty quantification of the linear parameter set.
    modelUncert : :ref:`UQResult`
        Uncertainty quantification of the fitted model.
    regparam : scalar
        Regularization parameter value used for the regularization of the linear parameters.
    plot : callable
        Function to display the results. It will display the fitted data.
        The function returns the figure object (``matplotlib.figure.Figure``)
        object as output, which can be modified. Using ``fig = plot(show=False)`` 
        will not render the figure unless ``display(fig)`` is called. 
    stats : dict
        Goodness of fit statistical estimators

        * ``stats['chi2red']`` - Reduced \chi^2 test
        * ``stats['r2']`` - R^2 test
        * ``stats['rmsd']`` - Root-mean squared deviation (RMSD)
        * ``stats['aic']`` - Akaike information criterion
        * ``stats['aicc']`` - Corrected Akaike information criterion
        * ``stats['bic']`` - Bayesian information criterion
    cost : float
        Value of the cost function at the solution.
    """

    if not isinstance(model,Model):
        raise TypeError('The input model must be a valid deerlab.Model object.')

    if len(constants)>0:
        constants = np.atleast_1d(constants)
        
    # Get boundaries and conditions for the linear and nonlinear parameters
    ubl,ub = model._split_linear(model._vecsort(model._getvector('ub')))
    lbl,lb = model._split_linear(model._vecsort(model._getvector('lb')))
    frozenl,frozen = model._split_linear(model._vecsort(model._getvector('frozen')))
    valuesl,values = model._split_linear(model._vecsort(model._getvector('value')))

    # Check the initial conditions and whether they are defined
    if par0 is None:
        _,par0 = model._split_linear(model._vecsort(model._getvector('par0')))
    if np.any(par0==None):
        raise RuntimeError(f"It appears some start values (par0) have not been specified. Either specify them in the model definition or using the keyword.")

    linfrozen = np.full(model.Nlin,None)
    linfrozen[frozenl] = valuesl[frozenl]
    nonlinfrozen = np.full(model.Nnonlin,None)
    nonlinfrozen[frozen] = values[frozen]

    if len(linfrozen)==0: linfrozen = [1]

    if model._getNparents()==1:
        y = [y]            
    nprev = 0
    ysubsets = []
    for yset in y:
        ysubsets.append(np.arange(nprev,nprev+len(yset)))
        nprev = nprev+len(yset)
    y = np.concatenate(y)


    if model.Nlin==0 and model.Nnonlin==0:
        raise AssertionError(f'The model has no parameters to fit.')    

    # Prepare the separable non-linear least-squares solver
    Amodel_fcn = lambda param: model.nonlinmodel(*constants,*param)
    fitfcn = lambda y: snlls(y,Amodel_fcn,par0,lb,ub,lbl,ubl,subsets=ysubsets,lin_frozen=linfrozen,nonlin_frozen=nonlinfrozen,**kwargs)        

    # Run the fitting algorithm 
    fitresults = fitfcn(y)

    # If requested, perform a bootstrap analysis
    if bootstrap>0: 
        def bootstrap_fcn(ysim): 
            fit = fitfcn(ysim)
            return fit.param,fit.model
        # Bootstrapped uncertainty quantification
        param_uq = bootan(bootstrap_fcn,y,fitresults.model,samples=bootstrap)
        # Substitute the fitted values by the bootsrapped median estimate
        fitresults.param = param_uq[0].median
        fitresults.paramUncert = param_uq[0]
        fitresults.model = param_uq[1].median
        fitresults.modelUncert = param_uq[1]

    # Get some basic information on the parameter vector
    keys = model._parameter_list(order='vector')
    param_idx = [np.atleast_1d(getattr(model,param).idx) for param in model._parameter_list('vector')]
    # Dictionary of parameter names and fitted values
    FitResult_param = {key : fitvalue for key,fitvalue in zip(keys,[fitresults.param[idx] for idx in param_idx])}
    # Dictionary of parameter names and fit uncertainties
    FitResult_paramuq = {f'{key}Uncert': model._getparamuq(fitresults.paramUncert,idx) for key,idx in zip(keys,param_idx)}
    # Dictionary of other fit quantities of interest
    FitResult_dict = {key: getattr(fitresults,key) for key in ['model','modelUncert','cost','plot','residuals']}
    # Generate FitResult object from all the dictionaries
    fit = FitResult({**FitResult_param,**FitResult_paramuq, **FitResult_dict }) 

    return fit
#==============================================================================================

def _importparameter(parameter):
    return {
        'lb' : parameter.lb,
        'ub' : parameter.ub,
        'par0' : parameter.par0,
        'description' : parameter.description,
        'units' : parameter.units,
        'frozen' : parameter.frozen,
        'value' : parameter.value,
    }

def _aresame(obj1,obj2):
    a = obj1.__dict__
    a = {key:val for key, val in a.items() if key not in ['_parent','idx']}
    b = obj2.__dict__
    b = {key:val for key, val in b.items() if key not in ['_parent','idx']}
    return a==b
# ==============================================================================
def link(model,**links):

    def _linkparameter(model,parameters,newname):
    # ---------------------------------------------------------------------  
        # Get a list of parameter names in the model
        model_parameters = model._parameter_list(order='vector')

        link_parameters,link_indices = [],[]
        for param in parameters:
            if not isinstance(param,Parameter) and isinstance(param,str):
                param = getattr(model,param)
            # Make list of parameter objects
            link_parameters.append(param)
            # Make list of parameter indices
            link_indices.append(param.idx)

        # Get the names of the parameters to be linked
        link_names = []
        for param in parameters:
            if isinstance(param,Parameter):
                for mparam in model_parameters:
                    if _aresame(param,getattr(model,mparam)):
                        link_names.append(mparam)
            else: 
                link_names.append(param)
        # Remove the first from the list as it will be kept
        link_names.pop(0)

        # Get the vector index of the parameter to be linked to
        nlink = link_indices[0]

        nnew = 0
        mapping = np.zeros(model.Nnonlin,dtype=int)
        mapping_linear = np.zeros(model.Nlin,dtype=int)
        Nnl = model.Nnonlin
        # Loop over all parameters in the model
        for n,param in enumerate(model_parameters):
            # If the parameter is to be linked...
            if param in link_names:
                # Update the number of parameters in the model
                model.Nparam -= 1
                if getattr(model,param).linear:
                    model.Nlin -= 1
                else: 
                    model.Nnonlin -= 1
                # Delete the linked parameter from the model
                if not getattr(model,param).linear:
                    # Update the parameter vector map
                    mapping[n] = nlink
                else: 
                    mapping_linear[n-Nnl] = nlink - Nnl                    
                delattr(model,param)

            # Otherwise...
            else:
                # Update the index of the parameter in the new vector 
                getattr(model,param).idx = nnew
                if not getattr(model,param).linear:
                    mapping[n] = nnew
                else: 
                    mapping_linear[n-Nnl] = nnew - Nnl  
                nnew += 1

        # Delete the old copy with the old name
        paramobj = getattr(model,model_parameters[nlink])
        delattr(model,model_parameters[nlink])
        # Create a copy of the linked parameter with the new name
        setattr(model,newname,paramobj)

        # Wrap the evaluation function         
        nonlinfcn = model.nonlinmodel

        linear_reduce_idx = [np.where(mapping_linear==n)[0].tolist() for n in np.unique(mapping_linear) ]
        Nconstants = len(model._constantsInfo)
        # ---------------------------------------------------------------------
        def linked_model_with_constants(*inputargs):
            # Redistribute the input parameter vector according to the mapping vector
            constants = inputargs[:Nconstants]
            θ = inputargs[Nconstants:]
            θ = np.atleast_1d(θ)[mapping]
            args = list(θ)
            if model._constantsInfo is not None:
                for info,constant in zip(model._constantsInfo,constants):
                    args.insert(info['argidx'],constant)                
            A = nonlinfcn(*args)
            Amapped = np.vstack([np.sum(np.atleast_2d(A[:,idx]),axis=1) for idx in linear_reduce_idx]).T
            return Amapped
        # ---------------------------------------------------------------------
        model.nonlinmodel = linked_model_with_constants

        # Return the updated model with the linked parameters
        return model
    # ---------------------------------------------------------------------

    if not isinstance(model,Model):
        raise TypeError('The first argument must be a Model object.')
    newmodel = deepcopy(model)
    for link_name in links: 
        newmodel = _linkparameter(newmodel,links[link_name],link_name)
    return newmodel
#==============================================================================================


#==============================================================================================
def combine(*inputmodels,**links): 

    # Initialize empty containers
    subsets_nonlin,arguments,arelinear = [],[],[]
    nprev = 0

    # Make deep-copies of the models to avoid modifying them
    models = [deepcopy(model) for model in inputmodels]

    # Loop over all models to be combined
    for n,model in enumerate(models): 
        
        # If one of the models has linear parameters, but not the others
        # add a dummy unity linear parameter 
        if model.Nlin==0:
            model.addlinear('scale',par0=1,lb=0)

        # Determine the subset of parameters for the current model
        subset = np.arange(nprev,nprev+model.Nnonlin,1)
        nprev += model.Nnonlin
        # From that subset, determine the non-linear subset
        subset_nonlin = subset[np.arange(model.Nnonlin)]
        subsets_nonlin.append(subset_nonlin)

        # Determine which parameters are linear
        arelinear = np.concatenate([arelinear,model._vecsort(model._getvector('linear'))])

        newarguments = model._parameter_list(order='vector') 
        # If there is more than one model, append a string to identify the origin
        if len(models)>1: 
            newarguments = [arg+f'_{n+1}' for arg in newarguments] 

        # Map the link parameter objects to the new model names for later
        oldargs = inputmodels[n]._parameter_list(order='vector')
        newarguments_ = newarguments.copy()
        for linkkey in links: 
            for i,param in enumerate(links[linkkey]):        
                for oldarg,newarg in zip(oldargs,newarguments_):
                    if getattr(inputmodels[n],oldarg)==param:
                        links[linkkey][i] = newarg 
                        newarguments_.remove(newarg)
                        oldargs.remove(oldarg)


        # Add the submodel arguments to the combined model signature
        arguments += newarguments


    # Preparation of the combined model signature    
    arelinear = np.asarray(arelinear).astype(bool)
    arguments = np.array(arguments)
    lin_params = arguments[arelinear]
    nonlin_params = arguments[~arelinear]

    # Account for the axis argument
    constants = []
    const_subsets = []
    for n,model in enumerate(models):
        subset = []
        for info in model._constantsInfo:
            constants.append(f"{info['argkey']}_{n+1}")
            subset.append(len(constants)-1)
        const_subsets.append(subset)
    signature = np.insert(nonlin_params,0,constants).tolist()
    Nconst = len(constants)
    Nsubsets = len(models)

    #---------------------------------------------------------------------
    def _combined_nonlinmodel(*inputargs):
        """Evaluates the nonlinear functions of the submodels and 
        concatenates them into a single design matrix"""
        constants = inputargs[:Nconst]
        param = inputargs[Nconst:]

        param = np.atleast_1d(param)
        constants = np.atleast_1d(constants)
        # Loop over the submodels in the model
        Amatrices = []
        for n,model in enumerate(models):
            # Evaluate the submodel
            Amatrix = np.atleast_2d(model.nonlinmodel(*constants[const_subsets[n],:],*param[subsets_nonlin[n]]))
            if np.shape(Amatrix)[1]!=model.Nlin: Amatrix = Amatrix.T
            Amatrices.append(Amatrix)
        
        model._subsets = 'adfsdfdafds'

        Anonlin_full = block_diag(Amatrices).toarray()

        if not any(arelinear):
            Anonlin_full = np.sum(Anonlin_full,1)

        return Anonlin_full
    #---------------------------------------------------------------------
    
    #---------------------------------------------------------------------
    def _split_output(y,*inputargs):
        constants = inputargs[:Nconst]
        param = inputargs[Nconst:]

        param = np.atleast_1d(param)
        constants = np.atleast_1d(constants)
        # Loop over the submodels in the model
        Amatrices = []
        for n,model in enumerate(models):
            # Evaluate the submodel
            Amatrix = np.atleast_2d(model.nonlinmodel(*constants[const_subsets[n],:],*param[subsets_nonlin[n]]))
            if np.shape(Amatrix)[1]!=model.Nlin: Amatrix = Amatrix.T
            Amatrices.append(Amatrix)         
        nprev = 0
        ysubsets = []
        for A in Amatrices:
            ysubsets.append(np.arange(nprev,nprev+A.shape[0]))
            nprev = nprev+A.shape[0]
        return [y[ysubsets[n]] for n in range(len(Amatrices))]
    #---------------------------------------------------------------------

    # Create the model object
    combinedModel = Model(_combined_nonlinmodel,constants=constants,signature=signature)

    # Add parent models 
    combinedModel.parents = models

    # Add post-evalution function for splitting of the call outputs
    setattr(combinedModel,'_posteval_fcn',_split_output) 

    # Add the linear parameters from the subset models    
    for lparam in lin_params:
        combinedModel.addlinear(lparam)

    parameters = np.concatenate([arguments,lin_params])
    # Import all parameter information from the subset models
    for name,param in zip(combinedModel._parameter_list(order='vector'),parameters):
        if '_' in name:
            param,n = name.split('_')
            n = int(n)-1
        else:
            param,n = name,0
        getattr(combinedModel,name).set(**_importparameter(getattr(models[n],param)))

    # Perform links of parameters if requested
    combinedModel = link(combinedModel,**links)

    # Return the new combined model object
    return combinedModel    
#==============================================================================================