# model.py - DeerLab's modelling interface
# ---------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2021: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
from scipy.sparse.construct import block_diag
from scipy.optimize import fminbound
from deerlab.solvers import snlls
from deerlab.correctphase import correctphase
from deerlab.classes import FitResult, UQResult
from deerlab.bootan import bootan
from deerlab.noiselevel import noiselevel
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
    def __init__(self, name=None, parent=None, idx=None, description=None, par0=None, frozen=False, lb=-np.inf, ub=np.inf,value=None, units=None, linear=False): 
        # Attributes
        self.name = name
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
    r"""Represents a model.

    Attributes
    ----------
    <parameter_name> : :ref:`Parameter`
        Model parameter object. 
    description : string 
        Description of the model.
    signature : string 
        Call signature (keyword arguments in order) of the model. 
    nonlinmodel : callable
        Function of the non-linear part of the model.
    Nnonlin : int scalar
        Number of non-linear parameters in the model.
    Nlin : int scalar
        Number of linear parameters in the model.
    Nparam : int scalar
        Number of parameters in the model.

    Methods
    -------

    """


    #=======================================================================================
    #                                         Constructor
    #=======================================================================================

    #---------------------------------------------------------------------------------------
    def __init__(self,Amodel,constants=None,signature=None): 
        """
        Construct a new model from a non-linear function. 

        Parameters
        ----------

        Amodel : callable 
            Function that takes a set of non-linear parameters and 
            returns either a the full model response or the design matrix 
            of the model response. A parameter will be added to the new model for each 
            input argument defined in the function signature.   
        
        constants : string or list thereof
            Names of the constant variables (non-parameter variables) taken
            by the ``Amodel`` function. These will not be added as parameters to the new model.

        signature : list of strings
            Signature of the ``Amodel`` function to manually specify the names of the input arguments.

        Returns
        -------

        model : `Model` object 
            Model object instance that takes the parameters defined for ``Amodel`` and returns the output of ``Amodel``.

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

        self._nonlinsignature = signature.copy()
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
            newparam = Parameter(parent=self, idx=n, name=param)
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
            n = 0
            for key in keylist: 
                if isinstance(getattr(self,key).idx,np.ndarray):
                    # ...insert the key string for the same number of linear parameters in that vector 
                    keylist = np.insert(keylist,n*np.ones(len(np.atleast_1d(getattr(self,key).idx))-1,dtype=int),key)  
                n += len(np.atleast_1d(getattr(self,key).idx))

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
    def _merge_linear(self,variable_nonlin,variable_lin):
        "Merge a vector's non-linear and linear parameter subset vectors"
        variable = np.zeros(len(variable_nonlin)+len(variable_lin))
        linear = np.concatenate([np.atleast_1d(getattr(self,param).linear) for param in dir(self) if isinstance(getattr(self,param),Parameter)])
        linear = self._vecsort(linear)
        variable[~linear] = variable_nonlin
        variable[linear] = variable_lin
        return variable
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
            if uq_full.type=='covariance':
                param_uq = uq_full.propagate(subset_model,lb=param_lb, ub=param_ub)
            elif uq_full.type=='bootstrap':
                param_uq = UQResult('bootstrap',data=uq_full.samples[:,paramidx],lb=param_lb, ub=param_ub)
            else:
                param_uq = UQResult('void')
        return param_uq
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    def _check_if_already_exists(self,key):
        if hasattr(self,key):
            raise KeyError(f'The model already has a "{key}" parameter.')
    #-----------------------------------------------------------------------------
    
    
    #=======================================================================================
    #                                         Methods
    #=======================================================================================

    #---------------------------------------------------------------------------------------
    def addnonlinear(self, key, lb=-np.inf, ub=np.inf, par0=None, name=None, units=None, description=None):
        """
        Add a new non-linear :ref:`Parameter` object. 

        Parameters
        ----------
        key : string
            identifier of the parameter.

        lb : float or array_like, optional
            Lower bound of the parameter. For vectorized parameters, must be a vector with ``vec`` elements. 

        ub : float or array_like, optional
            Lower bound of the parameter. For vectorized parameters, must be a vector with ``vec`` elements. 

        description : string, optional 
            Name/descriptrion of the parameter. 

        units : string, optional
            Physical units of the parameter.
        """
        self._check_if_already_exists(key)
        idx = self.Nparam
        self.Nparam += 1
        self.Nnonlin += 1
        newparam = Parameter(name=key, linear=False, parent=self, idx=idx, par0=par0, lb=lb, ub=ub, units=units, description=description)
        setattr(self,key,newparam)
        Nconstants = len(self._constantsInfo)
        Amodel = self.nonlinmodel
        topop = self.Nnonlin-1
        #------------------------------------------------
        def model_with_constants_and_added_nonlin(*inputargs):
            constants = inputargs[:Nconstants]
            θ = inputargs[Nconstants:]
            args = list(θ)
            args.pop(topop)
            if self._constantsInfo is not None:
                for info,constant in zip(self._constantsInfo,constants):
                    args.insert(info['argidx'],constant)
            return Amodel(*args)
        #------------------------------------------------
        self.nonlinmodel = model_with_constants_and_added_nonlin
        self.signature.append(key)
    #---------------------------------------------------------------------------------------


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

        description : string, optional 
            Name/descriptrion of the parameter. 

        units : string, optional
            Physical units of the parameter.
        """
        self._check_if_already_exists(key)
        if vec>1: 
            idx = np.arange(self.Nparam,self.Nparam+vec) 
            self.Nparam += vec        
            self.Nlin += vec
            newparam = Parameter(name=key, linear=np.full(vec,True), parent=self, idx=idx, par0=np.full(vec,par0), lb=np.full(vec,lb), ub=np.full(vec,ub), value=np.full(vec,None),frozen=np.full(vec,False), units=units, description=description)
        else:
            idx = self.Nparam
            self.Nparam += 1
            self.Nlin += 1
            newparam = Parameter(name=key, linear=True, parent=self, idx=idx, par0=par0, lb=lb, ub=ub, units=units, description=description)
        setattr(self,key,newparam)
        self.signature.append(key)
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def copy(self): 
        "Generate a deep copy of the model."
        return deepcopy(self)
    #---------------------------------------------------------------------------------------

    def __call__(self,*args,**kargs):
    #---------------------------------------------------------------------------------------
        """
        Evaluate the model for a given set of parameters. 

        Takes the constant variables and (non-linear and linear) parameter variables as positional
        or keyword arguments and evaluateds the model.
        """
        # Check that the correct number of arguments have been specified
        Nrequired = len(self._parameter_list())
        Nrequired += len(self._constantsInfo)
        if (len(args)+len(kargs))!=Nrequired:
            raise SyntaxError(f'The model requires {Nrequired} arguments, but {len(args)+len(kargs)} have been specified.')

        # Positional arguments       
        args_constants= [np.atleast_1d(args[info['argidx']]) for info in self._constantsInfo if info['argidx']<len(args)]
        args_list = [np.atleast_1d(arg) for idx,arg in enumerate(args) if idx not in [info['argidx'] for info in self._constantsInfo]]

        # Keywords arguments
        kargs_constants = [np.atleast_1d(kargs[info['argkey']]) for info in self._constantsInfo  if info['argkey'] in kargs] 
        kargs_list = [np.atleast_1d(kargs[param]) for param in self._parameter_list(order='vector')[len(args_list):]]
        
        constants = args_constants + kargs_constants
        param_list = args_list + kargs_list

        # Concatente all parameter into a single vector
        θ = np.concatenate(param_list)

        # Check that all parameters have been passed
        if len(θ)!=self.Nparam:
            raise SyntaxError(f'The model requires {self.Nparam} parameters, but {len(args_list)} were specified.')   

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

        
    #---------------------------------------------------------------------------------------
    def addregularization(self,functional='aic',description=None):
        """ 
        Add a new :ref:`Regularization` object to impose Tikhonov regularization upon 
        the linear parameters of the model.

        Parameters
        ----------
        key : string
            Identifier of the parameter.

        functional : string 
            Method for the automatic selection of the optimal regularization parameter/weight:

            * ``'lr'`` - L-curve minimum-radius method (LR)
            * ``'lc'`` - L-curve maximum-curvature method (LC)
            * ``'cv'`` - Cross validation (CV)
            * ``'gcv'`` - Generalized Cross Validation (GCV)
            * ``'rgcv'`` - Robust Generalized Cross Validation (rGCV)
            * ``'srgcv'`` - Strong Robust Generalized Cross Validation (srGCV)
            * ``'aic'`` - Akaike information criterion (AIC)
            * ``'bic'`` - Bayesian information criterion (BIC)
            * ``'aicc'`` - Corrected Akaike information criterion (AICC)
            * ``'rm'`` - Residual method (RM)
            * ``'ee'`` - Extrapolated Error (EE)
            * ``'ncp'`` - Normalized Cumulative Periodogram (NCP)
            * ``'gml'`` - Generalized Maximum Likelihood (GML)
            * ``'mcl'`` - Mallows' C_L (MCL)
        
            The default ``'aic'``.

        description : string, optional 
            Description of the penalty. 
        """
        if np.any([isinstance(getattr(self,attr),Regularization) for attr in dir(self)]):
            raise RuntimeError('The model already has a Regularization attribute.')
        self.regularization = Regularization(functional,description)
    #---------------------------------------------------------------------------------------
        

    #---------------------------------------------------------------------------------------
    def addpenalty(self,key,penaltyfcn,functional,**kargs):
        """
        Add a new :ref:`Penalty` object. 

        Parameters
        ----------
        key : string
            Identifier of the parameter.

        penaltyfcn : callable 
            Penalty functiona taking model parameters as inputs. Must return a vector which will be squared and appended 
            to the residual vector during the fitting of the model. A penalty weight will be added automatically to the 
            output vector.

        functional : string 
            Selection functional for optimization of the penalty weight. Must be a string from the following: 

            * ``'aic'`` - Akaike information criterion (AIC)
            * ``'aicc'`` - Corrected Akaike information criterion (AICc)
            * ``'bic'`` - Bayesian complexity criterion (BIC)
            * ``'icc'`` - Informational complexity criterion (ICC) 

        description : string, optional 
            Description of the penalty. 
        """        

        # Create penalty object
        penalty = Penalty(penaltyfcn,functional,**kargs)
        for arg in penalty.signature: 
            if arg not in self._parameter_list():
                raise KeyError(f'The penalty argument {arg} does not correspond to any valid model parameter.')

        # Get the parameter names of the model
        modelparam = self._parameter_list('vector')
        # Determine the indices of the subset of parameters the model depends on
        subsets = [getattr(self,modelparam[np.where(np.asarray(modelparam)==param)[0][0]]).idx for param in penalty.signature]
        # Adapt the signature of penaltyfcn for snlls
        penaltyfcn_ = penalty.penaltyfcn
        penalty.penaltyfcn = lambda pnonlin,plin,weight: penaltyfcn_(weight,*[np.concatenate([pnonlin,plin])[subset] for subset in subsets])

        # Store penalty object as model attribute
        setattr(self,key,penalty)
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def _parameter_table(self):
        string = inspect.cleandoc(f"""
    Model information 
    -----------------

    Model description: {self.description}
    Model call signature: {','.join(self.signature)}
    Constants: {[entry['argkey'] for entry in self._constantsInfo]}

    Parameter Table 
    ---------------

    ============ ========= ========== =========== ======== ========== ==========================
        Name       Lower     Upper      Type       Frozen   Units     Description  
    ============ ========= ========== =========== ======== ========== ==========================""")
        for n,paramname in enumerate(self._parameter_list(order='vector')): 
            string += f'\n   {paramname:7s}'
            string += f'     {np.atleast_1d(getattr(self,paramname).lb)[0]:5.3g}'
            string += f'     {np.atleast_1d(getattr(self,paramname).ub)[0]:5.3g}'
            string += f'      {"linear" if np.all(getattr(self,paramname).linear) else "nonlin"}'
            string += f'      {"Yes" if np.all(getattr(self,paramname).frozen) else "No":3s}'
            string += f'       {str(getattr(self,paramname).units):6s}'
            string += f'   {str(getattr(self,paramname).description):s}'
        string += f'\n============ ========= ========== =========== ======== ========== =========================='
        return string
    #---------------------------------------------------------------------------------------

    def __str__(self):
        return self._parameter_table()
    def __repr__(self):
        return self._parameter_table()        

#===================================================================================

#==============================================================================
class Regularization():
    def __init__(self,functional,description):
        self.weight = Parameter(parent=self, idx=0, name='weight',description='Regularization parameter/weight')
        self.weight.lb = np.finfo(float).eps
        self.weight.ub = 1/np.finfo(float).eps
        delattr(self.weight,'par0')
        delattr(self.weight,'linear')

        self.selection = functional 
        self.description = description
#==============================================================================


#==============================================================================
class Penalty():

    #--------------------------------------------------------------------------
    def __init__(self,penaltyfcn,selection,description=None,signature=None):

        #-------------------------------------------------------------------------------
        def selectionfunctional(fitfcn,y,sigma,log10weight):
            # Penalty weight: linear-scale -> log-scale
            weight = 10**log10weight

            # Run the fit
            fitresult = fitfcn(weight)

            if selection=='icc':
                # Get the fitted model
                yfit = fitresult.model

                # Get non-linear parameters covariance submatrix
                fitpars = fitresult.nonlin + 1e-16
                covmat = fitresult.nonlinUncert.covmat
                covmat = covmat/(fitpars[np.newaxis,:]*fitpars[:,np.newaxis])

                # Informational complexity criterion (ICC)
                if not np.all(covmat==0):
                    icc = np.sum((y - yfit)**2/sigma**2) + np.sum(np.log(np.diag(covmat))) + np.linalg.slogdet(covmat)[1]
                else:
                    icc = np.sum((y - yfit)**2/sigma**2)
                return icc 

            elif selection=='aic':
                aic = fitresult.stats['aic']
                return aic

            elif selection=='aicc':
                aicc = fitresult.stats['aicc']
                return aicc

            elif selection=='bic':
                bic = fitresult.stats['bic']
                return bic
        #-------------------------------------------------------------------------------

        # Set the weighted penalty function
        self.penaltyfcn = lambda weight,*args: weight*penaltyfcn(*args)
        # Set the selection functional
        self.selectionfcn = selectionfunctional
        self.selection = selection
        # Prepare empty attributes
        self.description = description

        # Get the penalty signature
        if signature is None:
            self.signature =  inspect.getfullargspec(penaltyfcn).args
        else: 
            self.signature = signature

        # Create parameter object for the penalty weight
        newparam = Parameter(parent=self, idx=0, name='weight')
        
        # Add to the penalty object
        setattr(self,'weight',newparam)
        self.weight.lb = np.finfo(float).eps

        # Remove useless attributes
        delattr(self.weight,'par0')
        delattr(self.weight,'linear')
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    def optimize(self,fitfcn,y,sigma):

        if not self.weight.frozen:
            # Extract optimization range from model penalty
            searchrange = np.array([np.log10(self.weight.lb), np.log10(self.weight.ub)])
            searchrange[np.isinf(searchrange)] = 20
            # Construct the selection functional
            selectionFunctional = lambda log10weight: self.selectionfcn(fitfcn,y,sigma,log10weight)

            # Minimization of the selection functional
            log10optweight = fminbound(selectionFunctional,*searchrange, xtol=0.1)

            # Logscale -> linear-scale
            optweight = 10**log10optweight
        else: 
            optweight = self.weight.value

        # Update optimized value to object
        self.optweight = optweight

        # Get the fit result
        fitresult = fitfcn(optweight)

        return fitresult
    #--------------------------------------------------------------------------
#==============================================================================

#--------------------------------------------------------------------------
def _outerOptimization(fitfcn,penalty_objects,y,sigma):

    # If there are no penalties
    if len(penalty_objects)==0:
        fitfcn_ = lambda y: fitfcn(y,[None])

    # Otherwise, prepare to solve multiobjective problem 
    elif len(penalty_objects)==3:
        thirdfcn = lambda y,*param: penalty_objects[2].optimize(lambda weight: fitfcn(y,[*param,weight]),y,sigma)
        secondfcn = lambda y,*param: penalty_objects[1].optimize(lambda weight: fitfcn(y,[*param,weight,thirdfcn(y,*param,weight)]),y,sigma)
        fitfcn_ = lambda y: penalty_objects[0].optimize(lambda weight: fitfcn(y,[weight,secondfcn(y,weight),thirdfcn(y,weight,secondfcn(weight))]),y,sigma)

    elif len(penalty_objects)==2:
        secondfcn = lambda y,*param: penalty_objects[1].optimize(lambda weight: fitfcn(y,[*param,weight]),y,sigma)
        fitfcn_ = lambda y: penalty_objects[0].optimize(lambda weight: fitfcn(y,[weight,secondfcn(y,weight)]),y,sigma) 

    elif len(penalty_objects)==1:
        fitfcn_ = lambda y: penalty_objects[0].optimize(lambda weight: fitfcn(y,[weight]),y,sigma)
    else: 
        raise RuntimeError('The fit() function can only handle up to three penalties.')

    return fitfcn_
#--------------------------------------------------------------------------


#==============================================================================================
def fit(model,y,*constants,par0=None,bootstrap=0, noiselvl=None,**kwargs):
    r"""
    Fit the input model to the data ``y`` via one of the three following approaches: 
    
    - Non-linear least-squares 
    - Regularized linear-least-squares 
    - Separable non-linear least-squares 

    The most appropiate solver is chosen automatically based on the model structure. 

    Parameters
    ----------
    model : :ref:`Model`
        Model object. 
    y : array_like 
        Data to be fitted. 
    par0 : array_like, optional 
        Value at which to initialize the parameter at the start of a fit routine. 
        Must be specified if not defined in the model. Otherwise, it overrides the definition in the model. 

    Any additional solver-specific keyword arguments can be specified. See :ref:`snlls` for a full reference. 

    Returns
    -------
    :ref:`FitResult` with the following fields defined:
    <parameter_name> : :ref:`Parameter`
        Fitted value of the <parameter_name> model parameter.
    <parameter_name>Uncert : :ref:`UQResult`
        Uncertainty quantification of the <parameter_name> model parameter.
    param : ndarray
        Fitted parameter vector ordered according to the model parameter indices.
    modelUncert : :ref:`UQResult`
        Uncertainty quantification of the parameter vector ordered according to the model parameter indices.
    model : ndarray
        Fitted model response.
    modelUncert : :ref:`UQResult`
        Uncertainty quantification of the fitted model response.        
    regparam : scalar
        Regularization parameter value used for the regularization of the linear parameters.
    plot : callable
        Function to display the results. It will display the fitted data.
        The function returns the figure object (``matplotlib.figure.Figure``)
        object as output, which can be modified. A vector for the x-axis and its label can
        be specified by calling ``FitResult.plot(axis=axis,xlabel='xlabel')``.
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
    else:
        model = model.copy()

    if len(constants)>0:
        constants = np.atleast_1d(constants)

    if model.Nlin==0:
        model.addlinear('scale',lb=-np.inf,ub=np.inf,description='Scaling factor')

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

    if not isinstance(y,list):
        y = [y]            
    nprev = 0
    ysubsets = []
    sigmas = []
    for n,yset in enumerate(y):
        ysubsets.append(np.arange(nprev,nprev+len(yset)))
        if noiselvl is None:
            sigmas.append(noiselevel(yset)*np.ones(len(yset)))
        else: 
            sigmas.append(noiselvl[n]*np.ones(len(yset)))
        nprev = nprev+len(yset)
    ysplit = y.copy()
    y = np.concatenate(y)
    sigmas = np.concatenate(sigmas)


    if model.Nlin==0 and model.Nnonlin==0:
        raise AssertionError(f'The model has no parameters to fit.')    

    # Get parameter indices in the order spitted out by the solver
    param_idx = [[]]*len(model._parameter_list('vector'))
    idxprev = 0
    for islinear in [False,True]:
        for n,param in enumerate(model._parameter_list('vector')):
            if np.all(getattr(model,param).linear == islinear):
                N = len(np.atleast_1d(getattr(model,param).idx))
                param_idx[n] = np.arange(idxprev,idxprev + N)
                idxprev += N  

    # If there are penalties in the model
    if np.any([isinstance(getattr(model,attr),Penalty) for attr in dir(model)]):
        penalty_names = [attr for attr in dir(model) if isinstance(getattr(model,attr),Penalty)]
        penalty_objects = [getattr(model,penalty) for penalty in penalty_names]
        # Get the penalty functions
        penalties = [penalty.penaltyfcn for penalty in penalty_objects]
        # Prepare the penalties to input to snlls
        extrapenalties = lambda *weights: [lambda nonlin,lin: penalty(nonlin,lin,weight) for penalty,weight in zip(penalties,weights)]
    else: 
        penalty_objects = []
        # If there are no penalties in the model
        extrapenalties = lambda *_: None

    # Unless specified by the model, disable regularization of linear parameters
    regparam,regparamrange,reg = None,None,False
    # If there is regularization added to the model
    if np.any([isinstance(getattr(model,attr),Regularization) for attr in dir(model)]):
        regularization_object = getattr(model,[attr for attr in dir(model) if isinstance(getattr(model,attr),Regularization)][0])
        reg = True
        if regularization_object.weight.frozen: 
            # Regularization parameter is fixed
            regparam = regularization_object.weight.value 
        else: 
            # Regularization parameter is optimized
            regparam = regularization_object.selection 
            regparamrange = [regularization_object.weight.lb, regularization_object.weight.ub]

    # Prepare the separable non-linear least-squares solver
    Amodel_fcn = lambda param: model.nonlinmodel(*constants,*param)
    fitfcn = lambda y,penweights: snlls(y, Amodel_fcn, par0, lb=lb, ub=ub, lbl=lbl, ubl=ubl, 
                                                subsets=ysubsets, lin_frozen=linfrozen, nonlin_frozen=nonlinfrozen,
                                                regparam=regparam, reg=reg, regparamrange=regparamrange,
                                                extrapenalty=extrapenalties(penweights), **kwargs)        

    # Prepare outer optimization of the penalty weights, if necessary
    fitfcn = _outerOptimization(fitfcn,penalty_objects,y,sigmas)

    # Run the fitting algorithm 
    fitresults = fitfcn(y)

    # If requested, perform a bootstrap analysis
    if bootstrap>0: 
        def bootstrap_fcn(ysim): 
            fit = fitfcn(np.concatenate(ysim))
            if not isinstance(fit.model,list): fit.model = [fit.model]
            return (fit.param,*fit.model)
        # Bootstrapped uncertainty quantification
        param_uq = bootan(bootstrap_fcn,ysplit,fitresults.model,samples=bootstrap)
        # Include information on the boundaries for better uncertainty estimates
        paramlb = model._vecsort(model._getvector('lb'))[np.concatenate(param_idx)] 
        paramub = model._vecsort(model._getvector('ub'))[np.concatenate(param_idx)] 
        fitresults.paramUncert = UQResult('bootstrap',data=param_uq[0].samples,lb=paramlb,ub=paramub)
        fitresults.param = fitresults.paramUncert.median
        # Get the uncertainty estimates for the model response
        fitresults.model = [param_uq[n].median for n in range(1,len(param_uq))]
        fitresults.modelUncert = [param_uq[n] for n in range(1,len(param_uq))]
        if len(fitresults.model)==1: 
            fitresults.model = fitresults.model[0]
            fitresults.modelUncert = fitresults.modelUncert[0]
    # Get some basic information on the parameter vector
    keys = model._parameter_list(order='vector')
 
    # Dictionary of parameter names and fitted values
    FitResult_param = {key : fitvalue if len(fitvalue)>1 else fitvalue[0] for key,fitvalue in zip(keys,[fitresults.param[idx] for idx in param_idx])}
    # Dictionary of parameter names and fit uncertainties
    FitResult_paramuq = {f'{key}Uncert': model._getparamuq(fitresults.paramUncert,idx) for key,idx in zip(keys,param_idx)}
    # Dictionary of other fit quantities of interest
    FitResult_dict = {key: getattr(fitresults,key) for key in ['param','paramUncert','model','modelUncert','cost','plot','residuals','stats','regparam']}

    _paramlist = model._parameter_list('vector')
    def propagate(model,*constants,lb=None,ub=None):
    # ----------------------------------------------------------------------------
        # Get the parameter names of the input model
        if isinstance(model,Model):
            modelparam = model._parameter_list('vector')
        elif callable(model):
            modelparam = inspect.getfullargspec(model).args
        else: 
            raise TypeError('The input must be a deerlab.Model object or a callable.')

        # Check that all parameters are in the fit object
        for param in modelparam:
            if not hasattr(fit,param): 
                raise KeyError(f'The fit object does not contain the {param} parameter.')
        # Determine the indices of the subset of parameters the model depends on
        subset = [param_idx[np.where(np.asarray(_paramlist)==param)[0][0]] for param in modelparam]
        subset = np.concatenate(subset)   
        # Propagate the uncertainty from that subset to the model
        modeluq = fitresults.paramUncert.propagate(lambda param: model(*constants,*param[subset]),lb,ub)
        return modeluq
    # ----------------------------------------------------------------------------
    def evaluate(model,*constants):
    # ----------------------------------------------------------------------------
        # Get the parameter names of the input model
        modelparam = model._parameter_list('vector')

        # Check that all parameters are in the fit object
        for param in modelparam:
            if not hasattr(fit,param): 
                raise KeyError(f'The fit object does not contain the {param} parameter.')
        # Determine the indices of the subset of parameters the model depends on
        subset = [np.where(np.asarray(_paramlist)==param)[0][0] for param in modelparam]
        # Evaluate the input model
        modeluq = model(*constants,*fitresults.param[subset])
        return modeluq
    # ----------------------------------------------------------------------------


    # Generate FitResult object from all the dictionaries
    fit = FitResult({**FitResult_param,**FitResult_paramuq, **FitResult_dict, 'propagate': propagate, 'evaluate': evaluate}) 

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
    try:
        np.testing.assert_equal(a,b)
        return True
    except Exception:
        return False
# ==============================================================================
def link(model,**links):
    """
    Link model parameters. 

    Parameters
    ----------
    model : :ref:`Model`
        Model object. 
    links : keyword-argument pairs 
        Keyword-argument pairs, where the arguments must be lists of model parameter names. 
        The corresponding model parameter will be assigned to new parameter whose name is given 
        by the keyword name. For example:: 

            newmodel = link(model,parC = ['parA','parB'])

        will return a new model where the values of ``parA`` and ``parB`` will be given by the
        new model parameter ``parC``. 

    Returns
    -------
    newmodel : :ref:`Model`
        New model object without the linked parameter and with the newly defined parameters. 
    """
    def _linkparameter(model,parameters,newname):
    # ---------------------------------------------------------------------  
        # Get a list of parameter names in the model
        model_parameters = model._parameter_list(order='vector')

        link_parameters,link_indices,link_param_idx = [],[],[]
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
                        break
            else: 
                link_names.append(param)

        # Remove the first from the list as it will be kept
        linked_name = link_names.pop(0)
        
        for n,name in enumerate(model_parameters):
            if name==linked_name: link_param_idx = n

        

        # Get the vector index of the parameter to be linked to
        link_indices = np.atleast_1d(link_indices[0])

        nnew = 0
        # Initialize the maps linked->unlinked
        mapping = np.zeros(model.Nnonlin,dtype=int)
        mapping_linear = np.zeros(model.Nlin,dtype=int)

        # Get the indices of the unlinked parameters in the maps
        unlinked_linear_idx = np.full(len(model_parameters),None)
        unlinked_nonlinear_idx = np.full(len(model_parameters),None)
        linked_linear_idx = np.full(len(model_parameters),None)
        linked_nonlinear_idx = np.full(len(model_parameters),None)
        q = 0
        nnew = 0
        for n,param in enumerate(model_parameters):
            if np.all(getattr(model,param).linear):
                m =  len(np.atleast_1d(getattr(model,param).idx)) 
                unlinked_linear_idx[n]= np.arange(q,q+m)
                q += m
                if param not in link_names:
                    linked_linear_idx[n]= np.arange(nnew,nnew+m)
            else:
                unlinked_nonlinear_idx[n] = np.array(n)
                if param not in link_names: 
                    linked_nonlinear_idx[n] = np.array(nnew)
                    m = 1
            if param not in link_names: 
                nnew += m
        Nnl = model.Nnonlin
        
        # Loop over all parameters in the model
        for n,param in enumerate(model_parameters):
            # If the parameter is to be linked...
            if param in link_names:
                # Update the number of parameters in the model
                Nremoved = len(np.atleast_1d(getattr(model,param).idx))
                model.Nparam -= Nremoved
                if np.all(getattr(model,param).linear):
                    model.Nlin -= Nremoved
                    # Update the parameter vector map
                    mapping_linear[unlinked_linear_idx[n]] = link_indices-Nnl                
                else: 
                    model.Nnonlin -= Nremoved
                    # Update the parameter vector map
                    mapping[unlinked_nonlinear_idx[n]] = link_indices
                # Delete the linked parameter from the model
                delattr(model,param)

            # Otherwise if the parameter is not linked...
            else:
                # Update the index of the parameter in the new vector 
                if not np.all(getattr(model,param).linear):
                    getattr(model,param).idx = linked_nonlinear_idx[n]
                    mapping[unlinked_nonlinear_idx[n]] = linked_nonlinear_idx[n]
                else: 
                    getattr(model,param).idx = linked_linear_idx[n]
                    mapping_linear[unlinked_linear_idx[n]] = linked_linear_idx[n]-Nnl

        # Delete the old copy with the old name
        paramobj = getattr(model,model_parameters[link_param_idx])
        delattr(model,model_parameters[link_param_idx])
        # Create a copy of the linked parameter with the new name
        setattr(model,newname,paramobj)

        # Monkey-patch the evaluation function         
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
    for link_newname in links: 
        to_link = [getattr(newmodel,parname) for parname in links[link_newname]]
        newmodel = _linkparameter(newmodel,to_link,link_newname)
    return newmodel
#==============================================================================================

# ---------------------------------------------------------------------
def _unique_ordered(vec):
    uniques = []
    for v in vec: 
        if v not in uniques:
            uniques.append(v)
    return uniques
# ---------------------------------------------------------------------


#==============================================================================================
def _combinemodels(mode,*inputmodels,addweights=False): 

    # Initialize empty containers
    subsets_nonlin,arguments,arelinear = [],[],[]
    nprev = 0

    if len(inputmodels)==1:
        return inputmodels[0]

    # Make deep-copies of the models to avoid modifying them
    models = [deepcopy(model) for model in inputmodels]

    if addweights:

        for n,(model,nonlinfcn) in enumerate(zip(models,[model.nonlinmodel for model in models])):
            signature = model._nonlinsignature
            signature.append('weight')
            def make_weighted_comb(nonlinfcn):
                def weighted_comb(*inputargs):
                    weight = inputargs[-1]
                    param = inputargs[:-1]
                    return weight*nonlinfcn(*param)
                return weighted_comb
            constants = [entry['argkey'] for entry in model._constantsInfo]
            weightedModel = Model(make_weighted_comb(nonlinfcn),constants=constants,signature=signature)
            for name in model._parameter_list(order='vector'):
                if np.any(getattr(model,name).linear):
                    weightedModel.addlinear(name,vec=len(np.atleast_1d(getattr(model,name).idx)))                    
                getattr(weightedModel,name).set(**_importparameter(getattr(model,name)))
            getattr(weightedModel,'weight').set(lb=0,par0=1,description='Weighting factor')
            models[n] = deepcopy(weightedModel)

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

        oldargs = models[n]._parameter_list(order='vector')
        i = 0
        for oldkey,newkey in zip(oldargs,newarguments): 
            if isinstance(getattr(model,oldkey).idx,np.ndarray):
                newarguments = np.insert(newarguments,i*np.ones(len(np.atleast_1d(getattr(model,oldkey).idx))-1,dtype=int),newkey).tolist()
            i += len(np.atleast_1d(getattr(model,oldkey).idx))

        # Add the submodel arguments to the combined model signature
        arguments += newarguments

    # Preparation of the combined model signature    
    arelinear = np.asarray(arelinear).astype(bool)
    arguments = np.array(arguments)
    lin_params = arguments[arelinear]
    nonlin_params = arguments[~arelinear]

    # Account for the constant arguments
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
    nonlinfcns = [model.nonlinmodel for model in models]
    Nlins = [model.Nlin for model in models]
    ysizes = [[]]*len(models)

    #---------------------------------------------------------------------
    def _combined_nonlinmodel(*inputargs):
        """Evaluates the nonlinear functions of the submodels and 
        concatenates them into a single design matrix"""

        nonlocal ysizes

        constants = inputargs[:Nconst]
        param = inputargs[Nconst:]

        param = np.atleast_1d(param)
        constants = np.atleast_2d(constants)
        # Loop over the submodels in the model
        Amatrices = []
        for n,nonlinfcn in enumerate(nonlinfcns):
            # Evaluate the submodel
            Amatrix = np.atleast_2d(nonlinfcn(*constants[const_subsets[n],:],*param[subsets_nonlin[n]]))
            if np.shape(Amatrix)[1]!=Nlins[n]: Amatrix = Amatrix.T
            Amatrices.append(Amatrix)
        if mode=='expand':
            ysizes = [A.shape[0] for A in Amatrices]
            Anonlin_full = block_diag(Amatrices).toarray()
            if not any(arelinear):
                Anonlin_full = np.sum(Anonlin_full,1)

        elif mode=='combine':
            Anonlin_full = np.hstack(Amatrices)

        return Anonlin_full
    #---------------------------------------------------------------------
    
    #---------------------------------------------------------------------
    def _split_output(y,*inputargs):
        nonlocal ysizes
        nprev = 0
        ysubsets = []
        for x in ysizes:
            ysubsets.append(np.arange(nprev,nprev+x))
            nprev = nprev+x
        return [y[ysubsets[n]] for n in range(len(ysizes))]
    #---------------------------------------------------------------------

    # Create the model object
    combinedModel = Model(_combined_nonlinmodel,constants=constants,signature=signature)

    # Add parent models 
    combinedModel.parents = models

    if mode=='expand':
        # Add post-evalution function for splitting of the call outputs
        setattr(combinedModel,'_posteval_fcn',_split_output) 

    # Add the linear parameters from the subset models   
    lin_param_set = []
    for param in _unique_ordered(lin_params):
        lin_param_set.append({'name':param, 'vec':np.sum(lin_params==param)})

    for lparam in lin_param_set:
        combinedModel.addlinear(lparam['name'], vec=lparam['vec'])

    parameters = np.concatenate([arguments,lin_params])
    # Import all parameter information from the subset models
    for name,param in zip(combinedModel._parameter_list(order='vector'),parameters):
        if '_' in name:
            param = name.rsplit('_',1)[0]
            n = name.rsplit('_',1)[-1]
            n = int(n)-1
        else:
            param,n = name,0
        getattr(combinedModel,name).set(**_importparameter(getattr(models[n],param)))

    # Return the new combined model object
    return combinedModel    
#==============================================================================================

#==============================================================================================
def expand(*inputmodels,addweights=False):
    """
    Create a multi-response model from multiple individual models. 

    Parameters
    ----------
    inputmodels : :ref:`Model` objects
        Model objects to be combined. If one of the models has no linear parameters, a linear 
        scaling factor parameters will be added. The names of the ``N``-th input model parameter will be 
        changed by a suffix ``_N`` in the new model. Example:: 

            newmodel = expand(model1,model2)
            newmodel.parA_1 # Originally parA from model1
            newmodel.parA_2 # Originally parA from model2


    addweights : boolean, optional 
        If true, the function will add a non-linear weight parameter for each model response
        even if the individual models have linear parameters. 

    Returns
    -------
    newmodel : :ref:`Model`
        New model object taking the combined parameter set and returning a list of model reponses
        correponding to each of the input models. 
    """
    return _combinemodels('expand',*inputmodels,addweights=addweights)
#==============================================================================================

#==============================================================================================
def combine(*inputmodels,addweights=False):
    """
    Create model whose response is a linear combination of multiple individual model responses. 

    Parameters
    ----------
    inputmodels : :ref:`Model` objects
        Model objects whose linear responses are to be linearly combined. If one of the models 
        has no linear parameters, a linear scaling factor parameters will be added. The names 
        of the ``N``-th input model parameter will be changed by a suffix ``_N`` in the new model. Example:: 

            newmodel = expand(model1,model2)
            newmodel.parA_1 # Originally parA from model1
            newmodel.parA_2 # Originally parA from model2

    addweights : boolean, optional 
        If true, the function will add a non-linear weight parameter for each model response
        even if the individual models have linear parameters. 

    Returns
    -------
    newmodel : :ref:`Model`
        New model object taking the combined parameter set and returning a response that is a linear
        combination of the input models.
    """
    return _combinemodels('combine',*inputmodels,addweights=addweights)
#==============================================================================================


#==============================================================================================
def relate(model,**functions):
    """
    Create functional relationships between model parameters. 

    Parameters
    ----------
    model : :ref:`Model`
        Model object. 
    functions : keyword-callable pairs 
        Functions describing the relationship between parameters. The keyword represents the parameter
        which will be funtionalized. The keyword argument must be a callable function taking a number 
        of model parameters (names must match any of the model parameter names) as input and returning
        the value to be assigned to the functionalized parameter. For example::

            newmodel = relate(model, parA = lambda parB: 2*parB)

        will create a new model ``newmodel`` based on ``model`` where the parameter ``parA`` is now 
        given by twice the value of parameter ``parB``. The model ``newmodel`` will no longer have ``parA`` 
        as a model parameter and will have a parameter less than ``model``. 
        Multiple parameters can be functionalized by specifying multiple keyword-callable pairs.

    Returns
    -------
    newmodel : :ref:`Model`
        New model object without the functionalized parameters. 
    """
    def _relate(model,function,dependent_name):
    # ---------------------------------------------------------------------  
        
        # Get a list of parameter names in the model
        model_parameters =np.array(model._parameter_list(order='vector'))

        # Get the index of the parameter which will be made dependent
        dependent_idx = np.where(model_parameters==dependent_name)

        # Get the names and indices of the parameters taken by the dependent's function
        arguments_names = inspect.getfullargspec(function).args

        if dependent_name not in model_parameters:
            raise KeyError(f"The assigned parameter '{dependent_name}' is not a parameter of the input model.")

        if getattr(model,dependent_name).linear:
            raise TypeError(f"Linear parameters cannot be used.")

        for arg in arguments_names:
            if arg not in model_parameters:
                raise KeyError(f"The function argument '{arg}' is not a parameter of the input model.")
            if getattr(model,arg).linear:
                raise TypeError(f"Linear parameters cannot be used.")

        param_idx = 0
        # Loop over all parameters in the model
        for param in model_parameters:
            Nidx = len(np.atleast_1d(getattr(model,param).idx))
            # Update the index of the parameter in the new vector 
            if not np.all(getattr(model,param).linear):
                getattr(model,param).idx = np.array(param_idx) 
            else: 
                getattr(model,param).idx = np.arange(param_idx,param_idx+Nidx)

            if param != dependent_name:
                param_idx += Nidx

        # If the parameter is to be linked...
        Nremove = len(np.atleast_1d(getattr(model,dependent_name).idx))
        # Update the number of parameters in the model
        model.Nparam -= Nremove
        if np.all(getattr(model,dependent_name).linear):
            model.Nlin -= Nremove
        else: 
            model.Nnonlin -= Nremove
        # Delete the linked parameter from the model
        delattr(model,dependent_name)

        # Monkey-patch the evaluation function         
        nonlinfcn = model.nonlinmodel
        Nconstants = len(model._constantsInfo)
        nonlinparams = np.array([param for param in model._parameter_list('vector') if not np.all(getattr(model,param).linear)])
        arguments_idx = np.concatenate([np.where(nonlinparams==name)[0] for name in arguments_names])
        dependent_idx = dependent_idx[0]
        # ---------------------------------------------------------------------
        def dependency_model_with_constants(*inputargs):
            # Redistribute the input parameter vector according to the mapping vector
            constants = inputargs[:Nconstants]
            θ = np.atleast_1d(inputargs[Nconstants:]).astype(float)
            θ = np.insert(θ,dependent_idx,function(*θ[arguments_idx]))  
            args = list(θ)
            if model._constantsInfo is not None:
                for info,constant in zip(model._constantsInfo,constants):
                    args.insert(info['argidx'],constant)                
            A = nonlinfcn(*args)
            return A
        # ---------------------------------------------------------------------
        model.nonlinmodel = dependency_model_with_constants

        # Return the updated model with the linked parameters
        return model
    # ---------------------------------------------------------------------

    if not isinstance(model,Model):
        raise TypeError('The first argument must be a deerlab.Model object.')
    newmodel = deepcopy(model)

    # Get the dependent's names and their function arguments
    dependents = [dependent for dependent in functions]
    arguments = [inspect.getfullargspec(functions[dependent]).args for dependent in dependents]

    # Check and correct the order to of functionalization
    maxtrials = 2*len(dependents)
    # Run for a maximum number of trials
    for _ in range(maxtrials):
        # Loop over all dependent variables... 
        for n,dependent in enumerate(dependents):
            # ...and check wheter there is a conflict with the other arguments
            for m in range(n,len(arguments)):
                # If a dependent is to be used as an argument later, there is a conflict
                conflict = dependent in arguments[m]
                if conflict: break
            if conflict: break
        else: 
            # If there are no conflicts, proceed with the functionalization
            break
        # If there is a conflict, swap the order of the dependents and check again
        dependents[n],dependents[m] = dependents[m],dependents[n] 
        arguments[n],arguments[m] = arguments[m],arguments[n] 
    else: 
        # If no solution could be found (due to cyclic relations), raise an error
        raise RuntimeError('There are cyclic relationships in the parameter definitions that could not be resolved.')

    # Loop over all parameter functionalizations...
    for dependent in dependents: 
        # ...and functionalize them one by one
        newmodel = _relate(newmodel,functions[dependent],dependent)

    return newmodel
#==============================================================================================
