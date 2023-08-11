# model.py - DeerLab's modeling interface
# ---------------------------------------------------------------------------
# This file is a part of DeerLab. License is MIT (see LICENSE.md). 
# Copyright(c) 2019-2023: Luis Fabregas, Stefan Stoll and other contributors.

import numpy as np
from scipy.sparse import block_diag
from scipy.optimize import fminbound
from deerlab.classes import UQResult
from deerlab.bootstrap_analysis import bootstrap_analysis
from deerlab.utils import formatted_table, parse_multidatasets
import inspect 
from copy import copy,deepcopy
from functools import partial
import difflib
#===================================================================================
class Parameter(): 
    r"""
    Represents a model parameter or a single parameter vector. 

    Attributes
    ----------
    name : string 
        Name of the parameter

    description : string 
        Description of the parameter

    unit : string 
        Physical unit of the parameter

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



    """

    #=======================================================================================
    #                                         Constructor
    #=======================================================================================

    #---------------------------------------------------------------------------------------
    def __init__(self, name=None, parent=None, idx=None, description=None, par0=None, frozen=False, lb=-np.inf, ub=np.inf,value=None, unit=None, linear=False): 
        """
        Construct a new parameter object.

        Parameters
        ----------
        name : str, optional
            Name of the parameter.
        
        parent : object, optional
            Reference to the parent object that contains the parameter.
        
        idx : int or list of int, optional
            Indices of the parameter in the parent object.
        
        description : str, optional
            Description of the parameter.
        
        par0 : float or array_like, optional
            Value at which to initialize the parameter at the start of a fit routine. 
            Must be specified in the model or latest upon fitting.
        
        frozen : bool, optional
            Whether the parameter will be frozen at a specific value during a fit.
            Default is False.
        
        lb : float or array_like, optional
            Lower bound of the parameter. If not specified, it is assumed to be unbounded.
            Default is -inf.
        
        ub : float or array_like, optional
            Upper bound of the parameter. If not specified, it is assumed to be unbounded.
            Default is inf.
        
        value : float, optional
            Value at which the parameter will be frozen during a fit.
        
        unit : str, optional
            Physical unit of the parameter.
        
        linear : bool, optional
            Whether the model behaves linearly with respect to the parameter.
            Default is False.
        """
        # Attributes
        self.name = name      # Name of parameter
        self._parent = parent # Parent 
        self.idx = idx        # Index of parameter
        self.description = description # Description
        self.unit = unit      # Unit
        self.par0 = par0      # Start values
        self.lb = lb          # Lower bounds
        self.ub = ub          # Upper bounds
        self.frozen = frozen  # Frozen
        self.value = value    # Value at which it is frozen
        self.linear = linear  # Linearity
    #---------------------------------------------------------------------------------------

    #=======================================================================================
    #                                         Methods
    #=======================================================================================

    #---------------------------------------------------------------------------------------
    def set(self,**attributes):
        """
        Set one or multiple attributes for a parameter. See list of attributes for a reference list. 

        This method allows users to set the value of one or more attributes for a ``Parameter`` object. 
        The list of attributes that can be set are the same as the attributes of the ``Parameter`` class, 
        namely: ``name``, ``description``, ``unit``, ``par0``, ``lb``, ``ub``, ``linear``, ``frozen``, and ``value``.

        Parameters
        ----------
        attributes : keyword/values pairs
            Pairs of keywords defining the parameter attribute and the value to be assigned.

        Examples
        --------
        Setting a parameter's start value ::

            parameter.set(par0=0.5)

        Setting both the lower bound and upper bound values :: 

            parameter.set(lb=0, ub=1)

        """
        for key in attributes:
            if not hasattr(self,key):
                raise AttributeError(f"'{key}' is not a valid parameter attribute.")
            setattr(self,key,attributes[key])
        return 
    #---------------------------------------------------------------------------------------
    
    #---------------------------------------------------------------------------------------
    def setas(self,parameter):
        """
        Copy all attributes from an input parameter to the current parameter. 

        This method allows the user to copy all attributes from an input parameter object
        to the current parameter object. This can be useful when initializing multiple
        parameters with the same values.

        Parameters  
        ----------
        parameter : ``Parameter`` object 
            Model parameters from which to extract the attributes. 
        """
        self.set(**_importparameter(parameter))
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def freeze(self,value):
        """
        Freeze a parameter during a fit/optimization to a given value.

        This method allows the user to freeze a parameter during optimization to a
        specific value. The parameter will remain fixed at this value during the
        optimization process, and will not be updated. This method does not affect
        model evaluation.

        Parameters
        ----------
        value : float or array_like
            Value at which to freeze the parameter during optimization. This value
            must be within the lower and upper bounds of the parameter.
        """

        if np.any(value>self.ub) or np.any(value<self.lb):
            if len(np.atleast_1d(value))>1:
                raise ValueError(f"Frozen values are outside of the bounds.")
            else: 
                raise ValueError(f"Frozen value {value} is outside of the bounds {self.lb} and {self.ub}.")

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
        Release a frozen parameter's value during a fit/optimization.

        This method allows the user to release a previously frozen parameter during
        optimization. The parameter will no longer be fixed at a specific value and
        will be updated according to the optimization algorithm. This method does not
        affect model evaluation.
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
    r"""
    Represents a model.

    The ``Model`` class provides a way to create a new model object from a
    non-linear function. The ``Model`` class also provides access to the model
    parameters and allows them to be modified.

    Attributes
    ----------
    <parameter_name> : :ref:`Parameter` 
        Model parameter. One :ref:`Parameter` instance is assigned for each
        parameter (with name ``<parameter_name>``) in the model.  
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
    """


    #=======================================================================================
    #                                         Constructor
    #=======================================================================================
    
    # Use a wrapper function to facilitate internal arguments manipulation        
    #---------------------------------------------------------------------------------------
    def _model_with_constants(self,nonlinfcn,constantsInfo,*inputargs):
        """
        Evaluate a non-linear model with constant arguments.

        This internal method is used to evaluate a non-linear model with constant
        arguments. The method takes the non-linear function, a list of constant
        arguments, and the input arguments and returns the output of the non-linear
        model with the constant arguments inserted at the appropriate positions.

        Parameters
        ----------
        nonlinfcn : callable
            Non-linear function to be evaluated.
        constantsInfo : list
            List of dictionaries containing information about the constant arguments
            to be used with the non-linear function. Each dictionary must contain
            two keys: "argkey" and "argidx". "argkey" must be the name of the
            argument, and "argidx" must be the index of the argument in the function
            signature.
        inputargs : variable length argument list
            Input arguments to be passed to the non-linear function.

        Returns
        -------
        output : variable type
            Output of the non-linear function with the constant arguments inserted
            at the appropriate positions.
        """
        # Get the number of constants
        Nconstants = len(constantsInfo)
        # Separate the constants from the other arguments
        constants = inputargs[:Nconstants]
        θ = inputargs[Nconstants:]
        # Create a list of arguments for the nonlinear function
        args = list(θ)
        # If there are constants, insert them into the list of arguments at the appropriate positions
        if constantsInfo is not None:
            for info,constant in zip(constantsInfo,constants):
                args.insert(info['argidx'],constant)
        # Call the nonlinear function with the arguments and return the result
        return nonlinfcn(*args)
    #---------------------------------------------------------------------------------------
    
    #---------------------------------------------------------------------------------------
    def __init__(self,nonlinfcn,constants=None,signature=None): 
        """
        Construct a new model from a non-linear function. 

        Parameters
        ----------

        nonlinfcn : callable 
            Function that takes a set of non-linear parameters and 
            returns either the full model response or the design matrix 
            of the model response. A parameter will be added to the new model for each 
            input argument defined in the function signature.   
        
        constants : string or list thereof, optional
            Names of the arguments taken by the ``nonlinfcn`` function to be defined as
            constants. These will not be added as parameters to the new model.

        signature : list of strings, optional
            Signature of the ``nonlinfcn`` function to manually specify the names
            of the input arguments. For internal use (mostly).

        Returns
        -------

        model : ``Model`` object 
            Model object instance that takes the parameters defined for ``nonlinfcn``
            and returns the output of ``nonlinfcn``.

        """
        # Check if nonlinfcn is callable and create a wrapper function if it's not
        if not callable(nonlinfcn):
            Amatrix = nonlinfcn.copy()
            nonlinfcn = lambda *_: Amatrix 
        self.nonlinmodel = nonlinfcn
        self.description = None
        self._constantsInfo = []
        self.parents = None

        # Check if signature is specified, and if not use the function's signature
        if signature is None:
            signature = inspect.getfullargspec(nonlinfcn).args

        # Create a copy of the signature for future reference
        self._nonlinsignature = signature.copy()
        parameters = signature.copy()
        self.signature = signature

        # Check if constants are specified and remove them from the parameters list
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

        # Construct the callable non-linear part of the model
        self.nonlinmodel = partial(self._model_with_constants,nonlinfcn,self._constantsInfo)

        # Update the number of parameters in the model
        self.Nparam = len(parameters)
        self.Nnonlin = len(parameters)
        self.Nlin = 0

        # Construct the parameters and add them to the model
        for n,param in enumerate(parameters): 
            newparam = Parameter(parent=self, idx=n, name=param)
            setattr(self,param,newparam)
    #---------------------------------------------------------------------------------------

    # Gets called when an attribute is accessed
    #--------------------------------------------------------------------------------
    def __getattribute__(self, attr):
        """
        This method is called when an attribute of this class instance is accessed.
        It first tries to access the attribute using the ``super()`` method, which
        will go up the inheritance chain to look for the attribute. If the attribute
        is not found, it checks for similar attributes that exist in the class 
        instance, and suggests those as potential matches in the error message. 
        If no matches are found, it raises an ``AttributeError`` with a custom message
        indicating that the requested attribute does not exist.
        
        Parameters
        ----------
        attr : str
            The name of the attribute being accessed.
        
        Raises
        ------
        AttributeError
            If the requested attribute does not exist in the class instance.
        """
        # First try to access the attribute as usual
        try:
            return super(Model, self).__getattribute__(attr)
        
        # If an AttributeError is raised, handle it by printing a custom error message
        except AttributeError:
            # Create the error message with the name of the attribute that was accessed
            errstr = f"The model has no attribute '{attr}'."
            # Get a list of all the attributes in the model
            attributes = [key for key in self.__dict__]
            # Use difflib to get a list of similar attributes to the one that was accessed
            proposal = difflib.get_close_matches(attr, attributes)
            # If there are any similar attributes, add a suggestion to the error message
            if len(proposal)>0:
                errstr += f' \n Did you mean: {proposal} ?'
            # Raise the AttributeError with the custom error message
            raise AttributeError(errstr)
    #--------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def _parameter_list(self, order='alphabetical'):
        """
        (Private method)
        
        Get the list of parameters defined in the model sorted alphabetically or by vector definition.

        Parameters
        ----------
        order : str, optional
            Specifies how the parameters should be sorted. Can be either ``'alphabetical'`` or ``'vector'``.
            Default is ``'alphabetical'``.

        Returns
        -------
        keylist : list
            List of the parameters in the model sorted according to the ``order`` parameter.
        """
        keylist = [param for param in dir(self) if isinstance(getattr(self,param),Parameter)]
        # Sort the list alphabetically if specified
        if order=='alphabetical':
            pass
        # Sort the list by vector definition if specified
        elif order=='vector':
            # If there are any parameters in vector form...
            n = 0
            for key in keylist: 
                if isinstance(getattr(self,key).idx,np.ndarray):
                    # ...insert the key string for the same number of linear parameters in that vector 
                    keylist = np.insert(keylist,n*np.ones(len(np.atleast_1d(getattr(self,key).idx))-1,dtype=int),key)  
                n += len(np.atleast_1d(getattr(self,key).idx))
            keylist = self._vecsort(keylist)
        else: 
            raise KeyError("The 'order' argument muse be either 'alphabetical' or 'vector'.")
        # Remove any duplicates
        keylist = list(dict.fromkeys(keylist))
        return keylist
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def _vecsort(self, paramlist):
        """
        (Private method)

        Sort vectorized parameters attributes from alphabetical ordering to vector indexing.

        This internal method sorts the list of parameters defined in the model according to their index in the parameter vector.

        Parameters
        ----------
        paramlist : list
            rameter names to be sorted.
        
        Returns
        -------
        orderedlist : ndarray
            List of parameter names sorted according to their index in the parameter vector.
        """
        # Convert the input list to a numpy array and squeeze any singleton dimensions
        paramlist = np.squeeze(np.atleast_1d(paramlist))

        # Get a list of indices for all parameters in the model
        indices = np.concatenate([np.atleast_1d(getattr(self,param).idx) for param in dir(self) if isinstance(getattr(self,param),Parameter)])

        # Create a copy of the input list and sort it according to the indices
        orderedlist = np.atleast_1d(paramlist.copy())
        orderedlist[indices] = paramlist

        return orderedlist
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def _split_linear(self,variable):
        """
        (Private method)

        Split a vector into non-linear and linear parameter subsets.

        This method is used to split a vector of model parameters into two subsets:
        the subset of non-linear parameters and the subset of linear parameters.

        Parameters
        ----------
        variable : ndarray
            Vector of model parameters to split into non-linear and linear subsets.

        Returns
        -------
        variable_lin : ndarray
            Subset of `variable` that contains only the linear parameters.
        variable_nonlin : ndarray
            Subset of `variable` that contains only the non-linear parameters.

        """
        # Convert `variable` to a numpy array with at least one dimension
        variable = np.atleast_1d(variable)

        # Get a list of all model parameters that are linear
        linear = np.concatenate([np.atleast_1d(getattr(self,param).linear) for param in dir(self) if isinstance(getattr(self,param),Parameter)]) 

        # Sort the list of linear parameters
        linear = self._vecsort(linear)

        # Split `variable` into a non-linear subset and a linear subset
        variable_nonlin = variable[~linear]
        variable_lin = variable[linear]

        return variable_lin, variable_nonlin
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def _merge_linear(self,variable_nonlin,variable_lin):
        """
        (Private method)

        Merge a vector's non-linear and linear parameter subset vectors.

        Parameters
        ----------
        variable_nonlin : array-like
            Non-linear part of the variable vector.
        variable_lin : array-like
            Linear part of the variable vector.

        Returns
        -------
        variable : ndarray
            Merged variable vector.
        """
        variable = np.zeros(len(variable_nonlin)+len(variable_lin))
        # Get the list of linear parameters in the model
        linear = np.concatenate([np.atleast_1d(getattr(self,param).linear) for param in dir(self) if isinstance(getattr(self,param),Parameter)])
        # Sort the list of linear parameters in the model
        linear = self._vecsort(linear)
        # Merge the non-linear and linear parts of the variable vector
        variable[~linear] = variable_nonlin
        variable[linear] = variable_lin
        return variable
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def _getvector(self,attribute):
        """
        (Private method)

        Get a vector of the given attribute for all parameters defined in the model. 

        Parameters
        ----------

        attribute : string
            Attribute to get for each parameter in the model.

        Returns
        -------

        vector : ndarray
            Vector of the specified attribute for all parameters in the model.

        """
        return np.concatenate([np.atleast_1d(getattr(getattr(self,param),attribute)) for param in dir(self) if isinstance(getattr(self,param),Parameter)])    
    #---------------------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    def _getparamuq(self,uq_full,paramidx):
        """
        (Private method)

        Get the uncertainty quantification of an individual parameter.

        Parameters
        ----------
        paramidx : int
            Index of the parameter to get the uncertainty quantification for.
        uq_full : UQResult
            Result of a full uncertainty quantification on the model.

        Returns
        -------
        param_uq : UQResult
            Result of the uncertainty quantification for the selected parameter.
        """
        # Create a subset model that only includes the selected parameter
        subset_model = lambda x: x[paramidx]

        # Get the bounds and freezing of the selected parameter
        param_lb =  self._vecsort(self._getvector('lb'))[paramidx]
        param_ub =  self._vecsort(self._getvector('ub'))[paramidx]
        frozen = self._vecsort(self._getvector('frozen'))[paramidx]

        # If the parameter is frozen, return a void UQResult
        if np.all(frozen): 
            param_uq = UQResult('void')
        else:
            # Propagate the full UQ through the subset model
            if uq_full.type=='moment':
                param_uq = uq_full.propagate(subset_model,lb=param_lb, ub=param_ub)
            elif uq_full.type=='bootstrap':
                param_uq = UQResult('bootstrap',data=uq_full.samples[:,paramidx],lb=param_lb, ub=param_ub)
            else:
                param_uq = UQResult('void')
        return param_uq
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    def _error_if_already_exists(self, key):
        """
        (Private method)

        Check if the model has a parameter with a given name.

        If a parameter with the given name exists in the model, raise a KeyError.

        Parameters
        ----------
        key : str
            Name of the parameter to check for.

        Raises
        ------
        KeyError
            If a parameter with the given name exists in the model.
        """
        if hasattr(self, key):
            raise KeyError(f'The model already has a "{key}" parameter.')
    #-----------------------------------------------------------------------------
    
    
    #=======================================================================================
    #                                         Methods
    #=======================================================================================

    #---------------------------------------------------------------------------------------
    def addnonlinear(self, key, lb=-np.inf, ub=np.inf, par0=None, name=None, unit=None, description=None):
        """
        Add a new non-linear parameter (:ref:`Parameter` object) to the model. 

        Parameters
        ----------
        key : string
            Identifier of the parameter. This name will be used to refer
            to the parameter in the model.

        lb : float or array_like, optional
            Lower bound of the parameter. If not specified, it is set to ``-np.inf``.

        ub : float or array_like, optional
            Lower bound of the parameter. If not specified, it is set to ``+np.inf``.

        description : string, optional 
            Description of the parameter. 

        unit : string, optional
            Physical unit of the parameter.
        """

        # Check that the parameter does not exist
        self._error_if_already_exists(key)

        # Update number of parameters
        idx = self.Nparam
        self.Nparam += 1
        self.Nnonlin += 1

        # Construct the new parameter object and added to the model
        newparam = Parameter(name=key, linear=False, parent=self, idx=idx, par0=par0, lb=lb, ub=ub, unit=unit, description=description)
        setattr(self,key,newparam)
        Nconstants = len(self._constantsInfo)
        Amodel = self.nonlinmodel
        topop = self.Nnonlin-1

        # Define a new model function with added non-linear parameter
        #------------------------------------------------
        def model_with_constants_and_added_nonlin(*inputargs):
            """
            A helper function that takes in a set of input arguments, 
            extracts constants and non-linear parameters from the input
            arguments and returns the model response.
            """
            # Extract constants from the input arguments
            constants = inputargs[:Nconstants]
            # Extract non-linear parameters from the input arguments
            θ = inputargs[Nconstants:]
            # Create a list of all arguments
            args = list(θ)
            args.pop(topop)
            # If there are any constant parameters
            if self._constantsInfo is not None:
                # Insert the constant parameters into the list of arguments at the specified indexes
                for info,constant in zip(self._constantsInfo,constants):
                    args.insert(info['argidx'],constant)
            # Evaluate and return the model response
            return Amodel(*args)
        #------------------------------------------------

        # Update the model function
        self.nonlinmodel = model_with_constants_and_added_nonlin
        self.signature.append(key)
    #---------------------------------------------------------------------------------------
    
    #---------------------------------------------------------------------------------------
    def rename_parameter(self, old, new):
        """
        Rename a parameter in the model. 

        This method renames a parameter in the model. If the
        old parameter name does not exist, a ``KeyError`` is raised.

        Parameters
        ----------
        old : string
            Old parameter name.
        new : string
            New parameter name. 

        """
        # Check if the old parameter exists in the model
        if not hasattr(self,old):
            raise KeyError(f'The model does not have a "{old}" parameter.')
        # Rename the parameter
        setattr(self, new, getattr(self, old))
        # Delete the old parameter
        delattr(self, old)
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def addlinear(self, name, vec=1, normalization=None, lb=-np.inf, ub=np.inf, par0=None, unit=None, description=None):
        """
        Add a new linear parameter (:ref:`Parameter` object) to the model. 

        Parameters
        ----------
        name : string
            Identifier of the parameter. This name will be used to refer
            to the parameter in the model.

        vec : int scalar, optional
            Number of elements in the parameter. If ``vec>1`` then the parameter will represent a 
            vector of linear parameters of length ``vec``. By default, a scalar parameter is defined.

        normalization : callable, optional
            Normalization function of the parameter. If specified, upon fitting the parameter will be normalized and 
            an the normalization factor will be reported separately. Does not add an additional normalization factor parameter. 
            Must be function taking the parameter and returning the normalized parameter.
            
        lb : float or array_like, optional
            Lower bound of the parameter. For vectorized parameters, must be a 
            vector with ``vec`` elements. If not specified, it is set to ``-np.inf``.

        ub : float or array_like, optional
            Lower bound of the parameter. For vectorized parameters, must be a
            vector with ``vec`` elements.  If not specified, it is set to ``+np.inf``.


        description : string, optional 
            Description of the parameter. 

        unit : string, optional
            Physical unit of the parameter.
        """
        # Check if the parameter already exists
        self._error_if_already_exists(name)

        # Check if the parameter is a vector
        if vec>1:
            idx = np.arange(self.Nparam,self.Nparam+vec) 
            self.Nparam += vec
            self.Nlin += vec
            newparam = Parameter(name=name, linear=np.full(vec,True), parent=self, idx=idx, par0=np.full(vec,par0), lb=np.full(vec,lb), ub=np.full(vec,ub), value=np.full(vec,None),frozen=np.full(vec,False), unit=unit, description=description)
        else:
            idx = self.Nparam
            self.Nparam += 1
            self.Nlin += 1
            newparam = Parameter(name=name, linear=True, parent=self, idx=idx, par0=par0, lb=lb, ub=ub, unit=unit, description=description)

        # Check if the normalization function is valid
        newparam.normalization = normalization 
        if newparam.normalization is not None:
            if not callable(normalization):
                raise TypeError('The normalization function must be a callable function.')
                newparam.description = str(newparam.description) + ' (normalized)' 

        # Add linear parameter to the model
        setattr(self,name,newparam)
        self.signature.append(name)
    #---------------------------------------------------------------------------------------


    def __call__(self,*args,**kargs):
    #---------------------------------------------------------------------------------------
        """
        Evaluate the model for a given set of parameters. 

        Takes the constant variables and (non-linear and linear) parameter 
        variables as positional or keyword arguments and evaluates the model. 

        Parameters
        ----------
        args : 
            Positional arguments representing the constant and parameter variables in the model. 

        kargs :
            Keyword arguments representing the constant and parameter variables in the model.

        Returns
        -------
        y : ndarray 
            The output of the model for the given set of parameters.

        """
        # Check that the correct number of arguments have been specified
        paramlist = self._parameter_list(order='vector')
        constlist = [info['argkey'] for info in self._constantsInfo]
        Nrequired = len(paramlist)
        Nrequired += len(constlist)
        if (len(args)+len(kargs))!=Nrequired:
            raise SyntaxError(f'The model requires {Nrequired} arguments, but {len(args)+len(kargs)} have been specified.')

        # Extract positional arguments  
        args_constants= [np.atleast_1d(args[info['argidx']]) for info in self._constantsInfo if info['argidx']<len(args)]     
        args_list = [np.atleast_1d(arg) for idx,arg in enumerate(args) if idx not in [info['argidx'] for info in self._constantsInfo]] 

        # Check keyword arguments 
        paramlist_ = paramlist[len(args_list):]
        karg_keys = kargs.keys()
        for key in karg_keys:
            if key in paramlist and key not in paramlist_: 
                raise SyntaxError(f'The parameter "{key}" has been specified twice.')
            if key not in paramlist and not key in constlist: 
                errstr = f'The argument \"{key}\" is not part of the model signature.'
                proposal = difflib.get_close_matches(key, paramlist)
                print(paramlist,proposal)
                if len(proposal)>0:
                    errstr += f' \n\t\tDid you mean: {proposal}?'
                raise AttributeError(errstr)    
        for param in paramlist_:
            if param not in karg_keys:
                raise KeyError(f'The parameter "{param}" has not been specified.')

        # Extract keywords arguments      
        kargs_constants = [np.atleast_1d(kargs[info['argkey']]) for info in self._constantsInfo if info['argkey'] in kargs]  
        kargs_list = [np.atleast_1d(kargs[param]) for param in paramlist_]

        constants = args_constants + kargs_constants
        param_list = args_list + kargs_list

        # Concatente all parameter into a single vector
        θ = np.concatenate(param_list)

        # Check that all parameters have been passed
        if len(θ)!=self.Nparam:
            raise SyntaxError(f'The model requires {self.Nparam} parameters, but {len(args_list)} were specified.')   

        # Determine which parameters are linear and which nonlinear
        θlin, θnonlin = self._split_linear(θ)

        # Calculate the design matrix (possibly dependent on non-linear parameters)
        A = self.nonlinmodel(*constants, *θnonlin)
        A = np.atleast_2d(A)

        # Ensure adequate linear-parameter vector shape
        θlin = np.atleast_1d(np.squeeze(θlin))
        if len(θlin)==0: 
            θlin = np.array([1])
        if A.shape[1]!=len(θlin): 
            A = A.T

        # Full model evaluation 
        y = A@θlin

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

        Returns
        -------
        metadata : dict
            Dictionary containing all the model's metadata in ordered vectors. The model keys 
            correspond to the model attributes, e.g. ``metadata['lb']`` corresponds to the model's
            parameters lower boundaries. 
        """
        # Return a dictionary of metadata vectors, with keys corresponding to the metadata attributes
        return {
            'names': self._parameter_list(order='vector'),
            'ub' : self._vecsort(self._getvector('ub')),
            'lb' : self._vecsort(self._getvector('lb')),
            'par0' : self._vecsort(self._getvector('par0')),
            'frozen' : self._vecsort(self._getvector('frozen')),
            'linear' : self._vecsort(self._getvector('linear')),
            'values' : self._vecsort(self._getvector('value')),
            'units' : self._vecsort(self._getvector('unit')),
            }
    #---------------------------------------------------------------------------------------

        
    #---------------------------------------------------------------------------------------
    def _parameter_table(self):
        """
        (Private method)

        Construct a parameter table for the model that describes each of the model's
        parameters.

        Returns
        -------
        table : str
            Formatted string representation of the model's parameter table.
        """

        # Get the description and signature of the model
        string = inspect.cleandoc(f"""
    Description: {self.description}
    Signature: ({', '.join(self.signature)})
    Constants: [{', '.join([entry['argkey'] for entry in self._constantsInfo])}]
    Parameter Table: 
    """)
        string += '\n'
        table = []
        # Add the table header to the list
        table.append(['Name','Lower','Start','Upper','Type','Frozen','Unit','Description'])  
        # Alignment for the table columns
        alignment = ['<','^','^','^','^','^','^','<']
        # Get a list of the model's parameters, in vector order
        for paramname in self._parameter_list(order='vector'): 
            # Build strings with the parameter's metadata 
            param_str = paramname
            lb_str = f'{np.atleast_1d(getattr(self,paramname).lb)[0]:5.3g}'
            ub_str = f'{np.atleast_1d(getattr(self,paramname).ub)[0]:5.3g}'
            par0_str = f'{np.atleast_1d(getattr(self,paramname).par0)[0]:5.3g}' if np.atleast_1d(getattr(self,paramname).par0)[0] is not None else "-"
            linear_str = "linear" if np.all(getattr(self,paramname).linear) else "nonlin"
            frozen_str = f'{np.atleast_1d(getattr(self,paramname).value)[0]:5.3g}' if np.all(getattr(self,paramname).frozen) else "No"
            unit_str = str(getattr(self,paramname).unit)
            desc_str = str(getattr(self,paramname).description)

            # Add the parameter's information to the table
            table.append([param_str,lb_str,par0_str,ub_str,linear_str,frozen_str,unit_str,desc_str])

        # Convert the table list to a formatted string
        string += formatted_table(table,alignment)
        return string
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def __str__(self):
        """
        Return a string representation of the model's parameters and their values.
        """
        return self._parameter_table()      
    #---------------------------------------------------------------------------------------
#===================================================================================

#==============================================================================
class Penalty():
    r"""Represents a penalty term of the objective function.  

    Attributes
    ----------
    weight : ``Parameter`` instance  
        Penalty weight parameter (a ``Parameter`` object instance without the
        ``linear`` and ``par0`` attributes). 
    description : string 
        Description of the penalty. 
    selection : string 
        Name of the selection functional for the penalty weight optimization. 

    """
    #--------------------------------------------------------------------------
    def __init__(self,penaltyfcn,selection,description=None,signature=None):
        r"""
        Construct a new penalty object. 

        Parameters
        ----------
        penaltyfcn : callable 
            Function that takes a set of parameters and returns a vector that
            will internally be squared and appended to the least-squares 
            residual vector. The names of the arguments defined in the function
            signature must match the names of the parameter in the model used along
            the penalty. 
        
        selection : string 
            Selection functional for the outer optimization of the penalty weight. 

            - ``'aic'`` - Akaike information criterion
            - ``'bic'`` - Bayesian information criterion
            - ``'aicc'`` - COrrected Akaike information criterion
            - ``'icc'`` - Informational complexity criterion 

        description : string, optional 
            Description of the penalty.
        
        signature : list of strings
            Signature of the ``penaltyfcn`` function to manually specify the names
            of the input arguments. For internal use (mostly).

        """ 
        #-------------------------------------------------------------------------------
        def selectionfunctional(fitfcn,y,sigma,log10weight):
            """
            (Private function)

            Calculate the selection functional used to find the optimal penalty weight.

            Parameters
            ----------
            fitfcn : callable
                Function that performs the fit to the data. 
            y : array_like
                Data to fit the model to.
            sigma : array_like
                Standard deviation of the data.
            log10weight : float
                Penalty weight in logarithmic scale.

            Returns
            -------
            selection_functional : float
                Value of the selection functional calculated using the data, the standard deviation, 
                the fit function and the penalty weight.
            """
            # Penalty weight: linear-scale -> log-scale
            weight = 10**log10weight
            self._weight_value = weight
            # Run the fit
            fitresult = fitfcn(weight)

            if selection=='icc':
                yfit = fitresult.model
                if isinstance(yfit,list):
                    # Get the fitted model
                    yfit = np.concatenate(yfit)

                # Get non-linear parameters covariance submatrix
                fitpars = fitresult.nonlin + np.finfo(float).eps
                covmat = fitresult.nonlinUncert.covmat + np.finfo(float).eps
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
        r"""
        Optimize the penalty weight. 

        Parameters
        ----------
        fitfcn : callable
            Fit function taking a penalty weight value. Must return
            a :ref:`FitResult` object.
        y : array_like 
            Dataset being fitted 
        sigma : scalar 
            Estimated noise level (standard deviation). 
        
        Returns 
        -------
        fit : :ref:`FitResult`
            Fit at the optimized penalty weight. 
        weightopt : scalar 
            Optimized penalty weight.
        """
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
            self._weight_value = optweight

        # Update optimized value to object
        self.optweight = optweight

        # Get the fit result
        fitresult = fitfcn(optweight)

        return fitresult,optweight
    #--------------------------------------------------------------------------
#==============================================================================


# ---------------------------------------------------------------------
def _importparameter(parameter):
    """
    (Private function)

    Import relevant information from a parameter object.

    This method extracts the lower and upper bounds (``lb``, ``ub``), initial value
    (``par0``), description, unit, frozen state, and current value of the input
    parameter object, and returns them as a dictionary.

    Parameters
    ----------
    parameter : ``Parameter`` object
        The parameter object whose information is to be extracted.

    Returns
    -------
    parameter_info : dict
        A dictionary containing the extracted information of the input
        parameter object.

    """
    # Extract the relevant information from the parameter object
    parameter_info = {
        'lb': parameter.lb,
        'ub': parameter.ub,
        'par0': parameter.par0,
        'description': parameter.description,
        'unit': parameter.unit,
        'frozen': parameter.frozen,
        'value': parameter.value,
    }

    return parameter_info
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
def _aresame(obj1,obj2):
    """
    (Private function)

    Compare two objects and return True if they have the same attribute values,
    excluding the parent attribute and the index attribute.

    Parameters
    ----------

    obj1 : object
        First object to compare.

    obj2 : object
        Second object to compare.

    Returns
    -------

    bool
        True if the objects have the same attribute values, False otherwise.
    """
    # Get the attribute dictionaries for both objects, excluding the parent and index attributes
    a = obj1.__dict__
    a = {key:val for key, val in a.items() if key not in ['_parent','idx']}
    b = obj2.__dict__
    b = {key:val for key, val in b.items() if key not in ['_parent','idx']}

    # Compare the dictionaries and return True if they are equal, False otherwise
    try:
        np.testing.assert_equal(a,b)
        return True
    except Exception:
        return False
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
def _linked_model_with_constants(nonlinfcn,mapping,constantsInfo,linear_reduce_idx,*inputargs):
    """
    (Private function)

    Wrapper function that redistributes the input parameter vector and returns the mapped model response.

    Parameters
    ----------
    nonlinfcn : callable
        Function that takes a set of non-linear parameters and returns the model response.
    
    mapping : ndarray
        Mapping vector specifying how the input parameter vector should be redistributed.
    
    constantsInfo : list
        List of dictionaries with information on the constants in the function.
    
    linear_reduce_idx : list
        List of arrays with the indices of the model response to be summed together.

    Returns
    -------
    A : ndarray
        Model response with the input parameter vector redistributed according to the mapping vector.
    """
    # Get the number of constant parameters
    Nconstants = len(constantsInfo)
    # separate the constants and model in the input arguments
    constants = inputargs[:Nconstants]
    θ = inputargs[Nconstants:]
    # Rearrange the model parameters according to the mapping vector
    θ = np.atleast_1d(θ)[mapping]
    # Create a list of input arguments to the model's non-linear function
    args = list(θ)
    # Insert the constant parameters into their specified positions in the argument list
    if constantsInfo is not None:
        for info,constant in zip(constantsInfo,constants):
            args.insert(info['argidx'],constant)                
    # Evaluate the model's non-linear function with the updated argument list
    A = nonlinfcn(*args)
    # Return the output of the non-linear function if it is a tuple (for internal use, undocumented)
    if isinstance(A,tuple): return A
    # Make a matrix if model function returns a vector
    if len(A.shape)<2: A = np.expand_dims(A,1)
    # Sum the output matrix along the columns specified in the linear_reduce_idx array
    Amapped = np.vstack([np.sum(np.atleast_2d(A[:,idx]),axis=1) for idx in linear_reduce_idx]).T
    return Amapped
# ---------------------------------------------------------------------

# ==============================================================================
def link(model,**links):
    """
    Link parameters together in a model to create equality relationships 
    between parameters. 

    Link different sets of parameters to single parameters in a model, such that
    the values of the linked parameters are all equal to the value of those 
    single parameters. The linked parameters are removed from the model and 
    only the single parameters remain.

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
        """
        (Private function)

        Link parameters together in a model.
        
        Link a series of parameters to a single parameter in a model, such that
        the values of the linked parameters are all equal to the value of the 
        single parameter. The linked parameters are removed from the model and 
        only the single parameter remains.
        
        Parameters
        ----------
        model : object
            Model object instance.
        parameters : list
            List of parameters to be linked together. 
        newname : string
            Name of the linked parameter.

        Returns
        -------
        model : object
            Updated model object instance with the linked parameters removed.
        """
        
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

        # Redefine the non-linear part of the model function
        model.nonlinmodel = partial(_linked_model_with_constants,nonlinfcn,mapping,model._constantsInfo,linear_reduce_idx)

        # Return the updated model with the linked parameters
        return model
    # ---------------------------------------------------------------------

    if not isinstance(model,Model):
        raise TypeError('The first argument must be a Model object.')
    newmodel = deepcopy(model)
    # Perform the linking, one by one
    for link_newname in links: 
        to_link = [getattr(newmodel,parname) for parname in links[link_newname]]
        newmodel = _linkparameter(newmodel,to_link,link_newname)
    # Update the new model signature
    for key in links.keys():
        newmodel.signature = [key if arg==links[key][0] else arg for arg in newmodel.signature]
        for arg in links[key]:
            if arg in newmodel.signature: 
                newmodel.signature.remove(arg)
    return newmodel
#==============================================================================================

# ---------------------------------------------------------------------
def _unique_ordered(vec):
    """
    (Private function)

    Returns a list of the unique elements in vec, ordered as they appear in vec.

    Parameters
    ----------
    vec : list
        List of elements to be processed.

    Returns
    -------
    uniques : list
        List of the unique elements in vec, ordered as they appear in vec.

    """
    # Create a list to store the unique elements
    uniques = []
    # Loop through the elements of vec
    for v in vec: 
        # If the element is not already in uniques, add it
        if v not in uniques:
            uniques.append(v)
    # Return the list of unique elements
    return uniques
# ---------------------------------------------------------------------


#---------------------------------------------------------------------
def _combined_nonlinmodel(mode,nonlinfcns,Nlins,Nconst,arelinear,const_subsets,subsets_nonlin,*inputargs):
    """
    (Private function)

    Combine the output of multiple non-linear models.

    This function is used to combine the output of multiple non-linear models
    into a single design matrix or model response. The combination mode
    can be either 'merge' (merge the output of each model into a single block-
    diagonal matrix) or 'lincombine' (concatenate the output of each model
    horizontally).

    Parameters
    ----------
    mode : str
        Either 'merge' or 'lincombine', depending on the desired
        combination mode.
    nonlinfcns : list of callables
        List of non-linear functions to combine.
    Nlins : list of ints
        List of the number of linear parameters in each non-linear function.
    Nconst : int
        Number of constants for each non-linear function.
    arelinear : list of bools
        List of bools indicating whether each non-linear function has
        only linear parameters.
    const_subsets : list of lists
        List of constant subsets for each non-linear function.
    subsets_nonlin : list of lists
        List of non-linear parameter subsets for each non-linear function.
    *inputargs : list of arrays
        List of input arguments to the non-linear functions, consisting
        of constants followed by parameters.

    Returns
    -------
    Anonlin_full : ndarray
        Array with the combined output of the non-linear functions.
    """
    # Unpack the input arguments into constants and parameters
    constants = inputargs[:Nconst]
    param = inputargs[Nconst:]

    # Make sure that constants and parameters are arrays
    param = np.atleast_1d(param)
    constants = np.atleast_2d(constants)

    # Loop over the submodels in the model
    Amatrices = []
    for n,nonlinfcn in enumerate(nonlinfcns):
        # Evaluate the submodel
        Amatrix = np.atleast_2d(nonlinfcn(*constants[const_subsets[n],:],*param[subsets_nonlin[n]]))

        # Transpose the output matrix if necessary
        if np.shape(Amatrix)[1]!=Nlins[n]:
            Amatrix = Amatrix.T
        Amatrices.append(Amatrix)

    # Merge or concatenate the output of the submodels
    if mode=='merge':
        Anonlin_full = block_diag(Amatrices).toarray()

        # If all submodels are non-linear, sum the output matrices
        if not any(arelinear):
            Anonlin_full = np.sum(Anonlin_full,1)

    elif mode=='lincombine':
        Anonlin_full = np.hstack(Amatrices)

    # Return the combined output
    return Anonlin_full
#---------------------------------------------------------------------

#---------------------------------------------------------------------
def _split_output(nonlinfcns,Nlins,Nconst,const_subsets,subsets_nonlin,y,*inputargs):
    """
    (Private function)

    Split the output of a model into the outputs of its submodels.

    This method is used to evaluate the outputs of each of the submodels that
    make up a composite model. It is called internally by the ``Model`` class.

    Parameters
    ----------
    nonlinfcns : list of callables
        List of functions representing the submodels that make up the composite model.

    Nlins : list of int
        List of the number of linear parameters in each of the submodels.

    Nconst : int
        Total number of constant arguments for all submodels.

    const_subsets : list of array_like
        List of indices specifying which constant arguments each submodel uses.

    subsets_nonlin : list of array_like
        List of indices specifying which non-linear parameters each submodel uses.

    y : array_like
        Output of the full model.

    *inputargs
        List of constant arguments followed by the non-linear parameters of the full model.
        
    Returns
    -------
    output : list of array_like
        List of the outputs of each submodel.
    """
    # Split the input arguments into constants and parameters
    constants = inputargs[:Nconst]
    param = inputargs[Nconst:]
    # Ensure that variables are arrays
    param = np.atleast_1d(param)
    constants = np.atleast_2d(constants)
    
    # Loop over the submodels in the model
    ysizes = []
    for n,nonlinfcn in enumerate(nonlinfcns):
        # Evaluate the submodel
        Amatrix = np.atleast_2d(nonlinfcn(*constants[const_subsets[n],:],*param[subsets_nonlin[n]]))
        # Transpose the matrix if necessary
        if np.shape(Amatrix)[1]!=Nlins[n]: Amatrix = Amatrix.T
        # Store the size of the output from this submodel
        ysizes.append(Amatrix.shape[0])
    # Create subsets of y for each submodel output
    nprev = 0
    ysubsets = []
    for x in ysizes:
        ysubsets.append(np.arange(nprev,nprev+x))
        nprev = nprev+x
    # Return the subsets of y
    return [y[ysubsets[n]] for n in range(len(ysizes))]
#---------------------------------------------------------------------


#---------------------------------------------------------------------
def _combinemodels(mode,*inputmodels,addweights=False): 
    """
    (Private function)

    Combine multiple ``Model`` objects into a single model. This helper function
    is used by the ``merge`` and ``lincombine`` functions.

    Parameters
    ----------
    mode : str
        The mode of combination, which can be ``'parallel'`` or ``'sequential'``. In
        parallel mode, each model is assumed to be independent of each other
        and the combined model returns the sum of the output of the
        individual models. In sequential mode, the output of each model is
        passed on to the next model as input.

    inputmodels : ``Model`` objects
        Any number of models to be combined.

    addweights : bool, optional
        If ``True``, add weighting factors to the combined model. These weighting
        factors can be used to adjust the contribution of each individual
        model to the combined model. The default is ``False``.

    Returns
    -------
    combmodel : Model object
        The combined model.
    """
    # Initialize empty containers
    subsets_nonlin,arguments,arelinear,lin_normalizations = [],[],[],[]
    nprev = 0

    if len(inputmodels)==1:
        return inputmodels[0]

    # Make deep-copies of the models to avoid modifying them
    models = [deepcopy(model) for model in inputmodels]

    if addweights:

        for n,(model,nonlinfcn) in enumerate(zip(models,[model.nonlinmodel for model in models])):
            constants = [constant['argkey'] for constant in model._constantsInfo]
            signature = []
            for param in model.signature: 
                if param in constants:
                    signature.append(param)                        
                else:
                    if not np.any(getattr(model,param).linear):
                        signature.append(param)                        
    
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
            model.addlinear('scale',par0=1,lb=0, description='Scaling factor')

        # Determine the subset of parameters for the current model
        subset = np.arange(nprev,nprev+model.Nnonlin,1)
        nprev += model.Nnonlin
        # From that subset, determine the non-linear subset
        subset_nonlin = subset[np.arange(model.Nnonlin)]
        subsets_nonlin.append(subset_nonlin)

        # Determine which parameters are linear
        arelinear = np.concatenate([arelinear,model._vecsort(model._getvector('linear'))])
        lin_normalizations += [getattr(model,param).normalization for param in model._parameter_list() if hasattr(getattr(model,param),'normalization') ]

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

    # Create the model object
    combinedModel = Model(partial(_combined_nonlinmodel,mode,nonlinfcns,Nlins,Nconst,arelinear,const_subsets,subsets_nonlin),constants=constants,signature=signature)

    # Add parent models 
    combinedModel.parents = models

    if mode=='merge':
        # Add post-evalution function for splitting of the call outputs
        setattr(combinedModel,'_posteval_fcn',partial(_split_output,nonlinfcns,Nlins,Nconst,const_subsets,subsets_nonlin)) 

    # Add the linear parameters from the subset models   
    lin_param_set = []
    for param, lin_normalization in zip(_unique_ordered(lin_params),lin_normalizations):
        lin_param_set.append({'name':param, 'vec':np.sum(lin_params==param), 'normalization':lin_normalization})

    for lparam in lin_param_set:
        combinedModel.addlinear(lparam['name'], vec=lparam['vec'], normalization=lparam['normalization'])

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
#---------------------------------------------------------------------

#==============================================================================================
def merge(*inputmodels,addweights=False):
    """
    Create a multi-response model from multiple individual models. 

    Parameters
    ----------
    inputmodels : :ref:`Model` objects
        Model objects to be combined. If one of the models has no linear parameters, a linear 
        scaling factor parameters will be added. The names of the ``N``-th input model parameter will be 
        changed by a suffix ``_N`` in the new model. Example:: 

            newmodel = merge(model1,model2)
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
    return _combinemodels('merge',*inputmodels,addweights=addweights)
#==============================================================================================

#==============================================================================================
def lincombine(*inputmodels,addweights=False):
    """
    Create model whose response is a linear combination of multiple individual model responses. 

    Parameters
    ----------
    inputmodels : :ref:`Model` objects
        Model objects whose linear responses are to be linearly combined. If one of the models 
        has no linear parameters, a linear scaling factor parameters will be added. The names 
        of the ``N``-th input model parameter will be changed by a suffix ``_N`` in the new model. Example:: 

            newmodel = lincombine(model1,model2)
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
    return _combinemodels('lincombine',*inputmodels,addweights=addweights)
#==============================================================================================

# ---------------------------------------------------------------------
def _dependency_model_with_constants(function,nonlinfcn,constantsInfo,arguments_idx,dependent_idx,*inputargs):
    """
    (Private function)

    A helper function that maps the input arguments of a model and
    calls the nonlinear function with the mapped arguments.

    This function is used to evaluate models that have dependencies
    between their parameters, where one parameter is defined in terms
    of another.

    Parameters
    ----------
    function : callable
        The function defining the dependency between parameters.

    nonlinfcn : callable
        The nonlinear model function to be called with the mapped
        input arguments.

    constantsInfo : list
        A list of dictionaries containing information about the
        constants in the model.

    arguments_idx : ndarray
        An array of indices specifying which arguments of the
        nonlinear model function are arguments of the dependency
        function.

    dependent_idx : int
        The index of the dependent parameter in the input argument
        vector.

    *inputargs
        The input arguments of the model.

    Returns
    -------
    A : array
        The model (design) matrix output of the nonlinear function with the mapped input
        arguments.
    """
    # Redistribute the input parameter vector according to the mapping vector
    Nconstants = len(constantsInfo)
    constants = inputargs[:Nconstants]
    θ = np.atleast_1d(inputargs[Nconstants:]).astype(float)

    # Insert the output of the dependency function at the correct index
    θ = np.insert(θ,dependent_idx,function(*θ[arguments_idx]))  

    # Build the list of arguments to be passed to the nonlinear function
    args = list(θ)
    if constantsInfo is not None:
        for info,constant in zip(constantsInfo,constants):
            args.insert(info['argidx'],constant)                

    # Call the nonlinear function with the mapped arguments
    A = nonlinfcn(*args)
    return A
# ---------------------------------------------------------------------

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

        if np.all(getattr(model,dependent_name).linear):
            raise TypeError(f"Linear parameters cannot be used.")

        for arg in arguments_names:
            if arg not in model_parameters:
                raise KeyError(f"The function argument '{arg}' is not a parameter of the input model.")
            if np.all(getattr(model,arg).linear):
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
        nonlinparams = np.array([param for param in model._parameter_list('vector') if not np.all(getattr(model,param).linear)])
        arguments_idx = np.concatenate([np.where(nonlinparams==name)[0] for name in arguments_names])
        dependent_idx = dependent_idx[0]
        
        # Redefine the non-linear part of the model function
        model.nonlinmodel = partial(_dependency_model_with_constants,function,nonlinfcn,model._constantsInfo,arguments_idx,dependent_idx)

        # Return the updated model with the linked parameters
        return model
    # ---------------------------------------------------------------------

    if not isinstance(model,Model):
        raise TypeError('The first argument must be a deerlab.Model object.')
    newmodel = deepcopy(model)

    # Get the dependent's names and their function arguments
    dependents = [dependent for dependent in functions]
    arguments = [inspect.getfullargspec(functions[dependent]).args for dependent in dependents]

    # Update the new model signature
    for arg in [item for sublist in arguments for item in sublist]:
        if arg in newmodel.signature: 
            newmodel.signature.remove(arg)

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
