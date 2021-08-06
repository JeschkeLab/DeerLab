import numpy as np
import deerlab as dl
import matplotlib.pyplot as plt 
from deerlab.solvers import rlls,nlls,snlls
import inspect 


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
    def __init__(self, parent=None, idx=None, name=None, par0=None, frozen=False, lb=-np.inf, ub=np.inf,value=None, units=None, linear=False): 
        # Attributes
        self._parent = parent # Parent 
        self.idx = idx
        self.name = name      # Name
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
    def __init__(self,Amodel): 
        """
        Model object constructor
        """
        if not callable(Amodel):
            Amatrix = Amodel.copy()
            Amodel = lambda: Amatrix 
        self.nonlinmodel = Amodel
        parameters = inspect.getfullargspec(Amodel).args

        # Update the number of parameters in the model
        self.Nparam = len(parameters)
        self.Nnonlin = len(parameters)
        self.Nlin = 0

        for n,param in enumerate(parameters): 
            newparam = Parameter(parent=self, idx=n)
            setattr(self,param,newparam)
    #---------------------------------------------------------------------------------------

    #=======================================================================================
    #                                    Auxiliary functions
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

        if A.shape[1]!=len(θlin): 
            A = A.T

        # Full model calculation 
        return A@θlin
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def _parameter_list(self):
        "Get the list of parameters defined in the model sorted alphabetically"
        return [param for param in dir(self) if isinstance(getattr(self,param),Parameter)]
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def _vecsort(self,list):
        "Sort vectorized parameters attributes from alphabetical ordering to vector indexing"
        list = np.squeeze(np.atleast_1d(list))
        indices = np.concatenate([np.atleast_1d(getattr(self,param).idx) for param in dir(self) if isinstance(getattr(self,param),Parameter)])
        orderedlist = list.copy()
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
        return [getattr(getattr(self,param),attribute) for param in dir(self) if isinstance(getattr(self,param),Parameter)]
    #---------------------------------------------------------------------------------------

    #=======================================================================================
    #                                         Methods
    #=======================================================================================

    #---------------------------------------------------------------------------------------
    def addlinear(self, key, vec=None, lb=-np.inf, ub=np.inf, par0=None, name=None, units=None):
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
        if vec is not None: 
            idx = np.arange(self.Nparam,self.Nparam+vec)
            self.Nparam += vec        
            self.Nlin += vec
            newparam = Parameter(linear=np.full(vec,True), parent=self, idx=idx, par0=np.full(vec,par0), lb=np.full(vec,lb), ub=np.full(vec,ub), value=np.full(vec,None),frozen=np.full(vec,False), units=units, name=name)
        else:
            idx = self.Nparam
            self.Nparam += 1
            self.Nlin += 1
            newparam = Parameter(linear=True, parent=self, idx=idx, par0=par0, lb=lb, ub=ub, units=units, name=name)
        setattr(self,key,newparam)
    #---------------------------------------------------------------------------------------

    def __call__(self,*args,**kargs):
    #---------------------------------------------------------------------------------------
        if kargs and args: 
            raise SyntaxError('The model must be called either with positional or keyword arguments. Not both.')
        
        if args: 
            # Values are already ordered
            θ = np.concatenate([np.atleast_1d(arg) for arg in args]) 
        elif kargs:
            θ = np.concatenate([np.atleast_1d(kargs[param]) for param in self._parameter_list()])
            θ = self._vecsort(θ)

        # Determine which parameters are linear and which nonlinear
        θlin, θnonlin = self._split_linear(θ)
        # If there are no linear parameters defined
        if len(θlin)==0: 
            θlin = 1
        return self._core_model(self.nonlinmodel,θnonlin,θlin)
    #---------------------------------------------------------------------------------------


    def fit(self,y,par0=None,**kwargs):
    #---------------------------------------------------------------------------------------
        """
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
        # Get boundaries and conditions for the linear and nonlinear parameters
        ubl,ub = self._split_linear(self._vecsort(self._getvector('ub')))
        lbl,lb = self._split_linear(self._vecsort(self._getvector('lb')))
        frozenl,frozen = self._split_linear(self._vecsort(self._getvector('frozen')))
        valuesl,values = self._split_linear(self._vecsort(self._getvector('value')))

        # Check the initial conditions and whether they are defined
        if par0 is None:
            _,par0 = self._split_linear(self._vecsort(self._getvector('par0')))
        if np.any(par0==None):
            raise RuntimeError(f"It appears some start values (par0) have not been specified. Either specify them in the model definition or using the keyword.")

        linfrozen = np.full(self.Nlin,None)
        linfrozen[frozenl] = valuesl[frozenl]
        nonlinfrozen = np.full(self.Nnonlin,None)
        nonlinfrozen[frozen] = values[frozen]


        # Determine the class of least-squares problem to solve
        if self.Nlin>0 and self.Nnonlin==0:      
            # ---------------------------------------------------------
            # Linear LSQ  
            # ---------------------------------------------------------

            # Get the design matrix
            Amatrix = self.nonlinmodel()
            # Run penalized LSQ solver
            fitResult = rlls(y,Amatrix,lbl,ubl,frozen=linfrozen,**kwargs)

        elif self.Nlin>0 and self.Nnonlin>0:        
            # ---------------------------------------------------------
            # Separable non-linear LSQ  
            # ---------------------------------------------------------

            Amodel_fcn = lambda param: np.atleast_2d(self.nonlinmodel(*param))
            fitResult = snlls(y,Amodel_fcn,par0,lb,ub,lbl,ubl,lin_frozen=linfrozen,nonlin_frozen=nonlinfrozen,**kwargs)        

        elif self.Nlin==0 and self.Nnonlin>0:       
            # ---------------------------------------------------------
            # Non-linear LSQ  
            # ---------------------------------------------------------

            model_fcn = lambda param: self.nonlinmodel(*param)
            fitResult = nlls(y,model_fcn,par0,lb,ub,frozen=nonlinfrozen,**kwargs)
        else:
            raise AssertionError(f'The model has no parameters to fit.')
            
        return fitResult
    #---------------------------------------------------------------------------------------
    


#===================================================================================
