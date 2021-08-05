import numpy as np
import deerlab as dl
import matplotlib.pyplot as plt 
from deerlab.solvers import rlls,nlls,snlls
import inspect 


#===================================================================================
class Parameter(): 
    #-------------------------------------------------------------------------------
    def __init__(self, parent=None, idx=None, name=None, par0=None, lb=-np.inf, ub=np.inf, units=None, linear=False): 
        # Attributes
        self._parent = parent # Parent 
        self.idx = idx
        self.name = name      # Name
        self.units = units    # Units
        self.par0 = par0      # Start values
        self.lb = lb          # Lower bounds
        self.ub = ub          # Upper bounds
        self.frozen = False   # Frozen
        self.linear = linear  # Linearity
    #-------------------------------------------------------------------------------

    #-------------------------------------------------------------------------------
    def set(self,**attributes):
        for key in attributes:
            setattr(self,key,attributes[key])
        return 
    #-------------------------------------------------------------------------------

#===================================================================================

#===================================================================================
class Model():
    #-------------------------------------------------------------------------------
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
    #-------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def addlinear(self, key, vec=None, lb=-np.inf, ub=np.inf, par0=None, name=None, units=None):
        """
        Add a new linear parameter.
        """
        if vec is not None: 
            idx = np.arange(self.Nparam,self.Nparam+vec)
            self.Nparam += vec        
            self.Nlin += vec
            newparam = Parameter(linear=np.full(vec,True), parent=self, idx=idx, par0=np.full(vec,par0), lb=np.full(vec,lb), ub=np.full(vec,ub), units=units, name=name)
        else:
            idx = self.Nparam
            self.Nparam += 1
            self.Nlin += 1
            newparam = Parameter(linear=True, parent=self, idx=idx, par0=par0, lb=lb, ub=ub, units=units, name=name)
        setattr(self,key,newparam)

    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def _core_model(self,Amodel,θnonlin,θlin): 
        """ 
        Core mathematical model 
        -----------------------

        Calculates the model ``y`` response based on the mathematical expression ``y = A(θnonlin)@θlin``.
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
        "Sort the parameters from alphabetical ordering to vector indexing"
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
        # If there are no linear parameters defined
        if len(variable_lin)==0: 
            variable_lin = 1
        return variable_lin, variable_nonlin
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

        return self._core_model(self.nonlinmodel,θnonlin,θlin)
    #---------------------------------------------------------------------------------------

    #---------------------------------------------------------------------------------------
    def _getvector(self,attribute):
        "Get the list of parameters attributes defined in the model sorted alphabetically"
        return [getattr(getattr(self,param),attribute) for param in dir(self) if isinstance(getattr(self,param),Parameter)]
    #---------------------------------------------------------------------------------------

    def fit(self,data,par0=None,**kwargs):
    #---------------------------------------------------------------------------------------

        # Get boundaries and conditions for the linear and nonlinear parameters
        ubl,ub = self._split_linear(self._vecsort(self._getvector('ub')))
        lbl,lb = self._split_linear(self._vecsort(self._getvector('lb')))
        if par0 is None:
            _,par0 = self._split_linear(self._vecsort(self._getvector('par0')))

        if np.any(par0==None):
            raise RuntimeError(f"It appears some start values (par0) have not been specified. Either specify them in the model definition or using the keyword.")

        if self.Nlin>0 and self.Nnonlin==0:      
            Amatrix = self.nonlinmodel()
            fitResult = rlls(data,Amatrix,lbl,ubl,**kwargs)

        elif self.Nlin>0 and self.Nnonlin>0:        
            Amodel = lambda param: np.atleast_2d(self.nonlinmodel(*param))
            fitResult = snlls(data,Amodel,par0,lb,ub,lbl,ubl,**kwargs)        

        elif self.Nlin==0 and self.Nnonlin>0:        
            model_fcn = lambda param: self.nonlinmodel(*param)
            fitResult = nlls(data,model_fcn,par0,lb,ub,**kwargs)
        else:
            raise AssertionError(f'The model has no parameters to fit.')
        return fitResult
    #---------------------------------------------------------------------------------------
    


#===================================================================================

