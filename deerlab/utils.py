from copy import Error
import warnings
import numpy as np
import scipy as scp
import scipy.optimize as opt
from types import FunctionType
from functools import wraps
import pickle

def parse_multidatasets(y_, A, weights, noiselvl, precondition=False, masks=None, subsets=None):
#===============================================================================
    """
    Parse one or multiple dataset(s) for local/global least-squares inference

    This function takes an arbitrary number of datasets, along with the model matrices and noise
    level for each dataset, and parses them into a common format for use in a global fit. It also
    performs some preprocessing steps such as scaling the signals and estimating the noise level
    if it is not provided.

    Parameters
    ----------
    y_ : array_like or list of array_like
        The input dataset(s). Must be a numpy array or a list of numpy arrays.
    
    A : array_like, list of array_like, or callable
        The model matrix(ces) for the dataset(s). Must be a numpy array or a list of numpy arrays, or a callable that returns such an object.
    
    weights : array_like or None, optional
        The weights for the global fit. If not provided, default weights will be used based on the noise level of each dataset.
    
    noiselvl : float or array_like, optional
        The noise level(s) of the dataset(s). If not provided, the noise level will be estimated using the `der_snr` function.
    
    precondition : bool, optional
        Whether to scale the signals to have maximum value of 1. Default is False.
    
    masks : array_like or list of array_like, optional
        Masks to be applied to the dataset(s). If not provided, all data will be used.
    
    subsets : list of array_like, optional
        The indices of the datasets in the concatenated global dataset. If not provided, the indices will be automatically determined.

    Returns
    -------
    y : array_like
        The concatenated dataset(s), with each dataset preprocessed as specified.
    
    Kmulti : array_like or callable
        The concatenated model matrix(ces), or a callable that returns the concatenated model matrix(ces) when given a set of parameters.
    
    subset : list of array_like
        The indices of the datasets in the concatenated global dataset.
    
    mask : array_like
        The concatenated mask(s), applied to the dataset(s).
    
    sigmas : array_like
        The noise level(s) of the dataset(s).
    
    weights : array_like
        The weights for the global fit.
    """
    # Make copies to avoid modifying the originals
    y = y_.copy()

    ylist = []
    # If multiple datasets are specified as a list...
    if type(y) is list and all([type(ys) is np.ndarray for ys in y]):
        nDatasets = len(y)

    elif type(y) is np.ndarray:
        nDatasets = 1
        y = [y]
    else:
        raise TypeError('The input dataset(s) must be numpy array or a list of numpy arrays.')
    
    # Pre-scale the signals, important when using global fits with arbitrary scales
    prescales = np.ones(nDatasets)
    for i in range(nDatasets):
        if precondition:
            prescales[i] = max(y[i])
            ylist.append(y[i]/prescales[i])
        else:
            ylist.append(y[i])

    # Get the indices to extract the subsets again
    Ns = [len(y) for y in ylist]
    if subsets is None:
        subset = [None]*nDatasets
        for i in  range(len(ylist)):
            if i==0:
                prev = 0
            else:
                prev = subset[i-1][-1]+1
            subset[i] = np.arange(prev,prev+Ns[i])
    else: 
        subset = subsets.copy()        

    # Parse the masks
    if masks is None: 
        masks = [np.full_like(y,True).astype(bool) for y in ylist]
    elif not isinstance(masks,list):
        masks = [masks]
    if len(masks)!= len(ylist): 
        raise SyntaxError('The number of masks does not match the number of signals.')
    mask = np.concatenate(masks, axis=0) 

    # Noise level estimation/parsing
    sigmas = np.zeros(len(subset))
    if noiselvl is None:
        for i in range(len(subset)):
            sigmas[i] = der_snr(ylist[i][masks[i]])
    else: 
        noiselvl = np.atleast_1d(noiselvl)
        sigmas = noiselvl.copy()


    # Concatenate the datasets along the list axis
    y = np.concatenate(ylist, axis=0)

    # -----------------------------------------------------------------
    def prepare_model_matrix(A,nDatasets):
        """
        Prepare model matrix for global fit

        This function takes an arbitrary number of model matrices, and parses
        them into a common format for use in a global fit. It concatenates 
        them along the list axis if they are specified as a list.

        Parameters
        ----------
        A : array_like, list of array_like, or callable
            The input model matrix(ces). Must be a numpy array or a list of numpy arrays, or a callable that returns such an object.
        nDatasets : int
            The number of datasets in the global fit.

        Returns
        -------
        A : array_like or callable
            The concatenated model matrix(ces), or a callable that returns the concatenated model matrix(ces) when given a set of parameters.
        """
        # If multiple kernels are specified as a list...
        if type(A) is tuple:
            A = [As for As in A]
        if type(A) is list and all([type(As) is np.ndarray for As in A]):
            nKernels = len(A)
            A = np.concatenate(A, axis=0) # ...concatenate them along the list 
        elif type(A) is np.ndarray:
            nKernels = 1
        else:
            raise TypeError('The input model matrix must be numpy array or a list of numpy arrays.')
        # Check that the same number of signals and kernel have been passed
        if nDatasets!=nKernels:
            raise KeyError('The same number of model matrices and datasets must be specified as lists.')
        return A
    # -----------------------------------------------------------------
    if type(A) is FunctionType:
        Kmulti = lambda p: prepare_model_matrix(A(p),nDatasets)
    elif A is None: 
        Kmulti = None
    else:
        Kmulti = prepare_model_matrix(A,nDatasets)

    # If global weights are not specified, set default based on noise levels
    if weights is None:
        if len(subset)==1:
            weights=1
        else:
            weights = np.zeros(len(subset))
            for i in range(len(subset)):
                weights[i] = 1/sigmas[i]
    
    # If multiple weights are specified as a list...
    if type(weights) is list or type(weights) is np.ndarray or not hasattr(weights, "__len__"):
        weights = np.atleast_1d(weights)
        if len(weights)==1:
                weights = np.repeat(weights,len(subset))
        if len(weights)!=len(ylist) and len(weights)!=np.sum([len(y) for y in ylist]):
            raise KeyError('If multiple signals are passed, the same number of weights are required.')
        if len(weights)!=np.sum([len(y) for y in ylist]):
            weights_ = []
            for i in range(len(weights)):
                weights_ = np.concatenate((weights_,weights[i]*np.ones(len(ylist[i]))))
            weights = weights_
    else:
        raise TypeError('The input weights(s) must be numpy array or a list of numpy arrays.')

    if precondition:
        return y, Kmulti, weights, mask, subset, sigmas, prescales
    else:
        return y, Kmulti, weights, mask, subset, sigmas
#===============================================================================


#===============================================================================
def formatted_table(table,align=None):
    """
    Generate formatted table in string form

    This function takes a table as a list of rows, and returns a string representation 
    of the table with automatic formatting. The table is formatted using a row-separator
    line, a header row, and the data rows.

    Parameters
    ----------
    table : list of list of objects
        The input table, as a list of rows, where each row is a list of objects representing
        the columns.
    align : list of str, optional
        The alignment of the columns in the table. Must be a list of strings, where each string is 
        either '<' for left alignment, '>' for right alignment, or '^' for center alignment. If 
        not provided, all columns will be left-aligned.

    Returns
    -------
    table2print : str
        The formatted table, as a string.
    """
    from pyparsing import Literal,Combine,Suppress,Optional,delimitedList,Word,nums,alphas,oneOf
    ESC = Literal('\x1b')
    escapeSeq = Combine(ESC + '[' + Optional(delimitedList(Word(nums),';')) + oneOf(list(alphas)))
    rm_ansi = lambda s : Suppress(escapeSeq).transformString(s)

    # Determine the maximal number of characters in each column
    N = []
    for column in range(len(table[0])):    
        row_lengths = [len(rm_ansi(str(row[column]))) for row in table]
        N.append(max(row_lengths))

    # If not specified, use left-alignment
    if align is None:
        align = ['<']*len(N)

    # Construct string for row-separator
    table_lines = ''
    for n,a in zip(N,align):
        table_lines += "="*(n+2) 
        table_lines += ' '

    # Construct row formatter string
    formatter = ''
    for n,a in zip(N,align):
        formatter = formatter + ' {:' + f'{a}' + f'{n}' + '}  ' 

    # Construct the table
    table2print = '' 
    table2print += table_lines + '\n'  # Add a line
    table2print += formatter.format(*table.pop(0)) + '\n' # Add a header 
    table2print += table_lines  + '\n' # Add a line
    for row in table:
        table2print += formatter.format(*row) + '\n' # Add data rows
    table2print += table_lines # Add a line
    return table2print
#===============================================================================


#===============================================================================
def goodness_of_fit(x,xfit,Ndof,noiselvl):
    """
    Goodness of Fit statistics

    Computes multiple statistical indicators of goodness of fit.

    Parameters
    ----------
    x : array_like
        Original data.
    xfit : array_like
        Fit.
    Ndof : int
        Number of degrees of freedom.
    noiselvl : float
        Standard deviation of the noise in x.

    Returns
    -------
    stats : dict
        Statistical indicators:
            stats['chi2red'] - Reduced chi-squared
            stats['rmsd'] - Root mean-squared deviation
            stats['R2'] - R-squared test
            stats['aic'] - Akaike information criterion
            stats['aicc'] - Corrected Akaike information criterion
            stats['bic'] - Bayesian information criterion
    """
    sigma = noiselvl
    Ndof = np.maximum(Ndof,1)
    residuals = x - xfit
    
    # Special case: no noise, and perfect fit
    if np.isclose(np.sum(residuals),0):
        return {'chi2red':1,'R2':1,'rmsd':0,'aic':-np.inf,'aicc':-np.inf,'bic':-np.inf,'autocorr':0}

    # Get number of xariables
    N = len(x)
    # Extrapolate number of parameters
    Q = Ndof - N

    # Reduced Chi-squared test
    chi2red = 1/Ndof*np.linalg.norm(residuals)**2/sigma**2

    # Autocorrelation based on Durbin–Watson statistic
    autocorr_DW = abs( 2 - np.sum((residuals[1:-1] - residuals[0:-2])**2)/np.sum(residuals**2) )

    # R-squared test
    R2 = 1 - np.sum((residuals)**2)/np.sum((xfit-np.mean(xfit))**2)

    # Root-mean square dexiation
    rmsd = np.sqrt(np.sum((residuals)**2)/N)

    # Log-likelihood
    loglike = N*np.log(np.sum((residuals)**2))

    # Akaike information criterion
    aic =  loglike + 2*Q

    # Corrected Akaike information criterion
    aicc = loglike + 2*Q + 2*Q*(Q+1)/(N-Q-1)

    # Bayesian information criterion
    bic =  loglike + Q*np.log(N)

    return {'chi2red':chi2red,'R2':R2,'rmsd':rmsd,'aic':aic,'aicc':aicc,'bic':bic,'autocorr':autocorr_DW}
#===============================================================================




#===============================================================================
def der_snr(y):
    """
    DER-SNR noise estimation 

    Estimate the noise level using the DER_SNR method [1].
    
    Parameters
    ----------
    y : array_like 
        Noisy dataset.

    Returns 
    -------
    sigma : float scalar 
        Noise standard deviation. 

    References
    ---------- 
    .. [1] F. Stoehr, R. White, M. Smith, I. Kamp, R. Thompson, D. Durand, W. Freudling,
       D. Fraquelli, J. Haase, R. Hook, T. Kimball, M. Kummel, K. Levay, M. Lombardi, A. Micol, T. Rogers 
       DERSNR: A Simple & General Spectroscopic Signal-to-Noise Measurement Algorithm
       Astronomical Data Analysis Software and Systems XVII, ASP Conference Series, Vol. 30, 2008, p5.4
    """

    n = len(y)
    sigma  = 1.482602/np.sqrt(6)*np.median(abs(2.0*y[2:n-2] - y[0:n-4] - y[4:n]))
    
    return sigma
#===============================================================================

#===============================================================================
def hccm(J,residual,mode='HC1'):
    """
    Heteroscedasticity Consistent Covariance Matrix (HCCM)

    Computes the heteroscedasticity consistent covariance matrix (HCCM) of
    a given LSQ problem given by the Jacobian matrix and the residual vector
    of a least-squares problem. The HCCM are valid for both heteroscedasticit and
    homoscedasticit residual vectors. 

    Parameters
    ----------
    J : NxM-element ndarray
        Jacobian matrix of the residual vector
    residual : N-element ndarray
        Vector of residuals
    mode : string, optional
        HCCM estimation method:

        * ``'HC0'`` - White, 1980 [1]_
        * ``'HC1'`` - MacKinnon and White, 1985 [2]_
        * ``'HC2'`` - MacKinnon and White, 1985 [2]_
        * ``'HC3'`` - Davidson and MacKinnon, 1993 [3]_
        * ``'HC4'`` - Cribari-Neto, 2004 [4]_
        * ``'HC5'`` - Cribari-Neto, 2007 [5]_

        If not specified, it defaults to ``'HC1'``

    Returns
    -------
    C : MxM-element ndarray
       Heteroscedasticity consistent covariance matrix 

    References
    ---------- 
    .. [1] White, H. (1980). A heteroskedasticity-consistent covariance matrix
       estimator and a direct test for heteroskedasticity. Econometrica, 48(4), 817-838
       DOI: 10.2307/1912934

    .. [2] MacKinnon and White, (1985). Some heteroskedasticity-consistent covariance
       matrix estimators with improved finite sample properties. Journal of Econometrics, 29 (1985), 
       pp. 305-325. DOI: 10.1016/0304-4076(85)90158-7

    .. [3] Davidson and MacKinnon, (1993). Estimation and Inference in Econometrics
       Oxford University Press, New York. 

    .. [4] Cribari-Neto, F. (2004). Asymptotic inference under heteroskedasticity of
       unknown form. Computational Statistics & Data Analysis, 45(1), 215-233
       DOI: 10.1016/s0167-9473(02)00366-3

    .. [5] Cribari-Neto, F., Souza, T. C., & Vasconcellos, K. L. P. (2007). Inference
       under heteroskedasticity and leveraged data. Communications in Statistics –
       Theory and Methods, 36(10), 1877-1888. DOI: 10.1080/03610920601126589
    """

    # Hat matrix
    H = J@np.linalg.pinv(J.T@J)@J.T
    # Get leverage
    h = np.diag(H)
    # Number of parameters (k) & Number of variables (n)
    n,k = np.shape(J)

    # IF the number of parameters and variables are equal default to the HC0 mode to avoid zero-division
    if n==k: mode='HC0'

    # Select estimation method using established nomenclature
    if mode.upper() == 'HC0': # White,(1980),[1]
        # Estimate the data covariance matrix
        V = np.diag(residual**2)
        
    elif mode.upper() == 'HC1': # MacKinnon and White,(1985),[2]
        # Estimate the data covariance matrix
        V = n/(n-k)*np.diag(residual**2)
        
    elif mode.upper() == 'HC2': # MacKinnon and White,(1985),[2]
        # Estimate the data covariance matrix
        V = np.diag(residual**2/(1-h))
        
    elif mode.upper() == 'HC3': # Davidson and MacKinnon,(1993),[3]
        # Estimate the data covariance matrix
        V = np.diag(residual/(1-h))**2
        
    elif mode.upper() == 'HC4': # Cribari-Neto,(2004),[4]
        # Compute discount factor
        delta = np.minimum(4,n*h/k)
        # Estimate the data covariance matrix
        V = np.diag(residual**2./((1 - h)**delta))
        
    elif mode.upper() == 'HC5': # Cribari-Neto,(2007),[5]
        # Compute inflation factor
        k = 0.7
        alpha = np.minimum(np.maximum(4,k*max(h)/np.mean(h)),h/np.mean(h))
        # Estimate the data covariance matrix
        V = np.diag(residual**2./(np.sqrt((1 - h)**alpha)))
            
    else:
        raise KeyError('HCCM estimation mode not found.')

    # Heteroscedasticity Consistent Covariance Matrix (HCCM) estimator
    C = np.linalg.pinv(J.T@J)@J.T@V@J@np.linalg.pinv(J.T@J)

    # Ensure that the covariance matrix is positive semi-definite
    C = nearest_psd(C)

    return C
#===============================================================================

#===============================================================================
def Jacobian(fcn, x0, lb, ub):
    """ 
    Finite difference Jacobian estimation 
     
    Estimates the Jacobian matrix of a vector-valued function `f(x)` at the 
    point `x_0` taking into consideration box-constraints of the vector-valued 
    variable ``x`` defined by its lower and upper bounds.

    Parameters
    ----------
    fcn : callable 
        Function `f(x)` to be differentiated.     
    x0 : ndarray
        Point `x_0` at which to differentiate. 
    lb : ndarray 
        Lower bounds of `x`.  
    ub : ndarray 
        Upper bounds of `x`. 

    Notes
    -----

    This is a wrapper around the ``scipy.optimize._numdiff.approx_derivative``
    function of the Scipy package.

    """
    J = opt._numdiff.approx_derivative(fcn,x0,method='2-point',bounds=(lb,ub))
    J = np.atleast_2d(J) 
    return J
#===============================================================================

#===============================================================================
def movmean(x, N):
    """
    Moving mean filter

    Returns an array of local `N`-point mean values, where each mean is calculated
    over a sliding window of length `k` across neighboring elements of `x`.

    Parameters
    ----------
    x : ndarray
        Array to be filtered
    N : integer scalar
        Window size
    
    Returns
    -------
    xfilt : ndarray
        Filtered array
    
    """
    xfilt = np.convolve(x, np.ones(N)/N, mode='same')
    return xfilt
#===============================================================================

#===============================================================================
def ovl(A,B):
    """
    Overlap index

    Returns the overlap index between two vectors A and B with respect to zero.

    Parameters
    ----------
    A : N-element ndarray
        First vector
    B : N-element ndarray
        Second vector
    
    Returns
    -------
    metric : array
        Overlap metric

    """
    A /= np.sum(A)
    B /= np.sum(B)
    metric = np.sum(np.minimum(A,B))
    return metric
#===============================================================================

#===============================================================================
def nearest_psd(A):
    """
    Find the nearest positive semi-definite matrix

    This function takes a square matrix `A` and returns the nearest positive semi-definite matrix.

    Parameters
    ----------
    A : ndarray
        The input matrix.

    Returns
    -------
    Cpsd : ndarray
        The nearest positive semi-definite matrix.

    Notes
    -----

    Modified from the algorithm in [1].

    References
    ----------

    [1] tjiagoM (2020, July 28th), 
        How can I calculate the nearest positive semi-definite matrix? 
        StackOverflow, https://stackoverflow.com/a/63131250/16396391
    """
    # If matrix is empty, return it empty (scipy.linalg.eigh cannot deal with empty matrix)
    if A.size==0: 
        return A
    # Symmetrize the matrix
    Asym = (A + A.T)/2
    # Construct positive semi-definite matrix via eigenvalue decomposition
    eigval, eigvec = scp.linalg.eigh(Asym)
    eigval[eigval < 0] = 0
    Cpsd = np.real(eigvec.dot(np.diag(eigval)).dot(eigvec.T))
    # Avoid round-off errors
    Cpsd[abs(Cpsd)<=np.finfo(float).eps] = 0
    return Cpsd
#===============================================================================

#===============================================================================
def isnumeric(obj):
    """
    Check whether an object is numeric

    This function checks whether an object supports the basic operations of a numeric type,
    such as addition, subtraction, multiplication, and division.

    Parameters
    ----------
    obj : object
        The object to be checked.

    Returns
    -------
    result : bool
        `True` if `obj` is numeric, `False` otherwise.
    """
    attrs = ['__add__', '__sub__', '__mul__', '__truediv__', '__pow__']
    return all(hasattr(obj, attr) for attr in attrs)
#===============================================================================

#===============================================================================
def multistarts(n,x0,lb,ub):
    """
    Generate multiple starting points

    This function generates multiple starting points for optimization, based on a given initial point,
    lower bounds, and upper bounds. If `n` is positive, `n-1` new points are generated within the bounds,
    excluding the bounds themselves.

    Parameters
    ----------
    n : int
        The number of starting points to generate.
    x0 : array_like
        The initial point.
    lb : array_like
        The lower bounds for the starting points.
    ub : array_like
        The upper bounds for the starting points.

    Returns
    -------
    x0 : ndarray
        The starting points.
    """
    # Ensure a 1D-vector
    x0 = np.atleast_1d(np.squeeze(x0)).astype(float)
    lb = np.atleast_1d(np.squeeze(lb)).astype(float)
    ub = np.atleast_1d(np.squeeze(ub)).astype(float)

    if n<0:
        raise ValueError('The number of requested starting points must be n>0.') 

    if len(x0) != len(lb) or len(x0) != len(ub):
        raise ValueError('The lower/upper bound size(s) are not compatible with the initial guess vector x0.') 

    _x0 = x0.copy()
    # Generate n-1 new starting points within the bounds (excluding the bounds)
    if n>1:
        x0 = np.linspace(lb,ub,n+2)[1:-1,:]
    else:
        x0 = [x0]

    # If there is some NaN or inf value (just put the original start value)
    for x, in zip(x0):
        x[np.isnan(x)] = _x0[np.isnan(x)]
        x[np.isinf(x)] = _x0[np.isinf(x)]
    
    return x0
#===============================================================================

# This is only required for the developers version which installs pytest automatically
try:
    import pytest     
    #===============================================================================
    def skip_on(errclass, reason="Default reason"):
        """ 
        Decorator for skipping failed tests during pytest's execution

        Skip a unit test in Pytest if a particular exception occurs during the execution of the test.
        When an error of the type `errclass` occurs, the test will be skipped and the given `reason` will be reported.

        Examples
        --------
        @skip_on(TypeError, reason="Data is not a string")
        def test_parse_string():
            assert parse_string(1234) == "1234"

        Source: modified from https://stackoverflow.com/a/63522579/16396391
        """ 
        # Func below is the real decorator and will receive the test function as param
        def decorator_func(f):
            @wraps(f)
            def wrapper(*args, **kwargs):
                try:
                    # Try to run the test
                    return f(*args, **kwargs)
                except BaseException as error:
                    if errclass in str(type(error)):
                        # If exception of given type happens
                        # just swallow it and raise pytest.Skip with given reason
                        pytest.skip(reason)
                    else: raise Error(error)
            return wrapper

        return decorator_func
    #===============================================================================


    import inspect
    #===============================================================================
    def assert_docstring(function):
        """
        Check if the docstring of a given function contains documentation for all of its input arguments.
        
        Parameters
        ----------
        function : callable
            The function to check.
        
        Raises
        ------
        AssertionError
            If the docstring does not contain documentation for all of the function's input arguments.
        """
        input_args = inspect.getfullargspec(function).args + inspect.getfullargspec(function).kwonlyargs
        docstring = function.__doc__

        for arg in input_args:
            try:
                assert arg in docstring
            except: 
                raise AssertionError(f'The argument "{arg}" is not documented.')
    #===============================================================================
except: pass

import dill as pickle

# --------------------------------------------------------------------------------------
def store_pickle(obj, filename):
    """
    Save/export an object to a ``.pkl`` file serialized as bytes.
    
    Parameters
    ----------
    obj : object 
        Python object to be saved/exported. 

    filename : string 
        Name (and path) of the file to save the pickled object to. The object is saved as a ``.pkl`` file. If `filename` does not end in ".pkl", the extension is added automatically.

    Examples
    --------
    store_pickle(my_object, 'my_file.pkl')
    store_pickle(my_object, './data/my_file') # equivalent to './data/my_file.pkl'
    """ 
    if '.pkl' not in filename:
        filename = filename + '.pkl'
    with open(filename, 'wb') as outp:  # Overwrites any existing file.
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)
# --------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------
def read_pickle(filename):
    """
    Load a pickled object file ``.pkl`` and deserialize the bytes into a Python object.
    
    Parameters
    ----------
    filename : string 
        Path to the ``.pkl`` file to load.

    Returns
    -------
    object
        The deserialized Python object.

    .. warning:: It is possible to construct malicious pickle data which will execute arbitrary code during unpickling. Never unpickle data that could have come from an untrusted source, or that could have been tampered with. See `here <https://docs.python.org/3/library/pickle.html>`_ for more information.

    """ 
    if '.pkl' not in filename:
        filename = filename + '.pkl'
    with open(filename, "rb") as f:
        while True:
            try:
                return pickle.load(f)
            except EOFError:
                break
# --------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------
def sophegrid(octants,maxphi,size):
    """
    Construct spherical grid over spherical angles based on input parameters. 
    The grid implemented in this function is often called the SOPHE grid [1]_. 
    Adapted from Easyspin [2]_ ``sphgrid_`` source code. 

    Parameters
    ----------
    octants : integer 
        Number of "octants" of the sphere; for each increment in theta, octants additional points are added along phi; special cases: octants=0 and octants=-1.
    maxphi : float 
        Largest value of angle phi (radians).
    size : integer  
        Number of orientations between theta=0 and theta=pi/2. 

    Returns
    -------
    phi : ndarray
        Array of the phi angles of the grid points in radians.
    theta : ndarray
        Array of the theta angles of the grid points in radians.
    weights : ndarray
        Array of the weights corresponding to each grid point.

    References
    ----------

    .. [1] D. Wang, G. R. Hanson
       J.Magn.Reson. A, 117, 1-8 (1995)
       https://doi.org/10.1006/jmra.1995.9978

    .. [2] Stefan Stoll, Arthur Schweiger
       J. Magn. Reson. 178(1), 42-55 (2006)
       https://doi.org/10.1016/j.jmr.2005.08.013
    """ 
    dtheta = (np.pi/2)/(size-1); # angular increment along theta
    if octants > 0: # if not Dinfh or O3 symmetry
        
        # Initializations
        if octants==8:
            nOct = 4
        else:
            nOct = octants
        nOrientations = int(size + nOct*size*(size-1)/2)
        phi = np.zeros(nOrientations)
        theta = np.zeros(nOrientations)
        weights = np.zeros(nOrientations)
        
        sindth2 = np.sin(dtheta/2)
        w1 = 1.0
        
        # North pole (z orientation)
        phi[0] = 0
        theta[0] = 0
        weights[0] = maxphi*(1 - np.cos(dtheta/2))
        
        # All but equatorial slice
        Start = 1
        for iSlice in np.arange(2,size):
            nPhi = nOct*(iSlice-1) + 1
            dPhi = maxphi/(nPhi - 1)
            idx = Start + np.arange(0,nPhi)
            weights[idx] = 2*np.sin((iSlice - 1)*dtheta)*sindth2*dPhi*np.concatenate([[w1], np.ones(nPhi-2), [0.5]])
            phi[idx] = np.linspace(0,maxphi,nPhi)
            theta[idx] = (iSlice - 1)*dtheta
            Start = Start + nPhi
        
        # Equatorial slice
        nPhi = nOct*(size - 1) + 1
        dPhi = maxphi/(nPhi - 1)
        idx = Start + np.arange(0,nPhi)
        phi[idx] = np.linspace(0,maxphi,nPhi)
        theta[idx] = np.pi/2
        weights[idx] = sindth2*dPhi*np.concatenate([[w1], np.ones(nPhi-2), [0.5]])
        
        # Border removal
        rmv = np.cumsum(nOct*np.arange(1,size)+1)
        phi = np.delete(phi,rmv)
        theta = np.delete(theta,rmv)
        weights = np.delete(weights,rmv)

        # For C1, add lower hemisphere
        if octants==8:
            idx = len(theta)-nPhi + np.arange(1,1,-1) -1
            phi = np.concatenate([phi, phi[idx]])
            theta = np.concatenate([theta, np.pi-theta[idx]])
            weights[idx] = weights[idx]/2; # half of real value
            weights = np.concatenate([weights, weights[idx]])
        
        weights = 2*(2*np.pi/maxphi)*weights; # sum = 4*pi

    elif octants==0: # Dinfh symmetry (quarter of meridian in xz plane)

        phi = np.zeros(1,size)
        theta = np.linspace(0,np.pi/2,size)
        weights = -2*(2*np.pi)*np.diff(np.cos(np.concatenate([[0], np.arange(dtheta/2,np.pi/2,dtheta), [np.pi/2]]))); # sum = 4*pi

    elif octants==-1: # O3 symmetry (z orientation only)
        
        phi = 0
        theta = 0
        weights = 4*np.pi

    else:    
        raise ValueError('Unsupported value #d for octants.',octants)

    # Remove orientations with zero weight
    phi = phi[weights!=0]
    theta = theta[weights!=0]
    weights = weights[weights!=0]
    
    # Normalize to unity 
    weights = weights/np.sum(weights)

    return phi,theta,weights
# --------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------
def choleskycovmat(Q,cholfactors):
    """
    Compute the covariance matrix using Cholesky decomposition.

    Parameters
    ----------
    Q : int
        Dimension of the covariance matrix.
    cholfactors : array_like
        List of Cholesky decomposition factors. The first `Q` values correspond to the 
        diagonal values of the Cholesky decomposition matrix. The remaining values are used to
        fill the lower triangular matrix row by row. 

    Returns
    -------
    Σ : ndarray 
        Covariance matrix.

    Examples
    --------
    >>> choleskycovmat(3, [3.0, 2.0, 1.0, 0.5, 0.5, 0.5])
    array([[3.  , 1.5 , 0.75],
           [1.5 , 2.25, 1.25],
           [0.75, 1.25, 1.5 ]])
    """
    
    #Check that the correct number of Cholesky factors has been specified
    N = Q*(Q+1)/2
    if len(cholfactors)!=N: 
        raise ValueError(f'The Cholesky decomposition factors vector must have {N} elements')

    # Initialize Cholesky decomposition matrix
    L = np.zeros((Q,Q))

    # Fill the Cholesky factors as a lower triangular matrix
    for q in range(Q):
        tril_idx = np.where(np.eye(Q,k=q) == 1)
        L[tril_idx] = cholfactors[:(Q-q)]
        cholfactors = cholfactors[(Q-q):]
        
    # Cholesky decomposition of the covariance matrix
    Σ = L@L.T

    return Σ 
# ----------------------------------------------------------------------------------
