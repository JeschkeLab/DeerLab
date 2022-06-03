from copy import Error
import warnings
import numpy as np
import scipy as scp
import scipy.optimize as opt
from types import FunctionType
from functools import wraps
import pickle

def parse_multidatasets(V_, K, weights, noiselvl, precondition=False, masks=None, subsets=None):
#===============================================================================
    
    # Make copies to avoid modifying the originals
    V = V_.copy()

    Vlist = []
    # If multiple signals are specified as a list...
    if type(V) is list and all([type(Vs) is np.ndarray for Vs in V]):
        nSignals = len(V)

    elif type(V) is np.ndarray:
        nSignals = 1
        V = [V]
    else:
        raise TypeError('The input signal(s) must be numpy array or a list of numpy arrays.')
    
    # Pre-scale the signals, important when using global fits with arbitrary scales
    prescales = np.ones(nSignals)
    for i in range(nSignals):
        if precondition:
            prescales[i] = max(V[i])
            Vlist.append(V[i]/prescales[i])
        else:
            Vlist.append(V[i])

    # Get the indices to extract the subsets again
    Ns = [len(V) for V in Vlist]
    if subsets is None:
        subset = [None]*nSignals
        for i in  range(len(Vlist)):
            if i==0:
                prev = 0
            else:
                prev = subset[i-1][-1]+1
            subset[i] = np.arange(prev,prev+Ns[i])
    else: 
        subset = subsets.copy()        

    # Noise level estimation/parsing
    sigmas = np.zeros(len(subset))
    if noiselvl is None:
        for i in range(len(subset)):
            sigmas[i] = der_snr(Vlist[i])
    else: 
        noiselvl = np.atleast_1d(noiselvl)
        sigmas = noiselvl.copy()

    if masks is None: 
        masks = [np.full_like(V,True).astype(bool) for V in Vlist]
    elif not isinstance(masks,list):
        masks = [masks]
    if len(masks)!= len(Vlist): 
        raise SyntaxError('The number of masks does not match the number of signals.')
    mask = np.concatenate(masks, axis=0) 

    # Concatenate the datasets along the list axis
    V = np.concatenate(Vlist, axis=0)

    # -----------------------------------------------------------------
    def prepareKernel(K,nSignals):
        # If multiple kernels are specified as a list...
        if type(K) is tuple:
            K = [Ks for Ks in K]
        if type(K) is list and all([type(Ks) is np.ndarray for Ks in K]):
            nKernels = len(K)
            K = np.concatenate(K, axis=0) # ...concatenate them along the list 
        elif type(K) is np.ndarray:
            nKernels = 1
        else:
            raise TypeError('The input kernel(s) must be numpy array or a list of numpy arrays.')
        # Check that the same number of signals and kernel have been passed
        if nSignals!=nKernels:
            raise KeyError('The same number of kernels and signals must be specified as lists.')
        return K
    # -----------------------------------------------------------------
    if type(K) is FunctionType:
        Kmulti = lambda p: prepareKernel(K(p),nSignals)
    elif K is None: 
        Kmulti = None
    else:
        Kmulti = prepareKernel(K,nSignals)

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
        if len(weights)!=len(Vlist) and len(weights)!=np.sum([len(V) for V in Vlist]):
            raise KeyError('If multiple signals are passed, the same number of weights are required.')
        if len(weights)!=np.sum([len(V) for V in Vlist]):
            weights_ = []
            for i in range(len(weights)):
                weights_ = np.concatenate((weights_,weights[i]*np.ones(len(Vlist[i]))))
            weights = weights_
    else:
        raise TypeError('The input weights(s) must be numpy array or a list of numpy arrays.')


    if precondition:
        return V, Kmulti, weights, mask, subset, sigmas, prescales
    else:
        return V, Kmulti, weights, mask, subset, sigmas
#===============================================================================

#===============================================================================
def formatted_table(table,align=None):
    """Generate auto-formatted table in string form from a list of rows."""
    
    # Determine the maximal number of characters in each column
    N = []
    for column in range(len(table[0])):    
        N.append(max([len(str(row[column])) for row in table]))

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
def der_snr(y):
    """
    DER-SNR noise estimation 

    Estimate the noise level using the DER_SNR method [1]_.
    
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


# =================================================================
def metadata(**kwargs):
    """
    Decorator: Set model metadata as function attributes 
    """
    attributes = list(kwargs.keys())
    metadata = list(kwargs.values())

    def _setmetadata(func):
        for attribute,data in zip(attributes,metadata):
            setattr(func,attribute,data)
        return func
    return _setmetadata
# =================================================================

def gsvd(A,B):
#===============================================================================
    m,p = A.shape
    n = B.shape[0]

    # Economy-sized.
    useQA = m > p
    useQB = n > p
    if useQA:
        QA,A = scp.linalg.qr(A)
        A = A[0:p,0:p]
        QA = QA[:,0:p]
        m = p

    if useQB:
        QB,B = scp.linalg.qr(B)
        B = B[0:p,0:p]
        QB = QB[:,0:p]
        n = p

    Q,_ = np.linalg.qr(np.vstack((A,B)), mode='reduced')
    Q1 = Q[0:m,0:p]
    Q2 = Q[m:m+n,0:p]
    C,S = csd(Q1,Q2)

    # Vector of generalized singular values.
    q = min(m+n,p)
    # Supress divide by 0 warning        
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        U = np.vstack((np.zeros((q-m,1),'double'), np.diag(C,max(0,q-m)).reshape(len(np.diag(C,max(0,q-m))),1)))/np.vstack((np.diag(S,0).reshape(len(np.diag(S,0)),1), np.zeros((q-n,1),'double') ))


    return U
#===============================================================================


def csd(Q1,Q2):
#===============================================================================
    """
    Cosine-Sine Decomposition
    -------------------------
    
    Given Q1 and Q2 such that Q1'* Q1 + Q2'* Q2 = I, the
    C-S Decomposition is a joint factorization of the form
          Q1 = U1*C*V' and Q2=U2*S*V'
    where U1,U2,V are orthogonal matrices and C and S are diagonal
    matrices (not necessarily square) satisfying
             C'* C + S'* S = I
    The diagonal entries of C and S are nonnegative and the
    diagonal elements of C are in nondecreasing order.
    The matrix Q1 cannot have more columns than rows.

    Based on the Octave code by Artiste (submitted by S.J.Leon): 
    http://www.ar-tiste.com/m-fun/m-fun-index.html

    """
    m,n = Q1.shape
    p,_ = Q2.shape
    if m < p:
        s,c = csd(Q2,Q1)
        j = np.flip(np.arange(n)) 
        c = c[:,j] 
        s = s[:,j] 
        m = np.minimum(m,p) 
        i =  np.flip(np.arange(m)) 
        c[np.arange(m),:] = c[i,:] 
        n = np.minimum(n,p) 
        i =  np.flip(np.arange(n)) 
        s[np.arange(n),:] = s[i,:] 
        return c,s

    _,sdiag,v = np.linalg.svd(Q1)
    c = np.zeros((m, n))
    np.fill_diagonal(c, sdiag)
    v = v.T.conj()
    z = np.eye(n,n)
    z = scp.linalg.hankel(z[:,n-1])
    c[0:n,:] = z@c[0:n,:]@z
    v = v@z
    Q2 = Q2@v
    k=0
    for j in range(1,n):
        if c[j,j] <= 1/np.sqrt(2): k=j
    b = Q2[:,0:k]
    u2,r = np.linalg.qr(b,mode='complete')
    s = u2.T@Q2
    t = np.minimum(p,n)
    tt = np.minimum(m,p)
    if k<t:
        r2 = s[np.ix_(range(k,p),range(k,t))]
        _,sdiag,vt = np.linalg.svd(r2)
        ss= np.zeros(r2.shape)
        np.fill_diagonal(ss, sdiag)
        vt = vt.T.conj()
        s[k:p,k:t] = ss
        c[:,k:t] = c[:,k:t]@vt
        w = c[k:tt,k:t]
        z,r = np.linalg.qr(w,mode='complete')
        c[k:tt,k:t] = r
    for j in range(n):
        if c[j,j]<0:
            c[j,j]  = -c[j,j]
    for j in range(t):
        if s[j,j]<0:
            s[j,j]  = -s[j,j]

    return c,s
#===============================================================================

#===============================================================================
def diagf(X):
    """
    Diagonal force

    X = diagf(X) zeros all the elements off the main diagonal of X.
    """
    X = np.triu(np.tril(X))
    return X
#===============================================================================

#===============================================================================
def diagp(Y,X,k):
    """
    Diagonal positive matrix.
    Y,X = diagp(Y,X,k) scales the columns of Y and the rows of X by
    unimodular factors to make the k-th diagonal of X real and positive.
    """
    D = np.diag(X,k)
    j = np.where((D.real < 0) | (D.imag != 0))
    D = np.diag(np.conj(D[j])/abs(D[j]))
    Y[:,j] = Y[:,j]@D.T
    X[j,:] = D@X[j,:]
    X = X+0 # use "+0" to set possible -0 elements to 0
    return Y,X
#===============================================================================

#===============================================================================
def Jacobian(fcn, x0, lb, ub):
    """ 
    Finite difference Jacobian estimation 
     
    Estimates the Jacobian matrix of a vector-valued function `f(x)` at the 
    point `x_0` taking into consideration box-constraints of the vector-valued variable ``x`` defined by its lower
    and upper bounds.

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

    Parameters
    ----------
    A : ndarray 
        Matrix

    Returns
    -------
    Cpsd : ndarray 
        Nearest positive semi-definite matrix 

    Notes
    -----

    Modified from the algorithm in [1]_.

    References
    ----------

    .. [1] tjiagoM (2020, July 28th), 
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

def isempty(A):
#===============================================================================Q
    A = np.atleast_1d(A)
    boolean = np.size(A)==0
    return boolean
#===============================================================================

def isnumeric(obj):
#===============================================================================
    attrs = ['__add__', '__sub__', '__mul__', '__truediv__', '__pow__']
    return all(hasattr(obj, attr) for attr in attrs)
#===============================================================================

def multistarts(n,x0,lb,ub):
#===============================================================================
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
        Skip a unit test in Pytest if a particular exception occurs during the execution of the test.
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
        input_args = inspect.getfullargspec(function).args
        docstring = function.__doc__

        for arg in input_args:
            try:
                assert arg in docstring
            except: 
                raise AssertionError(f'The argument {arg} is not documented.')
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
        Name (and path) of the file to save the pickled object to. The object is saved as a ``.pkl`` file. 
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
