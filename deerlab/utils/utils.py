import warnings
import numpy as np
import cmath as math
import scipy as scp
import scipy.optimize as opt

from types import FunctionType 


def parse_multidatasets(V,K,weights,precondition=False):
#===============================================================================
    
    # Identify if the signals have already been processed by this function
    if type(V) is not list:
        if V.size == np.atleast_1d(weights).size:
            # If so, just return without doing anything
            if precondition:
                return V,K,weights,[np.arange(0,len(V))],[1]
            else:
                return V,K,weights,[np.arange(0,len(V))]

    # If multiple signals are specified as a list...
    if type(V) is list and all([type(Vs) is np.ndarray for Vs in V]):
        nSignals = len(V)
        prescales = np.zeros(nSignals)
        Vlist = []
        # Pre-scale the signals, important for fitregmodel when using global fits with arbitrary scales
        for i in range(nSignals):
            if precondition:
                prescales[i] = max(V[i])
                Vlist.append(V[i]/prescales[i])
            else:
                Vlist.append(V[i])
        V = np.concatenate(Vlist, axis=0) # ...concatenate them along the list 
    elif type(V) is np.ndarray:
        nSignals = 1
        prescales = [1]
        Vlist = [V]
    else:
        raise TypeError('The input signal(s) must be numpy array or a list of numpy arrays.')
    
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

    if type(K) is FunctionType:
        Kmulti = lambda p: prepareKernel(K(p),nSignals)
    else:
        Kmulti = prepareKernel(K,nSignals)

    # If multiple weights are specified as a list...
    if type(weights) is list or not hasattr(weights, "__len__"):
        weights = np.atleast_1d(weights)
        if len(weights)==1:
                weights = np.repeat(weights,nSignals)
        weights = weights/sum(weights)
        if len(weights)!=nSignals:
            raise KeyError('If multiple signals are passed, the same number of weights are required.')
        weights_ = []
        for i in range(len(weights)):
            weights_ = np.concatenate((weights_,weights[i]*np.ones(len(Vlist[i]))))
        weights = weights_
    else:
        raise TypeError('The input weights(s) must be numpy array or a list of numpy arrays.')

    # Get the indices to extract the subsets again
    Ns = [len(V) for V in Vlist]
    subset = [None]*nSignals
    for i in  range(nSignals):
        if i==0:
            prev = 0
        else:
            prev = subset[i-1][-1]+1
        subset[i] = np.arange(prev,prev+Ns[i])

    if precondition:
        return V,Kmulti,weights,subset,prescales
    else:
        return V,Kmulti,weights,subset
#===============================================================================


def hccm(J,*args):
    """
    Heteroscedasticity Consistent Covariance Matrix (HCCM)
    ======================================================

    Computes the heteroscedasticity consistent covariance matrix (HCCM) of
    a given LSQ problem given by the Jacobian matrix (J) and the covariance
    matrix of the data (V). If the residual (res) is specified, the
    covariance matrix is estimated using some of the methods specified in
    (mode). The HCCM are valid for both heteroscedasticit and
    homoscedasticit residual vectors. 

    Usage:
    ------
    C = hccm(J,V)
    C = hccm(J,res,mode)

    Arguments:
    ----------
    J (NxM-element array)
        Jacobian matrix of the residual vector
    res (N-element array)
        Vector of residuals
    mode (string)
        HCCM estimator, options are:
            'HC0' - White, H. (1980)
            'HC1' - MacKinnon and White, (1985)
            'HC2' - MacKinnon and White, (1985)
            'HC3' - Davidson and MacKinnon, (1993)
            'HC4' - Cribari-Neto, (2004)
            'HC5' - Cribari-Neto, (2007)

    Returns:
    --------
    C (MxM-element array) 
       Heteroscedasticity consistent covariance matrix 

    References:
    ------------ 
    [1] 
    White, H. (1980). A heteroskedasticity-consistent covariance matrix
    estimator and a direct test for heteroskedasticity. Econometrica, 48(4), 817-838
    DOI: 10.2307/1912934

    [2] 
    MacKinnon and White, (1985). Some heteroskedasticity-consistent covariance
    matrix estimators with improved finite sample properties. Journal of Econometrics, 29 (1985), 
    pp. 305-325. DOI: 10.1016/0304-4076(85)90158-7

    [3] 
    Davidson and MacKinnon, (1993). Estimation and Inference in Econometrics
    Oxford University Press, New York. 

    [4] 
    Cribari-Neto, F. (2004). Asymptotic inference under heteroskedasticity of
    unknown form. Computational Statistics & Data Analysis, 45(1), 215-233
    DOI: 10.1016/s0167-9473(02)00366-3

    [5] 
    Cribari-Neto, F., Souza, T. C., & Vasconcellos, K. L. P. (2007). Inference
    under heteroskedasticity and leveraged data. Communications in Statistics –
    Theory and Methods, 36(10), 1877-1888. DOI: 10.1080/03610920601126589
    """

    # Unpack inputs
    if len(args)==2:
        res,mode = args
        V = []
    elif len(args)==1:
        V = args[0]

    # Hat matrix
    H = J@np.linalg.pinv(J.T@J)@J.T
    # Get leverage
    h = np.diag(H)
    # Number of parameters (k) & Number of variables (n)
    n,k = np.shape(J)

    if isempty(V):
        # Select estimation method using established nomenclature
        if mode.upper() == 'HC0': # White,(1980),[1]
            # Estimate the data covariance matrix
            V = np.diag(res**2)
            
        elif mode.upper() == 'HC1': # MacKinnon and White,(1985),[2]
            # Estimate the data covariance matrix
            V = n/(n-k)*np.diag(res**2)
            
        elif mode.upper() == 'HC2': # MacKinnon and White,(1985),[2]
            # Estimate the data covariance matrix
            V = np.diag(res**2/(1-h))
            
        elif mode.upper() == 'HC3': # Davidson and MacKinnon,(1993),[3]
            # Estimate the data covariance matrix
            V = np.diag(res/(1-h))**2
            
        elif mode.upper() == 'HC4': # Cribari-Neto,(2004),[4]
            # Compute discount factor
            delta = np.minimum(4,n*h/k)
            # Estimate the data covariance matrix
            V = np.diag(res**2./((1 - h)**delta))
            
        elif mode.upper() == 'HC5': # Cribari-Neto,(2007),[5]
            # Compute inflation factor
            k = 0.7
            alpha = np.minimum(np.maximum(4,k*max(h)/np.mean(h)),h/np.mean(h))
            # Estimate the data covariance matrix
            V = np.diag(res**2./(np.sqrt((1 - h)**alpha)))
                
        else:
            raise KeyError('HCCM estimation mode not found.')


    # Heteroscedasticity Consistent Covariance Matrix (HCCM) estimator
    C = np.linalg.pinv(J.T@J)@J.T@V@J@np.linalg.pinv(J.T@J)

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
    DIAGP  Diagonal positive.
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
    Estimates the Jacobian matrix of a vector-valued function ``fcn`` at the 
    point ``x0`` taking into consideration box-constraints defined by the lower
    and upper bounds ``lb`` and ``ub``.

    This is a wrapper around the ``scipy.optimize._numdiff.approx_derivative`` function.

    """
    J = opt._numdiff.approx_derivative(fcn,x0,method='2-point',bounds=(lb,ub))
    J = np.atleast_2d(J)
    return J
#===============================================================================

#===============================================================================
def movmean(x, N):
    """
    Moving mean
    ===========

    Returns an array of local N-point mean values, where each mean is calculated over a sliding window of length k across neighboring elements of x.

    Usage:
    ------
        xfilt = movmean(x,N)

    Arguments:
    ----------
    x (array)
        Array to be filtered
    N (scalar)
        Window size
    
    Returns:
    --------
    xfilt (array)
        Filtered array
    
    """
    xfilt = np.convolve(x, np.ones(N)/N, mode='same')
    return xfilt
#===============================================================================

#===============================================================================
def ovl(A,B):
    """
    Overlap metric
    ==============

    Returns the overlap between two vectors A and B.

    Usage:
    ------
        metric = ovl(A,B)

    Arguments:
    ----------
    A (N-element array)
        First vector
    B (N-element array)
        Second vector
    
    Returns:
    --------
    metric (array)
        Overlap metric

    """
    A /= np.sum(A)
    B /= np.sum(B)
    metric = np.sum(np.minimum(A,B))
    return metric
#===============================================================================


def isempty(A):
#===============================================================================Q
    A = np.atleast_1d(A)
    boolean = np.size(A)==0
    return boolean
#===============================================================================


def multistarts(n,x0,lb,ub):
#===============================================================================

    if n<0:
        raise ValueError('The number of requested starting points must be n>0.') 

    if len(x0) != len(lb) or len(x0) != len(ub):
        raise ValueError('The lower/upper bound size(s) are not compatible with the initial guess vector x0.') 

    # Generate n-1 new starting points within the bounds
    if n>1:
        x0 = np.linspace(lb,ub,n-1)
    else:
        x0 = [x0]
    return x0
#===============================================================================
