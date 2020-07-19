from builtins import str

def isstring(x):
    bool = isinstance(x,str)
    return bool


import warnings
import matplotlib.pyplot as plt
import numpy as np
import cmath as math
import scipy as scp
import random
from scipy.sparse import coo_matrix

def parse_multidatasets(V,K,weights):
#===============================================================================
    # If multiple signals are specified as a list...
    if type(V) is list and all([type(Vs) is np.ndarray for Vs in V]):
        nSignals = len(V)
        Vlist = V
        V = np.concatenate(V, axis=0) # ...concatenate them along the list 
    elif type(V) is np.ndarray:
        nSignals = 1
    else:
        raise TypeError('The input signal(s) must be numpy array or a list of numpy arrays.')

    # If multiple kernels are specified as a list...
    if type(K) is list and all([type(Ks) is np.ndarray for Ks in K]):
        nKernels = len(K)
        K = np.concatenate(K, axis=0) # ...concatenate them along the list 
    elif type(K) is np.ndarray:
        nKernels = 1
    else:
        raise TypeError('The input kernel(s) must be numpy array or a list of numpy arrays.')

    # If multiple weights are specified as a list...
    if type(weights) is list and all([not hasattr(w, "__len__") for w in weights]):
        weights = np.asarray(weights)
        weights = weights/sum(weights)
        if len(weights)!=nSignals:
            raise KeyError('If multiple signals are passed, the same number of weights are required.')
        weights_ = []
        for i in range(len(weights)):
            weights_ = np.concatenate((weights_,weights[i]*np.ones(len(Vlist[i]))))
        weights = weights_
    elif not hasattr(weights, "__len__"):
        if weights!=1:
            raise ValueError('If multiple signals are passed, the same number of weights are required.')
    else:
        raise TypeError('The input signal(s) must be numpy array or a list of numpy arrays.')

    if nSignals!=nKernels:
        raise KeyError('The same number of kernels and signals must be specified as lists.')

    return V,K,weights
#===============================================================================


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
    U,_,_,C,S = csd(Q1,Q2)

    # Vector of generalized singular values.
    q = min(m+n,p)
    # Supress divide by 0 warning        
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        U = np.vstack((np.zeros((q-m,1),'double'), np.diag(C,max(0,q-m)).reshape(len(np.diag(C,max(0,q-m))),1))) /  np.vstack((np.diag(S,0).reshape(len(np.diag(S,0)),1), np.zeros((q-n,1),'double') ))


    return U
#===============================================================================


def csd(Q1,Q2):
#===============================================================================
    """
    Cosine-Sine Decomposition

    U,V,Z,C,S = csd(Q1,Q2)

    Given Q1 and Q2 such that Q1'@Q1 + Q2'@Q2 = I, the
    C-S Decomposition is a joint factorization of the form
    Q1 = U@C@Z' and Q2=V@S@Z'
    where U, V, and Z are orthogonal matrices and C and S
    are diagonal matrices (not necessarily square) satisfying
    C'@C + S'@S = I
    """
    [m,p] = np.shape(Q1)
    n = np.shape(Q2)[0]


    if m < n:
        V,U,Z,S,C = csd(Q2,Q1)
        j = np.flip(np.arange(p)) 
        C = C[:,j] 
        S = S[:,j] 
        Z = Z[:,j]
        m = min(m,p) 
        i =  np.flip(np.arange(m)) 
        C[np.arange(m),:] = C[i,:] 
        U[:,np.arange(m)] = U[:,i]
        n = min(n,p) 
        i =  np.flip(np.arange(n)) 
        S[np.arange(n),:] = S[i,:] 
        V[:,np.arange(n)] = V[:,i]
        return U,V,Z,C,S
    
    U,sdiag,VH = np.linalg.svd(Q1)
    C= np.zeros((m, p))
    np.fill_diagonal(C, sdiag)
    Z = VH.T.conj()

    q = min(m,p)
    i = np.arange(0, q, 1)
    j = np.arange(q-1, -1, -1)
    C[i,i] = C[j,j]
    U[:,i] = U[:,j]
    Z[:,i] = Z[:,j]
    S = Q2@Z

    if q == 1:
        k = 0
    elif m < p:
        k = n
    else:
        k = max(0, np.max(np.where(np.diag(C) <= 1/math.sqrt(2))))
    V,_ = scp.linalg.qr(S[:,0:k])

    S = V.T@S
    r = min(k,m)
    S[:,0:r-1] = diagf(S[:,0:r-1])

    if (m == 1 and p > 1):
        S[0,0] = 0

    if k < min(n,p):
        r = min(n,p)
        i = np.arange(k, n, 1)
        j = np.arange(k, r, 1)
        [UT,STdiag,VH] = np.linalg.svd(S[np.ix_(i,j)])
        ST = np.zeros((len(i), len(j)))
        np.fill_diagonal(ST, STdiag)
        VT = VH.T.conj()
        if k > 0:
            S[0:k,j] = 0
        S[np.ix_(i,j)] = ST
        C[:,j] = C[:,j]@VT
        V[:,i] = V[:,i]@UT
        Z[:,j] = Z[:,j]@VT
        i = np.arange(k, q, 1)
        [Q,R] = scp.linalg.qr(C[np.ix_(i,j)])
        C[np.ix_(i,j)] = diagf(R)
        U[:,i] = U[:,i]@Q

    if m < p:
        # Diagonalize final block of S and permute blocks.
        eps = np.finfo(float).eps
        q = min([np.count_nonzero(abs(np.diag(C,0))>10*m*eps),
                np.count_nonzero(abs(np.diag(S,0))>10*n*eps),
                np.count_nonzero(np.amax(abs(S[:,m+1:p]),axis = 1)<np.sqrt(eps))])

        # maxq: maximum size of q such that the expression used later on,
        #        i = [q+1:q+p-m, 1:q, q+p-m+1:n],
        # is still a valid permutation.
        maxq = m + n - p
        q = q + np.count_nonzero(np.amax(abs(S[:,q+1:maxq+1]),axis = 1)>np.sqrt(eps))

        i = np.arange(q,n,1)
        j = np.arange(m,p,1)
        # At this point, S(i,j) should have orthogonal columns and the
        # elements of S(:,q+1:p) outside of S(i,j) should be negligible.
        Q,R = scp.linalg.qr(S[np.ix_(i,j)])
        S[:,q+1:p] = 0
        S[np.ix_(i,j)] = diagf(R)
        V[:,i] = V[:,i]@Q
        if n > 1:
            i = np.concatenate((np.arange(q,q+p-m,1), np.arange(0,q,1), np.arange(q+p-m,n,1)))
        else:
            i = 1
        j = np.concatenate((np.arange(m,p,1), np.arange(0,m,1)))
        C = C[:,j]
        S = S[np.ix_(i,j)]
        Z = Z[:,j]
        V = V[:,i]

    if n < p:
        #Final block of S is negligible.
        S[:,n+1:p] = 0

    # Make sure C and S are real and positive.
    U,C = diagp(U,C,max(0,p-m))
    C = C.real

    V,S = diagp(V,S,0)
    S = S.real

    return U,V,Z,C,S
#===============================================================================

def diagf(X):
#===============================================================================
    """
    Diagonal force

    X = diagf(X) zeros all the elements off the main diagonal of X.
    """
    X = np.triu(np.tril(X))
    return X
#===============================================================================


def diagp(Y,X,k):
#===============================================================================
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
