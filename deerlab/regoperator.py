import numpy as np


def regoperator(r,n=2):
    N = np.size(r)
    L = np.eye(N,N+n)
    L = np.diff(L,n)
    L = L.T
    L = L[0:len(r)-n,:]
    return L
