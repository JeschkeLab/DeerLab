import numpy as np
from deerlab import fnnls
from scipy.optimize import nnls
import deerlab as dl

def fitregmodel(V,K,r,alpha=1,regorder=2, solver='fnnls'):
    L = dl.regoperator(r,regorder)

    if solver == 'fnnls':
        KtKreg = np.transpose(K)@K + alpha**2 * np.transpose(L)@L
        KtV = np.transpose(K)@V
        Pfit = fnnls(KtKreg,KtV)
    elif solver == 'cvx':
        KtKreg = np.transpose(K)@K + alpha**2 * np.transpose(L)@L
        KtV = np.transpose(K)@V
        Pfit = dl.cvxnnls(KtKreg, KtV, K, V)
    elif solver == 'nnls':
        C = np.concatenate([K, alpha * L])
        d = np.concatenate([V, np.zeros(K.shape[1])])
        Pfit, _ = nnls(C, d)
    else:
        raise KeyError(f'{solver} is not a known non-negative least squares solver')


    return Pfit
