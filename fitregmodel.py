import numpy as np
from fnnls import fnnls
from scipy.optimize import nnls
from regoperator import regoperator

def fitregmodel(V,K,r,alpha=1,regorder=2, solver='fnnls'):
    L = regoperator(r,regorder)

    if solver == 'fnnls':
        KtKreg = np.transpose(K)@K + alpha**2 * np.transpose(L)@L
        KtV = np.transpose(K)@V
        Pfit = fnnls(KtKreg,KtV)
    elif solver == 'nnls':
        C = np.concatenate([K, alpha * L])
        d = np.concatenate([V, np.zeros(K.shape[1])])
        Pfit, _ = nnls(C, d)
    else:
        raise KeyError(f'{solver} is not a known non-negative least squares solver')


    return Pfit
