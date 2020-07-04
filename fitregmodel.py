import numpy as np
from fnnls import fnnls
from regoperator import regoperator

def fitregmodel(V,K,r,alpha=1,regorder=2):
    L = regoperator(r,regorder)
    KtKreg = np.transpose(K)@K + alpha**2 * np.transpose(L)@L
    KtV = np.transpose(K)@V
    Pfit = fnnls(KtKreg,KtV)
    return Pfit
