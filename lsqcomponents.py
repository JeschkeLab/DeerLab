import numpy as np

def lsqcomponents(V,K,L,alpha):

    Kt = K.T
    KtV = Kt@V
    KtK = Kt@K

    LtL = L.T@L;

    KtKreg = KtK + alpha**2*LtL

      
    return KtKreg, KtV

