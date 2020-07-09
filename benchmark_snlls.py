
#%%
import numpy as np
import time
import math
from deerlab import *
import matplotlib.pyplot as plt
import scipy
def Kmodel(p,t,r):

    # Unpack parameters
    lam = p[0]
    k = p[1]

    # Get background
    B = bg_exp(t,k)

    # Generate 4pDEER kernel
    K = dipolarkernel(t,r,lam,B)
    return K

M = 20
N = np.linspace(80,1000,M)
tocs1 = np.zeros(M)
tocs2 = np.zeros(M)
tocs3 = np.zeros(M)
for ii in range(len(N)):

    np.random.seed(1)
    t = np.linspace(-0.5,5,N[ii])
    r = np.linspace(2,6,N[ii])

    # Generate ground truth and input signal
    P = dd_gauss2(r,[3.5, 0.4, 0.4, 4.5, 0.7, 0.6])
    lam = 0.36
    k = 0.15 #uM
    K = Kmodel([lam,k],t,r)
    V = K@P + whitegaussnoise(t,0.01)

    #--------------------------
    # Non-linear parameters:
    #--------------------------
    #       lam  c0
    #--------------------------
    par0 = [0.5, 0.5 ] # Start values
    lb   = [ 0 , 0.05] # lower bounds
    ub   = [ 1 ,  1  ] # upper bounds

    #--------------------------
    # Linear parameters: 
    #--------------------------
    #          Pfit
    #--------------------------
    lbl = np.zeros(len(r)) # Non-negativity constraint of P
    ubl = []

    Amodel = lambda p: Kmodel(p,t,r)

    # Run SNLLS optimization
    #tic = time.clock()
    #parfit1, Pfit1, uq = snlls(V,Amodel,par0,lb,ub,lbl,ubl,linsolver='fnnls')
    #toc = time.clock()
    #tocs1[ii] = toc-tic
    #print('Processing time fnnls: ',toc-tic,'seconds')

    #tic = time.clock()
    #parfit2, Pfit2, uq = snlls(V,Amodel,par0,lb,ub,lbl,ubl,linsolver='nnls')
    #toc = time.clock()
    #tocs2[ii] = toc-tic
    #print('Processing time nnls: ',toc-tic,'seconds')

    for i in range(5):
        tic = time.clock()
        parfit3, Pfit3, uq = snlls(V,Amodel,par0,lb,ub,lbl,ubl,linsolver='cvx')
        toc = time.clock()
        tocs3[ii] = tocs3[ii] + (toc-tic)
    tocs3[ii] = tocs3[ii]/5
    print('Processing time cvx: ',toc-tic,'seconds')

# %%
