
#%%
import numpy as np
import time
import math
from deerlab import *
import matplotlib.pyplot as plt

def Kmodel(p,t,r):

    # Unpack parameters
    lam = p[0]
    k = p[1]

    # Get background
    B = bg_exp(t,k)

    # Generate 4pDEER kernel
    dr = r[2] - r[1]
    K = dipolarkernel(t,r)/dr
    K = (1-lam + lam*K)*B[:,np.newaxis]*dr
    return K


M = 20
N = np.linspace(80,600,M)
tocs = np.zeros(M)
for ii in range(len(N)):

    np.random.seed(1)
    t = np.linspace(-0.5,5,200)
    r = np.linspace(2,6,N[ii])

    # Generate ground truth and input signal
    P = dd_gauss2(r,[3.5, 0.4, 0.4, 4.5, 0.7, 0.6])
    lam = 0.36
    k = 0.5 #uM
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
    tic = time.clock()
    parfit, Pfit, uq = snlls(V,Amodel,par0,lb,ub,lbl,ubl)
    toc = time.clock()

    # Get fitted model
    Vfit = Kmodel(parfit,t,r)@Pfit
    Pci95 = uq.ci(95,'lin')
    Pci50 = uq.ci(50,'lin')

    tocs[ii] = toc-tic

    print('Fit: lambda=',parfit[0],', k =',parfit[1],'us-1')
    print('Processing time: ',toc-tic,'seconds')






# %%
