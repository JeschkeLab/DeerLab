#%%
import numpy as np
import time
import math
from deerlab import *
import matplotlib.pyplot as plt

# No randomness
np.random.seed(1)


M = 20
N = np.linspace(80,600,M)
tocs = np.zeros(M)
for ii in range(len(N)):

    tic = time.clock()

    # Simulation
    t = np.linspace(-0.5,5,300)
    r = np.linspace(2,6,N[ii])
    P = dd_gauss2(r,[4, 0.5, 0.4, 5, 0.6, 0.6])
    K = dipolarkernel(t,r)
    V = K@P + whitegaussnoise(t,0.00)

    # Fitting
    alpha = selregparam(V,K,r)
    Pfit  = fitregmodel(V,K,r,alpha)
    Vfit = K@Pfit

    toc = time.clock()
    tocs[ii] = toc - tic
    print(toc - tic, 'seconds')


# %%
plt.plot(N,tocs)

# %%
