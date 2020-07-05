#%%
import numpy as np
import time
import math
from deerlab import *
import matplotlib.pyplot as plt
import matrixtools as mat

# No randomness
np.random.seed(1)

tic = time.clock()

# Simulation
t = np.linspace(-0.5,5,300)
r = np.linspace(2,6,300)
P = dd_gauss2(r,[4, 0.5, 0.4, 5, 0.6, 0.6])
K = dipolarkernel(t,r)

print(np.shape(P))
print(np.shape(K))

V = K@P + whitegaussnoise(t,0.01)

# Fitting
alpha = selregparam(V,K,r)
Pfit  = fitregmodel(V,K,r,alpha)
Vfit = K@Pfit

toc = time.clock()
print(toc - tic, 'seconds')

plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t,V,'k.',t,Vfit,'b')
plt.subplot(2,1,2)
plt.plot(r,P,'k',r,Pfit,'b')
plt.show


# %%
