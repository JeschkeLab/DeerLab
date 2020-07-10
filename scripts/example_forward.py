#%%
import numpy as np
import time
import math
import deerlab as dl
import matplotlib.pyplot as plt

# No randomness
np.random.seed(1)

tic = time.clock()
# Simulation
t = np.linspace(-0.5,5,300)
r = np.linspace(2,6,300)
P = dl.dd_gauss2(r,[4, 0.5, 0.4, 5, 0.6, 0.6])
K = dl.dipolarkernel(t,r)
V = K@P + dl.whitegaussnoise(t,0.00)

# Fitting
alpha = dl.selregparam(V,K,r)
Pfit  = dl.fitregmodel(V,K,r,alpha)
Vfit = K@Pfit

toc = time.clock()
print(toc - tic, 'seconds')
np.argmax
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t,V,'k.',t,Vfit,'b')
plt.subplot(2,1,2)
plt.plot(r,P,'k',r,Pfit,'b')
plt.show


# %%
import deerlab.test 

deerlab.test.runall()

# %%
