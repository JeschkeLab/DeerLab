
#%%
import numpy as np
import time
import math
from snlls import snlls
from deerlab import *
import matplotlib.pyplot as plt
import time


def Kmodel(p,t,r):

    # Unpack parameters
    lam = p[0]
    k = p[1]

    # Get background
    B = bg_exp(t,k)

    # Generate 4pDEER kernel
    dr = r[2] - r[1]
    K = dipolarkernel(t,r)/dr
    K = 1-lam + lam*K*B[:,np.newaxis]*dr
    return K

t = np.linspace(-0.5,5,300)
r = np.linspace(2,6,200)

# Generate ground truth and input signal
P = dd_gauss2(r,[3.5, 0.4, 0.4, 4.5, 0.7, 0.6])
lam = 0.36
k = 0.1 #uM
K = Kmodel([lam,k],t,r)
V = K@P + whitegaussnoise(t,0.01)

#--------------------------
# Non-linear parameters:
#--------------------------
#       lam  c0
#--------------------------
par0 = [0.5,  50 ] # Start values
lb   = [ 0 , 0.05] # lower bounds
ub   = [ 1 , 1000] # upper bounds

#--------------------------
# Linear parameters: 
#--------------------------
#          Pfit
#--------------------------
lbl = [0]*len(r) # Non-negativity constraint of P
ubl = []

Amodel = lambda p: Kmodel(p,t,r)

# Run SNLLS optimization
tic = time.clock()
parfit,Pfit = snlls(V,Amodel,par0,lb,ub,lbl,ubl)
toc = time.clock()

# Get fitted model
Vfit = Kmodel(parfit,t,r)@Pfit

print('Fit: lambda=',parfit[0],', k =',parfit[1],'us-1')
print('Processing time: ',toc-tic,'seconds')
plt.figure(1)
plt.subplot(2,1,1)
plt.plot(t,V,'k.',t,Vfit,'b')
plt.subplot(2,1,2)
plt.plot(r,P,'k',r,Pfit,'b')
plt.show



# %%
