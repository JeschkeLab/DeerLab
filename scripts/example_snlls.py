
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
    K = dipolarkernel(t,r,lam,B)
    return K

t = np.linspace(-0.5,5,400)
r = np.linspace(2,6,400)

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
tic = time.clock()
parfit, Pfit, uq = snlls(V,Amodel,par0,lb,ub,lbl,ubl,linsolver='cvx')
toc = time.clock()

# Get fitted model
Vfit = Kmodel(parfit,t,r)@Pfit
Pci95 = uq.ci(95,'lin')
Pci50 = uq.ci(50,'lin')

print('Fit: lambda=',parfit[0],', k =',parfit[1],'us-1')
print('Processing time: ',toc-tic,'seconds')

plt.figure(1)

plt.subplot(2,1,1)
plt.plot(t,V,'k.',t,Vfit,'b')
plt.tight_layout()
plt.grid('on',alpha=0.4)
plt.xlabel('t [us]')
plt.ylabel('V(t)')
plt.legend(['data','fit'])

plt.subplot(2,1,2)
plt.plot(r,P,'k',r,Pfit,'b')
plt.fill_between(r,Pci95[:,0],Pci95[:,1],color='b',alpha=0.2)
plt.fill_between(r,Pci50[:,0],Pci50[:,1],color='b',alpha=0.3)
plt.tight_layout()
plt.grid('on',alpha=0.4)

plt.xlabel('r [nm]')
plt.ylabel('P(r) [nm^-1]')
plt.legend(['data','fit','CI'])




# %%
