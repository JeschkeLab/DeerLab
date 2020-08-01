# %% [markdown]
"""
Fitting a custom kernel model with a parameter-free distribution
=================================================================

How the use of SNLLS to fit a kernel model and a parameter-free 
distribution to a dipolar signal.
"""
import numpy as np
import matplotlib.pyplot as plt
from deerlab import *

# %% [markdown]
# Generating a dataset
#-----------------------------------------------------------------------------
# For this example we will simulate a simple 4pDEER signal

t = np.linspace(-0.5,5,300)
r = np.linspace(2,6,200)

# Generate ground truth and input signal
P = dd_gauss2(r,[3.5, 0.4, 0.4, 4.5, 0.7, 0.6])
lam = 0.36
c0 = 250 #uM
B = bg_hom3d(t,c0,lam)
K = dipolarkernel(t,r,lam,B)
V = K@P  + whitegaussnoise(t,0.01)

# %% [markdown]
# Fitting via SNLLS
#------------------
# Now in order to fit a non-linear dipolar kernel model ``Kmodel`` and a
# linear parameter-free distance distribution ``Pfit`` simultaneously, we
# can use the separable non-linear least squares ``SNLLS`` method. 
#
# First we define the function that contains the model for the dipolar kernel we want to fit. It 
# is a non-linear functon that accepts the parameter array ``p`` and returns the 
# fitted dipolar kernel ``K``. The linear parameters, in this case ``P``, are
# computed by solving a Tikhonov-regularized linear LSQ problem automatically in the ``snlls`` function. 

def Kmodel(p):

    # Unpack parameters
    lam,c0 = p
    # Get background
    B = bg_hom3d(t,c0,lam)
    # Generate 4pDEER kernel
    K = dipolarkernel(t,r,lam,B)

    return K

# %% [markdown]
# Next, there are two different parameter sets being fitted at the same time:
# linear and non-linear parameters. Therefore, the lower/upper bounds for
# the two sets need (or can) be specified.

#--------------------------
# Non-linear parameters:
#--------------------------
#       lam  c0
#--------------------------
par0 = [0.5, 50 ] # Start values
lb   = [ 0, 0.05] # lower bounds
ub   = [ 1, 1000] # upper bounds

#--------------------------
# Linear parameters: 
#--------------------------
#          Pfit
#--------------------------
lbl = np.zeros_like(r) # Non-negativity constraint of P
ubl = [] # Unconstrained upper boundary

# Run SNLLS optimization
parfit,Pfit,puq,stats = snlls(V,Kmodel,par0,lb,ub,lbl,ubl)

# Get non-linear parameters uncertainty
param95 = puq.ci(95,'nonlin')  #  95#-confidence interval

# Get linear parameters (distribution) uncertainty
Pci50 = puq.ci(50,'lin') #  50#-confidence interval
Pci95 = puq.ci(95,'lin') #  95#-confidence interval

# Print result
print('lambda = {:.2f}({:.2f}-{:.2f})'.format(parfit[0],param95[0,0],param95[0,1]))
print('c0 = {:.2f}({:.2f}-{:.2f})uM'.format(parfit[1],param95[1,0],param95[1,1]))

# Get fitted model
Kfit = Kmodel(parfit)
Vfit = Kfit@Pfit

# %% [markdown]
# Plots
#------

plt.subplot(211)
plt.plot(t,V,'k.',t,Vfit,'b')
plt.grid(alpha=0.3)
plt.xlabel('t [$\mu s$]')
plt.ylabel('V(t)')
plt.legend(['data','fit'])

plt.subplot(212)
plt.plot(r,P,'k',r,Pfit,'b')
plt.fill_between(r,Pci50[:,0],Pci50[:,1],color='b',alpha=0.4,linestyle='None')
plt.fill_between(r,Pci95[:,0],Pci95[:,1],color='b',alpha=0.2,linestyle='None')
plt.grid(alpha=0.3)
plt.xlabel('r [nm]')
plt.ylabel('P(r) [nm$^{-1}$]')
plt.legend(['truth','fit','50%-CI','95%-CI'])





# %%
