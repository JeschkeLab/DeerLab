# %% [markdown] 
"""
Comparing confidence intervals for regularization results
=========================================================

A simpe example of uncertainty estimation for Tikhonov regularization 
results. The example will cover the use of confidence intervals
obtained from curvature matrices and boostrap analysis.
""" 

import numpy as np 
import matplotlib.pyplot as plt
from deerlab import *

# %% [markdown] 
# Simulate the data
# ------------------
#
# Let's start by generating some data.


# Prepare signal components
t = np.linspace(-0.4,3.5,200)

# Use a distance-axis with less points to make analysis faster
r = np.linspace(2,5,200)

P = dd_gauss2(r,[3, 0.1, 0.6, 3.5, 0.2, 0.4])
B = bg_strexp(t,[0.04,1])
lam = 0.32

# Generate the dipolar kernel
K = dipolarkernel(t,r,lam,B)
# Simulate signal
V = K@P + whitegaussnoise(t,0.01)

# %% [markdown] 
# For the sake of simplicity, in this examples we will assume that we know the 
# background exactly. Our first step is to generate the proper dipolar kernel.
#
# Covariance-based confidence intervals
# -------------------------------------
#
# We now have all the elements required to fit our distance distribution via 
# regularization. We will use the AIC to select the regularization parameter in 
# the Tikhonov regularization.

# Fit data via regularization
fit = fitregmodel(V,K,r,'tikhonov','aic')
Pfit = fit.P
Puq = fit.uncertainty
# Obtain time-domain fit
Vfit = K@Pfit

plt.plot(r,P,'k',r,Pfit,'r',linewidth=1)
Pci95 = Puq.ci(95)
Pci50 = Puq.ci(50)
plt.fill_between(r,Pci50[:,0],Pci50[:,1],color='r',linestyle='None',alpha=0.45)
plt.fill_between(r,Pci95[:,0],Pci95[:,1],color='r',linestyle='None',alpha=0.25)
plt.grid()
plt.xlabel('r [nm]')
plt.ylabel('P(r) [nm$^{-1}$]')
plt.title('Curvature Matrix CI')
plt.legend(['Truth','Fit','50%-CI','95%-CI'])

# %% [markdown]
# Bootstrapped confidence intervals
# ---------------------------------
#
# Now we are interested in the bootstrap confidence intervals. For this, we
# need to define a boot function e.g. ``mybootfcn()`` which takes a signal as
# output and returns the outputs of interest (``Pfit`` in our example).


def mybootfcn(V):
    fit = fitregmodel(V,K,r,'tikhonov','aic')
    return fit.P

# Launch bootstrapping
Nsamples = 50
booci = bootan(mybootfcn,V,Vfit,Nsamples)
Pci95 = booci.ci(95)
Pci50 = booci.ci(50)

# %% [markdown]
# By plotting the results, one can see that the bootstrapped confidence intervals 
# are narrower in comparison to the ones obtained via the curvature
# matrices. This is due to the inherent accurate nature of bootstrapping. 

plt.plot(r,P,'k',r,Pfit,'b',linewidth=1)
plt.fill_between(r,Pci50[:,0],Pci50[:,1],color='b',linestyle='None',alpha=0.45)
plt.fill_between(r,Pci95[:,0],Pci95[:,1],color='b',linestyle='None',alpha=0.25)
plt.grid(alpha=0.3)
plt.xlabel('r [nm]')
plt.ylabel('P(r) [nm$^{-1}$]')
plt.title('Bootstrapped CI')
plt.legend(['Truth','Fit','50%-CI','95%-CI'])


# %%
