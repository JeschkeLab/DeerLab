
# %% [markdown]
"""
Uncertainty propagation from parameter fits using covariance-based uncertainty quantification
=============================================================================================

How to propagate the uncertainty of the fitted parameters to the models which depend on them.
"""

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %% [markdown]
# Generate data
# -------------

t = np.linspace(-0.2,4,300)   # µs
r = np.linspace(2,5,400)     # nm
center = 3.5 # Rician center distance, nm
width = 0.3 # Rician width, nm
lam = 0.27 # Modulation depth
conc = 150 # Spin concentration, µM
P = dl.dd_rice(r,[center, width])
B = dl.bg_hom3d(t,conc,lam)
K = dl.dipolarkernel(t,r,mod=lam,bg=B)
V = K@P + dl.whitegaussnoise(t,0.03,seed=0)

# %% [markdown]
# Fit the data
# ------------
# First we define the models for the different elements in our analysis
# (background, distribution and dipolar signal). For simplicity these
# models take the full parameter set
#
# ``par = [lambda center width conc]``
#
# and select the appropiate elements from the parameter set, i.e.
#
# ``Pmodel = f(center,width) -> par[1] & par[2]``
# ``Bmodel = f(conc,lambda)  -> par[3] & par[0]``
# ``Vmodel = f(par)          -> par[0] & par[1] & par[2] & par[3]``
#
# By defining the models like this, we can spare then the indexing of the
# parameters each time we call one of these model and can pass the full
# parameter set directly.

# Pre-calculate the elemental dipolar kernel (for speed)
K0 = dl.dipolarkernel(t,r)

Pmodel = lambda par: dl.dd_rice(r,par[1:3])
Bmodel = lambda par: dl.bg_hom3d(t,par[3],par[0])
Vmodel = lambda par: (1 - par[0] + par[0]*K0@Pmodel(par))*Bmodel(par)

# %% [markdown]
# Next since we are dealing with a custom-defined model we need to specify
# the start values as well as boundaries of the parameter set:

# Parameters:[lam center width conc]
par0  =      [0.35, 4.0,  0.4, 500 ] # start values
lower =      [0.10, 2.0,  0.1, 0.1 ] # lower bounds
upper =      [0.50, 7.0,  0.5, 1500] # upper bounds

# Finally we can run the fit and get the fitted parameters and their uncertainties
fit = dl.fitparamodel(V,Vmodel,par0,lower,upper)

parfit = fit.param
paruq = fit.paramUncert

# Forward-calculate the models with the fitted parameters
Vfit = Vmodel(parfit)
Pfit = Pmodel(parfit)
Bfit = Bmodel(parfit)
lamfit = parfit[0]

# %% [markdown]
# Uncertainty propagation
#------------------------
# In DeerLab, all uncertainty quantification objects contain a method
# ``.propagate()``, which has all the internal information on the 
# covariance matrices required to propagate the uncertainty from 
# the parameters to the fits. 
#
# Thus, all we neeed to do is call ``.propagate``` and pass the model function
# which we want to propagate the uncertainty to. It is important that if
# the uncertainty quantification structure is defined for N-parameters (N=4
# in this case) the model function must accept all N parameters. Since we
# defined our model function to accept all N parameters already we do not
# need to worry about it.

# %% [markdown]
# 1. Uncertainty of the dipolar signal fit: This case is easy, we already have the model and it is unconstrained
Vuq = paruq.propagate(Vmodel) # Uncertainty quantification for Vfit
Vci95 = Vuq.ci(95) # 95#-confidence intervals for Vfit

# %% [markdown]
# 2. Uncertainty of the distance distribution: In this case, the distribution has a non-negativity constraint which we
# can specify via the lb input. 
lb = np.zeros_like(r) # Non-negativity constraint
Puq = paruq.propagate(Pmodel,lb) # Uncertainty quantification for Pfit
Pci95 = Puq.ci(95) # 95#-confidence intervals for Pfit

# %% [markdown]
# 3. Uncertainty of the background: In this case, since we want to use this for plotting we need to evaluate
# the function (1-lambda)*Bfit instead of just Bfit in order to plot the\
# correct function.
Buq = paruq.propagate(lambda p:(1-p[0])*Bmodel(p)) # Uncertainty quantification for (1-lam)Bfit
Bci95 = Buq.ci(95) # 95#-confidence intervals for (1-lam)Bfit

# %% [markdown]
# Plots
# -----

plt.figure(figsize=(7,7))

# Time-domain
plt.subplot(211)
plt.plot(t,V,'k.',t,Vfit,'r',t,(1-lamfit)*Bfit,'b',linewidth=1.5)
plt.fill_between(t,Vci95[:,0],Vci95[:,1],color='r',alpha=0.3,linestyle='None')
plt.fill_between(t,Bci95[:,0],Bci95[:,1],color='b',alpha=0.3,linestyle='None')
plt.grid(alpha=0.3)
plt.xlabel('t (µs)')
plt.ylabel('V')
plt.legend(['data','Vfit','Bfit','Vfit 95%-CI','Bfit 95%-CI'])

# Distance-domain
plt.subplot(212)
plt.plot(r,P,'k',r,Pfit,'r',linewidth=1.5)
plt.fill_between(r,Pci95[:,0],Pci95[:,1],color='r',alpha=0.3,linestyle='None')
plt.xlabel('r (nm)')
plt.ylabel('P (nm⁻¹)')
plt.grid(alpha=0.3)
plt.legend(['truth','Pfit','Pfit 95%-CI'])

plt.show()

# %%


# %%
