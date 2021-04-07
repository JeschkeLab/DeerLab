# %% [markdown]
"""
Global model fits with global, local and fixed parameters
=========================================================

This example shows how to fit multiple signals to a global model, which
may depend on some parameters which need to be globally fitted, some
locally and some might be fixed and not fitted. 

"""

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %% [markdown]
# Generate two datasets
#-----------------------------------------------------------------------------
# For this example we will simulate a system containing two states A and B
# both havng a Gaussian distribution of known width but unknown mean
# distance. For this system we have two measurements V1 amd V2 measured
# under two different conditions leading to different fractions of states A
# and B. 

r = np.linspace(2,6,300)  # distance axis, in nm
t1 = np.linspace(0,4,200) # time axis of first measurement, in µs
t2 = np.linspace(0,6,150) # time axis of first measurement, in µs

# Parameters
rmeanA = 3.45 # mean distance state A, in nm
rmeanB = 5.05 # mean distance state B, in nm
sigmaA = 0.3  # standard deviation state A, in nm
sigmaB = 0.2   # standard deviation state B, in nm

fracA1 = 0.8 # Molar fraction of state A under conditions 1
fracA2 = 0.2 # Molar fraction of state A under conditions 2
# The molar fraction of state B is not required as it follows fracB = 1 - fracA

# Generate the two distributions for conditions 1 & 2
P1 = dl.dd_gauss2(r,[rmeanA, sigmaA, fracA1, rmeanB, sigmaB, 1-fracA1])
P2 = dl.dd_gauss2(r,[rmeanA, sigmaA, fracA2, rmeanB, sigmaB, 1-fracA2])

# Generate the corresponding dipolar kernels
K1 = dl.dipolarkernel(t1,r)
K2 = dl.dipolarkernel(t2,r)

# ...and the two corresponding signals
V1 = K1@P1 + dl.whitegaussnoise(t1,0.05,seed=0)
V2 = K2@P2 + dl.whitegaussnoise(t2,0.1,seed=1)
# (for the sake of simplicity no background and 100# modulation depth are assumed)

# %% [markdown]
# Global fit
# ----------
# Now when considering such systems is always important to (1) identify the
# parameters which must be fitted and (2) identify which parameters are the
# same for all signals (global) and which are specific for a individual
# signal (local). 
#
# In this examples we have the following parameters:
#   - fixed: ``sigmaA``, ``sigmaB`` (known paramters)
#   - global: ``rmeanA``, ``rmeanB`` (same for both signals)
#   - local: ``fracA1``, ``fracA2`` (different for both signals/conditions)
#
# The next step is to construct the model function which describes our
# system. This function models the signals in our A-B system, and it is used to
# simulate all signals passed to ``fitparamodel``. The function must
# return (at least) a list of simulations of all the signals
# passed to ``fitparamodel``.

# Model definition
def myABmodel(par):

    #Fixed parameters
    sigmaA = 0.5
    sigmaB = 0.3
    #Global parameters
    rmeanA = par[0]
    rmeanB = par[1]
    #Local parameters
    fracA1 = par[2]
    fracA2 = par[3]
    
    # Generate the signal-specific distribution
    Pfit1 = dl.dd_gauss2(r,[rmeanA, sigmaA, fracA1, rmeanB, sigmaB, max(1-fracA1,0)])
    Pfit2 = dl.dd_gauss2(r,[rmeanA, sigmaA, fracA2, rmeanB, sigmaB, max(1-fracA2,0)])

    # Generate signal #1
    V1fit = K1 @ Pfit1
    # Generate signal #2
    V2fit = K2 @ Pfit2
    # Return as a list
    Vfits = [V1fit,V2fit]

    return Vfits,Pfit1,Pfit2

#-----------------------------------------
#                Fit parameters 
#-----------------------------------------
#        [rmeanA rmeanB fracA1 fracA2]
#-----------------------------------------
par0 =   [2,       2,    0.5,    0.5]
lower =  [1,       1,     0,      0]
upper =  [20,     20,     1,      1]
#-----------------------------------------

# %% [markdown]
# Note that our model function ``myABmodel`` returns multiple outputs.
# This is advantegoud to later recover all fits directly, however, the fit
# function does only allow one output, specifically, the list of simulated signals.
# Therefore, we must create a lambda function which just takes the first ouput argument 
# of ``myABmodel``.

model = lambda par: myABmodel(par)[0] # call myABmodel with par and take the first output

# Collect data for global fit into cell arrays
Vs = [V1,V2]

# Fit the global parametric model to both signals
fit = dl.fitparamodel(Vs,model,par0,lower,upper,multistart=40)

# The use of the option 'multistart' will help the solver to find the
# global minimum and not to get stuck at local minima.

# Define individual models for each fitted distribution
P1_model = lambda param: myABmodel(param)[1]
P2_model = lambda param: myABmodel(param)[2]

# Get the fitted models 
Vfit1 = fit.model[0]
Vfit2 = fit.model[1]
Pfit1 = P1_model(fit.param)
Pfit2 = P2_model(fit.param)

# Get uncertainties of the fitted signals
Vfit1_uq = fit.modelUncert[0]
Vfit2_uq = fit.modelUncert[1]
# Propagate parameter uncertainty to the distribution models accounting for non-negativity
Pfit1_uq = fit.paramUncert.propagate(P1_model,lbm=np.zeros_like(r))
Pfit2_uq = fit.paramUncert.propagate(P2_model,lbm=np.zeros_like(r))

# Get their 95%-confidence intervals 
Vfit1_ci = Vfit1_uq.ci(95)
Vfit2_ci = Vfit2_uq.ci(95)
Pfit1_ci = Pfit1_uq.ci(95)
Pfit2_ci = Pfit2_uq.ci(95)

# %% [markdown]
# Plot results
# ------------
plt.subplot(221)
plt.plot(t1,V1,'.',color='grey')
plt.plot(t1,Vfit1,'r')
plt.fill_between(t1,Vfit1_ci[:,0],Vfit1_ci[:,1],color='r',alpha=0.3)
plt.grid(alpha=0.3)
plt.xlabel('t (µs)')
plt.ylabel('V')
plt.title('Conditions #1')

plt.subplot(222)
plt.plot(r,P1,'k',r,Pfit1,'r')
plt.fill_between(r,Pfit1_ci[:,0],Pfit1_ci[:,1],color='r',alpha=0.3)
plt.grid(alpha=0.3)
plt.xlabel('r (nm)')
plt.ylabel('P (nm⁻¹)')
plt.legend(['truth','fit'])

plt.subplot(223)
plt.plot(t2,V2,'.',color='grey')
plt.plot(t2,Vfit2,'b')
plt.fill_between(t2,Vfit2_ci[:,0],Vfit2_ci[:,1],color='b',alpha=0.3)
plt.grid(alpha=0.3)
plt.xlabel('t (µs)')
plt.ylabel('V')
plt.title('Conditions #2')

plt.subplot(224)
plt.plot(r,P2,'k',r,Pfit2,'b')
plt.fill_between(r,Pfit2_ci[:,0],Pfit2_ci[:,1],color='b',alpha=0.3)
plt.grid(alpha=0.3)
plt.xlabel('r (nm)')
plt.ylabel('P (nm⁻¹)')
plt.legend(['truth','fit'])

plt.tight_layout()
plt.show()
 # %%
