# %% [markdown]
"""
Global fitting of a two-state model to several DEER traces
=================================================================

This example shows how to fit a two-state model to several DEER traces. Each of the
two states has a one-Gauss distance distribution, and each DEER trace comes from a
sample with different fractional populations of the two states, which could be the
consequence of a chemical or conformational equilibrium.
"""

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %% [markdown]
# Generate datasets
#-----------------------------------------------------------------------------
# For this example, we generate synthetic data.

# Parameters for the distance distibutions for states A and B
rmeanA = 3.45 # mean distance state A, in nm
rmeanB = 5.05 # mean distance state B, in nm
sigmaA = 0.3  # standard deviation state A, in nm
sigmaB = 0.2   # standard deviation state B, in nm
r = np.linspace(2,6,300)  # distance axis, in nm

fracA = [0.8, 0.5, 0.1] # Molar fraction of state A under several conditions

# Model parameters
par0 = [rmeanA, rmeanB, sigmaA, sigmaB] + fracA

# Generate list of time axes
tmin = [-0.2, -0.2, -0.2]   # start times, in µs
tmax = [4, 5, 4]   # end times, in µs
nt = [200, 150, 150] # number of time-domain points, in µs
N = len(tmin)
t = [np.linspace(tmin[i],tmax[i],nt[i]) for i in range(N)]

# Generate the corresponding dipolar kernels and signals
# (for the sake of simplicity no background and 100% modulation depth are assumed)
K = [dl.dipolarkernel(t[i],r) for i in range(N)]

# Model functions for V and P (needs K and r)
def myABmodel(par):

    # Unpack parameters
    rmeanA = par[0]
    rmeanB = par[1]
    sigmaA = par[2]
    sigmaB = par[3]
    fracA = par[4:7]
    
    N = len(fracA) # number of traces
    
    # Generate the state distributions
    PA = dl.dd_gauss(r,[rmeanA, sigmaA])
    PB = dl.dd_gauss(r,[rmeanB, sigmaB])
    
    # Generate distributions and signals for each sample
    P = [fracA[i]*PA+max(1-fracA[i],0)*PB for i in range(N)]
    V = [K[i]@P[i] for i in range(N)]

    return V, P

# Generate noise-free synthetic traces and distributions
V0, P0 = myABmodel(par0)

# Add noise
noiselevel = [0.05, 0.1, 0.1]
Vexp = [V0[i] + dl.whitegaussnoise(t[i],noiselevel[i],seed=i) for i in range(N)]


# %% [markdown]
# Global fit
# ----------
# Now when considering such systems is always important to (1) identify the
# parameters which must be fitted and (2) identify which parameters are the
# same for all signals (global) and which are specific for a individual
# signal (local). 
#
# In this examples we have the following parameters:
#   - global: ``rmeanA``, ``rmeanB`` (same for all samples)
#   - global: ``sigmaA``, ``sigmaB`` (same for all samples)
#   - Local:  ``fracA1``, ``fracA2`` (different from sample to sample)
#
# The next step is to construct the model function which describes our
# system. This function models the signals in our A-B system, and it is used to
# simulate all signals passed to ``fitparamodel``. The function must
# return (at least) a list of simulations of all the signals
# passed to ``fitparamodel``.

# Set starting values and bounds of fit parameters 
#        [rmeanA rmeanB sigmaA sigmaB fracA1 fracA2]
par0 =   [2,       2,     0.3,   0.3,   0.5,  0.5,  0.5]
lower =  [1,       1,     0.05,  0.05,  0,    0,     0]
upper =  [10,     10,     0.6,   0.6,   1,    1,     1]

# %% [markdown]
# Note that our model function ``myABmodel`` returns multiple outputs.
# This is advantegoud to later recover all fits directly, however, the fit
# function does only allow one output, specifically, the list of simulated signals.
# Therefore, we must create a lambda function which just takes the first ouput argument 
# of ``myABmodel``.

model = lambda par: myABmodel(par)[0] # call myABmodel and take the first output (V)

# Fit the global parametric model to both signals
fit = dl.fitparamodel(Vexp,model,par0,lower,upper,multistart=40)

# The use of the option 'multistart' will help the solver to find the
# global minimum and not to get stuck at local minima.

# Get the fitted models and distributions
Vfit, Pfit = myABmodel(fit.param)

# Get 95% confidence intervals of the fitted signals
Vfit_ci = [fit.modelUncert[i].ci(95) for i in range(N)]

# Propagate parameter uncertainty to the distribution models accounting for non-negativity
Pfit_uq = [fit.paramUncert.propagate(lambda param: myABmodel(param)[1][i],lbm=np.zeros_like(r)) for i in range(N)]

# Get their 95% confidence intervals 
Pfit_ci = [Pfit_uq[i].ci(95) for i in range(N)]

# %% [markdown]
# Plot results
# ------------

for i in range(N):
    plt.subplot(N,2,2*i+1)
    plt.plot(t[i],V0[i],'.',color='grey')
    plt.plot(t[i],Vfit[i],'tab:red')
    plt.fill_between(t[i],Vfit_ci[i][:,0],Vfit_ci[i][:,1],color='tab:red',alpha=0.3)
    plt.grid(alpha=0.3)
    plt.xlabel('t (µs)')
    plt.ylabel('V')
    plt.title(f'Trace {i+1}')

    plt.subplot(N,2,2*i+2)
    plt.plot(r,P0[i],'k',r,Pfit[i],'tab:red')
    plt.fill_between(r,Pfit_ci[i][:,0],Pfit_ci[i][:,1],color='tab:red',alpha=0.3)
    plt.grid(alpha=0.3)
    plt.xlabel('r (nm)')
    plt.ylabel('P (nm⁻¹)')
    plt.legend(['truth','fit'])

plt.tight_layout()
plt.show()
 # %%
