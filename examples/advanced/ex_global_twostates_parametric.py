# %% [markdown]
"""
Global fitting of a two-state model to a series of DEER traces
=================================================================

This example shows how to fit a two-state model to a series DEER traces. Each of the
two states, A and B, has a one-Gauss distance distribution, and each DEER trace comes
from a sample with different fractional populations of the two states. This could be
the consequence of a chemical or conformational equilibrium. The model contains global
parameters needed for all samples traces (the distribution parameters) and local
parameters needed for individual samples traces (the fractional populations).
"""

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %% [markdown]
# Generate datasets
#-----------------------------------------------------------------------------
# For this example, we generate synthetic data.

# Parameters for the distance distibutions for states A and B
rmeanA = 3.45  # mean distance state A, in nm
rmeanB = 5.05  # mean distance state B, in nm
sigmaA = 0.3   # standard deviation state A, in nm
sigmaB = 0.2   # standard deviation state B, in nm
r = np.linspace(2,6,300)  # distance axis, in nm

fracA = [0.8, 0.5, 0.1] # molar fraction of state A for each sample

# Model parameters
par0 = [rmeanA, rmeanB, sigmaA, sigmaB] + fracA

# Generate list of time axes
tmin = [-0.2, -0.2, -0.2]   # start times, in µs
tmax = [4, 5, 4]            # end times, in µs
nt = [200, 150, 150]        # number of time-domain points
N = len(tmin)               # number of DEER traces
t = [np.linspace(tmin[i],tmax[i],nt[i]) for i in range(N)]

# Generate the corresponding dipolar kernels
# (for the sake of simplicity, no background and 100% modulation depth are assumed)
K = [dl.dipolarkernel(t_,r) for t_ in t]

# Model functions for V and P (needs K and r)
def ABmodel(par):

    # Unpack parameters
    rmeanA, rmeanB, sigmaA, sigmaB, *fracA = par
    
    N = len(fracA) # number of traces
    
    # Generate the state distributions
    PA = dl.dd_gauss(r,[rmeanA, sigmaA])
    PB = dl.dd_gauss(r,[rmeanB, sigmaB])
    
    # Generate distributions and signals for each sample
    P = [fA*PA+(1-fA)*PB for fA in fracA]
    V = [K[i]@P[i] for i in range(N)]

    return V, P

# Generate noise-free synthetic traces and distributions
V0, P0 = ABmodel(par0)

# Add noise
noiselevel = [0.05, 0.1, 0.1]
Vexp = [V0[i] + dl.whitegaussnoise(t[i],noiselevel[i],seed=i) for i in range(N)]


# %% [markdown]
# Global fit
# ----------
# We now want to fit the model to the generated data. The fit parameters are
# the distribution parameters for the two states (rmeanA, rmeanB, sigmaA, sigmaB)
# and the fractional population of state A in each sample (fracA)

# Set starting values and bounds of fit parameters 
#        [rmeanA rmeanB sigmaA sigmaB fracA1 fracA2 fracA3]
par0 =   [2,       2,     0.3,   0.3,   0.5,  0.5,  0.5]
lower =  [1,       1,     0.05,  0.05,  0,    0,     0]
upper =  [10,     10,     0.6,   0.6,   1,    1,     1]

# %% [markdown]
# Out model function ``ABmodel`` returns multiple outputs, V and P.
# This is useful for later obtaining boht V and P, but the fit function requires
# just one model output. Therefore, we  create a small wrapper function that just
# takes the first ouput argument of ``ABmodel``.

model = lambda par: ABmodel(par)[0] # call ABmodel and take the first output (V)

# Fit the model to all traces simultaneously (global fit)
fit = dl.nlls(Vexp,model,par0,lower,upper,multistart=40)
# The use of the option 'multistart' helps the solver to find the
# global minimum and not to get stuck in a local minimum.

# Display fitted parameters and their uncertainties (standard deviations)
list(zip(fit.param, fit.paramUncert.std))

# %%
# Get the fitted model traces and distributions
Vfit, Pfit = ABmodel(fit.param)

# Get 95% confidence intervals of the fitted traces
Vfit_ci = [fit.modelUncert[i].ci(95) for i in range(N)]

# Propagate parameter uncertainty to the distribution models, accounting for non-negativity
Pfit_uq = [fit.paramUncert.propagate(lambda param: ABmodel(param)[1][i],lbm=np.zeros_like(r)) for i in range(N)]

# Get their 95% confidence intervals 
Pfit_ci = [Pfit_uq[i].ci(95) for i in range(N)]

# %% [markdown]
# Plot results
# ------------

plt.figure(figsize=(8,6))
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
