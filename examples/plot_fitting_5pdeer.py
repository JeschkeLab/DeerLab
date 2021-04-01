# %% [markdown]
"""
Fitting a 5-pulse DEER signal with a parameter-free distribution
==================================================================

This example shows how to fit a 5-pulse DEER signal with a non-parametric
distribution, a background, and all pathways parameters.
""" 
# %%
import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %%
# Generate data
t = np.linspace(-0.1,6.5,200)      # time axis, µs
r = np.linspace(1.5,6,100)         # distance axis, nm
param0 = [3, 0.1, 0.2, 3.5, 0.1, 0.65, 3.8, 0.05, 0.15] # parameters for three-Gaussian model
P = dl.dd_gauss3(r,param0)         # model distance distribution
B = lambda t,lam: dl.bg_hom3d(t,300,lam) # background decay
exparam = [0.6, 0.3, 0.1, 3.2]     # parameters for 5pDEER experiment
pathways = dl.ex_5pdeer(exparam)   # pathways information

K = dl.dipolarkernel(t,r,pathways=pathways,bg=B)
Vexp = K@P + dl.whitegaussnoise(t,0.005,seed=1)

# %% [markdown]
# Now, 5pDEER data contain 3 additional parameters compared to 4pDEER (due
# to the additional dipolar pathway present in the signal). However, the
# refocusing time of the second dipolar pathway is very easy to constrain
# and strongly helps stabilizing the fit. 
# 
# This pathway ususally refocuses at around ``t = max(t)/2``, and usually can
# be even estimated from simple visual inspection of the signal. 
# Thus, we can strongly constraint this parameters while leaving the
# pathway amplitudes pretty much unconstrained.
# 

# %% [markdown]
# Define initial values and bounds for model parameters
# %%
ex_lb   = [ 0,   0,   0,  max(t)/2-1] # lower bounds
ex_ub   = [100, 100, 100, max(t)/2+1] # upper bounds
ex_par0 = [0.5, 0.5, 0.5, max(t)/2  ] # start values

# %% [markdown]
# Run the fit with a 5-pulse DEER signal model

# %%
fit = dl.fitmodel(Vexp,t,r,'P',dl.bg_hom3d,dl.ex_5pdeer,ex_par0=ex_par0,ex_lb=ex_lb,ex_ub=ex_ub)
fit.plot()
