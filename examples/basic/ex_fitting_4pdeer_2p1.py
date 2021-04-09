# %% [markdown]
"""
Basic fitting of a 4-pulse DEER signal, with 2+1 contribution
====================================================================

This example shows how to fit a 4-pulse DEER signal with a non-parametric
distribution, including the "2+1" pathway contribution.
""" 
# %%
import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %% [markdown]
# Load and pre-process data
# ---------------------------
#
# Uncomment and use the following lines if you have experimental data:
#
# t,Vexp = dl.deerload('my\path\4pdeer_data.DTA')
# Vexp = dl.correctphase(Vexp)
# t = dl.correctzerotime(Vexp,t)
#
# In this example we will use simulated data instead:

# %%
#Generate data
#--------------
#
# In this example we will use simulated data instead:

#%%

# Define a function that generates synthetic data
def generatedata():
    t = np.linspace(-0.2,4,200)                             # time axis, µs
    r = np.linspace(2,5,200)                                # distance axis, nm
    param0 = [3, 0.1, 0.2, 3.5, 0.1, 0.65, 3.8, 0.05, 0.15] # parameters for three-Gaussian model
    P = dl.dd_gauss3(r,param0)                              # model distance distribution
    B = lambda t,lam: dl.bg_hom3d(t,300,lam)                # background decay
    exparam = [0.6, 0.3, 0.1, 4.1]                          # parameters for the experiment model
    pathways = dl.ex_ovl4pdeer(exparam)                     # pathways information
    K = dl.dipolarkernel(t,r,pathways=pathways,bg=B)
    Vexp = K@P + dl.whitegaussnoise(t,0.01,seed=1)
    return t, Vexp
    
t, Vexp = generatedata()

# %% [markdown]
# Now, if we take the "2+1" contribution into account, the 4pDEER model
# contains 3 additional parameters compared to a simplified 4pDEER model (due
# to the additional dipolar pathway present in the signal). However, the
# refocusing time of the second dipolar pathway is very easy to constrain
# and strongly helps stabilizing the fit. 
# 
# This pathway can even estimated visually from the signal or estimated from the 
# pulse sequence timings. Thus, we can strongly constraint this parameters while leaving the
# pathway amplitudes pretty much unconstrained.
# 
# %%

# Define initial values and bounds for model parameters
ex_lb   = [ 0,   0,   0,  max(t)-1] # lower bounds
ex_ub   = [1, 1, 1, max(t)+1] # upper bounds
ex_par0 = [0.3, 0.3, 0.3, max(t)  ] # start values

# %% [markdown]
# Run the fit with a 4-pulse DEER model that includes the "2+1" pathway contribution

# %%
r = np.linspace(2,5,200)
fit = dl.fitmodel(Vexp,t,r,'P',dl.bg_hom3d,dl.ex_ovl4pdeer,
                  ex_par0=ex_par0,ex_lb=ex_lb,ex_ub=ex_ub,verbose=True)
fit.plot();

# %%
