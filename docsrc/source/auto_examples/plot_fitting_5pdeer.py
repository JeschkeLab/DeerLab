# %% [markdown]
"""
Fitting a 5-pulse DEER signal with a parameter-free distribution
==================================================================

This example shows how to fit a 5-pulse DEER signal with a parameter-
free distribution, a background, and all pathways parameters
""" 
# %%
import numpy as np
import matplotlib.pyplot as plt
from deerlab import *


# %%
# Generate data
np.random.seed(1)
t = np.linspace(-0.1,6.5,200)    # time axis, us
r = np.linspace(1.5,6,100)          # distance axis, ns
param0 = [3, 0.3, 0.2, 3.5, 0.3, 0.65, 3.8, 0.2, 0.15] # parameters for three-Gaussian model
P = dd_gauss3(r,param0)         # model distance distribution
B = lambda t,lam: bg_hom3d(t,300,lam) # background decay
exparam = [0.6, 0.3, 0.1, 3.2]     # parameters for 5pDEER experiment
pathinfo = ex_5pdeer(exparam)   # pathways information

K = dipolarkernel(t,r,pathinfo,B)
Vexp = K@P + whitegaussnoise(t,0.005)

# %% [markdown]
# Now, 5pDEER data contain 3 additional parameters compared to 4pDEER (due
# to the additional dipolar pathway present in the signal). However, the
# refocusing time of the second dipolar pathway is very easy to constrain
# and strongly helps stabilizing the fit. 
# 
# This pathway ususally refocuses at around ``t = max(t)/2``, and usually can
# be even estimated from simple visual inspection of the signal. 
# Thus, we can strongly constraint this parameters while leaving the
# pathway amplitudes pretty unconstrained.
# 

# %%
ex_lb   = [ 0,   0,   0,  max(t)/2-1] # lower bounds
ex_ub   = [100, 100, 100, max(t)/2+1] # upper bounds
ex_par0 = [0.5, 0.5, 0.5, max(t)/2  ] # start values

# %% [markdown]
# In this case we only want to set the bounds for the experiment
# parameters, so we can leave the rest empty:
# 

# %%
ub = [[],[],ex_ub]
lb = [[],[],ex_lb]
par0 = [[],[],ex_par0]

# %% [markdown]
# Run the fit with a 5-pulse DEER signal model

# %%
Vfit,Pfit,Bfit,parfit,paruq,moduq,stats = fitsignal(Vexp,t,r,'P',bg_hom3d,ex_5pdeer,par0,lb,ub,display=True)
plt.show()