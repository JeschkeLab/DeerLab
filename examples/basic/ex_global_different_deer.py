# %% [markdown]
""" 
Global fitting of multiple different DEER signals, non-parametric distribution
-------------------------------------------------------------------------------

How to fit multiple signals from different DEER experiments to a model with a non-parametric
distribution and a homogeneous background, using Tikhonov regularization.
""" 

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %% [markdown]
# Load and pre-process data
# ---------------------------
#
# All experimental data must be loaded and pre-processed::
#
# datasets = ('file1.DTA','file2.DTA','file3.DTA')
# data = [dl.deerload(ds) for ds in datasets]
# t = [_[0] for _ in data]
# V = [_[1] for _ in data]
#

# %% [markdown]
# Simulate data
# --------------
#
# In this example we will use two simulated signals (one 4-pulse DEER signal and one 5-pulse DEER signal) instead:

#%%

# Define a function that generates synthetic data
def generatedata():
    r = np.linspace(2,5,150)                # distance axis, nm
    param = [3, 0.1, 0.2, 3.5, 0.1, 0.65]   # parameters for three-Gaussian model
    P = dl.dd_gauss2(r,param)               # model distance distribution

    t1 = np.linspace(-0.2,3,200)             # time axis 1, µs
    t2 = np.linspace(-0.1,6,400)             # time axis 2, µs

    path4p = dl.ex_4pdeer(0.3)                      # 4-pulse DEER pathways
    path5p = dl.ex_5pdeer([0.6, 0.3, 0.1, 3.2])     # 5-pulse DEER pathways

    K1 = dl.dipolarkernel(t1,r,pathways=path4p,bg=lambda t,lam: dl.bg_hom3d(t1,100,lam))  # dipolar kernel 1
    K2 = dl.dipolarkernel(t2,r,pathways=path5p,bg=lambda t,lam: dl.bg_hom3d(t2,20,lam))   # dipolar kernel 2

    V1 = K1@P + dl.whitegaussnoise(t1,0.005,seed=1) # simulated signal 1
    V2 = K2@P + dl.whitegaussnoise(t2,0.01,seed=2) # simulated signal 2
    
    return [t1,t2], [V1,V2]

t, V = generatedata()

# %% [markdown]
# When doing global fitting, you must specify a list of the signals as well as a list of the corresponding time axes. 
# Then, a model type for the global distance distribution model and finally a list of models for the background and experiment.
# In this case we assume that both signals can be modelled as 3D-homogenous distributions of spins. 
# Since we know that the 

#%%

# Change the start values of the 5-pulse DEER model to reflect the experimental data
dl.ex_5pdeer.start = [0.3, 0.3, 0.3, 3.2]

# Run fit
r = np.linspace(2,5,150)
fit = dl.fitmodel(V,t,r,'P',[dl.bg_hom3d,dl.bg_hom3d],[dl.ex_4pdeer,dl.ex_5pdeer],verbose=True,weights=[2,1])
fit.plot();

# %%
