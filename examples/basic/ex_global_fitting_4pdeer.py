# %% [markdown]
""" 
Global fitting of multiple 4-pulse DEER signals, non-parametric distribution
-----------------------------------------------------------------------------

How to fit multiple 4-pulse DEER signals to a model with a non-parametric
distribution and a homogeneous background, using Tikhonov regularization.
""" 

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %% [markdown]
# Load and pre-process data
# ---------------------------
#
# Uncomment and use the following lines if you have experimental data:
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
# In this example we will use three simulated 4-pulse DEER signals instead:

#%%

# Define a function that generates synthetic data
def generatedata():
    r = np.linspace(2,5,200)                # distance axis, nm
    param = [3, 0.1, 0.2, 3.5, 0.1, 0.65]   # parameters for three-Gaussian model
    P = dl.dd_gauss2(r,param)               # model distance distribution

    t1 = np.linspace(-0.1,4,250)             # time axis 1, µs
    t2 = np.linspace(-0.2,3,200)             # time axis 2, µs
    t3 = np.linspace(-0.2,2,500)             # time axis 3, µs

    K1 = dl.dipolarkernel(t1,r,mod=0.3,bg=dl.bg_hom3d(t1,150,0.3))  # dipolar kernel 1
    K2 = dl.dipolarkernel(t2,r,mod=0.4,bg=dl.bg_hom3d(t2,80,0.4))   # dipolar kernel 2
    K3 = dl.dipolarkernel(t3,r,mod=0.4,bg=dl.bg_hom3d(t3,20,0.4))   # dipolar kernel 3

    V1 = K1@P + dl.whitegaussnoise(t1,0.01,seed=1) # simulated signal 1
    V2 = K2@P + dl.whitegaussnoise(t2,0.02,seed=2) # simulated signal 2
    V3 = K3@P + dl.whitegaussnoise(t3,0.01,seed=3) # simulated signal 3
    
    return [t1,t2,t3], [V1,V2,V3]

t, V = generatedata()

# %% [markdown]
# When doing global fitting, you must specify a list of the signals as well as a list of the corresponding time axes. 
# Then, a model type for the global distance distribution model and finally a list of models for the background and experiment.
# In this case we assume that both signals can be modelled as 3D-homogenous distributions of spins and that they are simple 4-pulse
# DEER experiments.

# %%

# Run fit
r = np.linspace(2,5,200)
fit = dl.fitmodel(V,t,r,'P',[dl.bg_hom3d,dl.bg_hom3d,dl.bg_hom3d],[dl.ex_4pdeer,dl.ex_4pdeer,dl.ex_4pdeer],verbose=True)
fit.plot();
