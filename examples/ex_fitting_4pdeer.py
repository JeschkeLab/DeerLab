# %% [markdown]
""" 
Basic fitting of a 4-pulse DEER signal, non-parametric distribution
-------------------------------------------------------------------

Fit a simple 4-pulse DEER signal with a model with a non-parametric
distribution and a homogeneous background, using Tikhonov regularization.
""" 

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %% [markdown]
#Generate data
#--------------

t = np.linspace(-0.1,4,250)        # time axis, µs
r = np.linspace(1.5,6,len(t))      # distance axis, ns
param = [3, 0.1, 0.2, 3.5, 0.1, 0.65, 3.8, 0.05, 0.15] # parameters for three-Gaussian model
P = dl.dd_gauss3(r,param)          # model distance distribution
lam = 0.5                          # modulation depth
B = dl.bg_hom3d(t,300,lam)         # background decay
K = dl.dipolarkernel(t,r,mod=lam,bg=B)    # kernel matrix
Vexp = K@P + dl.whitegaussnoise(t,0.01,seed=0)

# %% [markdown]
# Run fit
#---------
fit = dl.fitmodel(Vexp,t,r,'P',dl.bg_hom3d,dl.ex_4pdeer)
fit.plot()

# %%
