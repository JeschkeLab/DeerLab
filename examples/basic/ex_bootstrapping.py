# %% [markdown]
""" 
Bootstrapped confidence intervals in routine analysis
-------------------------------------------------------------------

How to obtain bootstrapped confidence intervals for simple routine operations with ``fitmodel``.
""" 

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %% [markdown]
# Experimental data must be loaded and pre-processed::
#
#        t,Vexp = dl.deerload('my\path\4pdeer_data.DTA')
#        Vexp = dl.correctphase(Vexp)
#        t = dl.correctzerotime(Vexp,t)
#


# %% [markdown]
# In this example we will use simulated data instead:

#%% 

t = np.linspace(-0.1,4,250)        # time axis, µs
r = np.linspace(2,5,200)           # distance axis, nm
param = [3, 0.1, 0.2, 3.5, 0.1, 0.65, 3.8, 0.05, 0.15] # parameters for three-Gaussian model
P = dl.dd_gauss3(r,param)          # model distance distribution
lam = 0.5                          # modulation depth
B = dl.bg_hom3d(t,300,lam)         # background decay
K = dl.dipolarkernel(t,r,mod=lam,bg=B)    # kernel matrix
Vexp = K@P + dl.whitegaussnoise(t,0.01,seed=0)

# %% [markdown]
# Unless specified otherwise, the function ``fitmodel`` will return asymptotic confidence intervals based on the covariance matrix 
# of the objective function used to fit the data. These are quick to calculate and therefore very comfortable for quick estimates of
# the uncertainty during routine analysis or testing. 
#
# However, for publication-level analysis, these confidence intervals might be inaccurate. It is strongly recommended to employ bootstrapped 
# confidence intervals to get accurate estimates of the uncertainty. 
# Conviniently, ``fitmodel`` integrates bootstrapping to make it accessible by only adjusting the option ``'uq'`` to ``'bootstrap'``. 
# The option also allows to specify the number of samples to analyze to estimate the uncertainty. The larger this number, the more accurate 
# the confidence intervals but the longer the analysis will be. The standard for publication is typically 1000 samples. 
#
# If the option ``'verbose'`` is enabled, a progress bar will be shown indicating the remaining time for the bootstrapping to finish.
#
# In this example, for the sake of time, we will just use 50 bootstrap samples.  

#%%
BootSamples = 50 # For publication-grade analysis, the standard is 1000
fit = dl.fitmodel(Vexp,t,r,'P',dl.bg_hom3d,dl.ex_4pdeer,uq=['bootstrap',BootSamples])
fit.plot();

# %% [markdown]
# Note that in the automatic plot, instead of the fit, the median of the bootstrapped distribution is shown. 