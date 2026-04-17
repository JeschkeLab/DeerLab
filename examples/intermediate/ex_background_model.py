# %% [markdown]
""" 
Generating a background model for a 4-pulse DEER experiment
-------------------------------------------------------------------------

Since DeerLab v1.2, the background model can be generated for a given experiment using the :func:`deerlab.dipolarbackgroundmodel` function, as is normally automatically calculated in the fit function as `FitResult.bg` and  `FitResult.bgUncert`.

However, it can still be useful to generate the background model separately, and so this example shows how to do that for a 4-pulse DEER experiment.

First we will fit a simple 4-pulse DEER dataset with multiple-pathways. 
""" 


import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl
# %%
# File location
path = '../data/'
file = 'example_4pdeer_1.DTA'

# Experimental parameters
tau1 = 0.3      # First inter-pulse delay, μs
tau2 = 4.0      # Second inter-pulse delay, μs
tmin = 0.1      # Start time, μs

# Load the experimental data
t,Vexp = dl.deerload(path + file)

# Pre-processing
Vexp = dl.correctphase(Vexp) # Phase correction
Vexp = Vexp/np.max(Vexp)     # Rescaling (aesthetic)
t = t - t[0]                  # Account for zerotime
t = t + tmin    
# Distance vector
r = np.arange(2.5,5,0.01) # nm

# Construct the model
experimentInfo = dl.ex_4pdeer(tau1,tau2, pathways=[1,2,3])
Vmodel = dl.dipolarmodel(t,r, experiment = experimentInfo, Bmodel=dl.bg_hom3d)

# Fit the model to the data
results = dl.fit(Vmodel,Vexp)

# Print results summary
print(results)


# %%
"""
Now using DeerLab > v1.2, the background model can be directly plotted from the fit results.
"""
plt.figure(figsize=[6,7])
violet = '#4550e6'
plt.plot(t,Vexp,'.',color='grey',label='Data')
plt.plot(t,results.model,linewidth=3,color=violet,label='Fit')
plt.plot(t,results.bg,'--',linewidth=3,color=violet,label='Unmodulated contribution')


# %%
"""
Alternatively we could generate the background model separately using the `dl.dipolarbackgroundmodel` function, which takes the experimental information as input.
"""

bg_model = dl.dipolarbackgroundmodel(experimentInfo, basis = dl.bg_hom3d)
bg = results.evaluate(bg_model,t)
bgUncert = results.propagate(bg_model,t)

"""
In DeerLab <1.2, the multi-pathway background model needs to be generated manually from a function.
"""

bf_func = lambda lam1,lam2,lam3,reftime1,reftime2,reftime3,conc: results.P_scale*(1-lam1-lam2-lam3)*dl.bg_hom3d(t-reftime1,conc,lam1)*dl.bg_hom3d(t-reftime2,conc,lam2)*dl.bg_hom3d(t-reftime3,conc,lam3)
bf = results.evaluate(bf_func, t)
bfUncert = results.propagate(bf_func, t)
