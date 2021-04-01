# %% [markdown]
"""
Selecting an optimal parametric model for fitting a dipolar signal
==================================================================

How to optimally select a parametric model for a given dipolar signal.
"""

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %% [markdown]
# Data Generation
# ----------------
#
# Let's start by constructing a simple dipolar signal with some noise arising 
# from a bimodal Gaussian distance distribution.

# Prepare the signal components
t = np.linspace(-0.3,3.5,300)                # time axis, µs
r = np.linspace(2,6,200)                     # distance axis, nm
P = dl.dd_gauss2(r,[3.8, 0.4, 0.7, 4.5, 0.2, 0.7])   # distance distribution
K = dl.dipolarkernel(t,r)                    # dipolar kernel matrix
V = K@P + dl.whitegaussnoise(t,0.02)         # DEER signal, with added noise

# %% [markdown]
# Selecting an optimal model
# --------------------------
#
# Even though we know the ground truth, in this example we will cosider the 
# following set of potential parametric models: 
#
# * Unimodal Rician distribution
# * Bimodal Rician distribution
# * Trimodal Rician distribution
# * Unimodal Gaussian distribution
# * Bimodal Gaussian distribution
# * Trimodal Gaussian distribution
# * Mixed bimodal Gaussian/Rician distribution
#
# The first six models have built-in parametric models which we can use directly. 
# The last model we can construct from built-in models using the ``mixmodels`` function.

# Prepare the mixed model
dd_rice_gauss = dl.mixmodels(dl.dd_rice,dl.dd_gauss)
 
# Prepare list of candidate parametric models
models = [dl.dd_rice,dl.dd_rice2,dl.dd_rice3,dl.dd_gauss,dl.dd_gauss2,dl.dd_gauss3,dd_rice_gauss]

# %% [markdown]
# In order to make an appropiate choice, we need some liklihood estimator. All fit functions is DeerLab returns a stats 
# dictionary which contains (amongst other estimators) likelihood estimators such as the Akaike information criterion (AIC).
# The model with the lowers AIC value can be considered to most likely to be the optimal model.
#
# To do this, we just have to evaluate the parametric models with ``fitparamodel`` while looping over all the distribution models
# we listed above, and collecting the AIC-values for each model.
 
aic = []
for model in models:
    info = model()
    # Prepare the signal model with the new distance model
    Vmodel = lambda par: K@model(r,par)
    # Fit the signal
    fit = dl.fitparamodel(V,Vmodel,par0=info['Start'],lb=info['Lower'],ub=info['Upper'])
    parfit = fit.param
    stats= fit.stats
    # Add current AIC value to the list
    aic.append(stats['aic'])

# %% [markdown]
# Since the absolute AIC values have no meaning, it is standard practice to look at the relative 
# changes in AIC values between the evaluated models.

daic = aic - min(aic)

# %% [markdown]
# Akaike Weights
#-----------------------------------------------------------------------------
# It is often more useful to look at these results from the perspective of
# Akaike weights, i.e. the probabilities of a model being the most optimal.

weights = 100*np.exp(-(daic/2))/sum(np.exp(-daic/2))

# %% [markdown]
# Plot results
# ------------

plt.figure(figsize=(9,8))

plt.subplot(2,2,1)
plt.plot(t,V,'k.')
plt.grid(alpha=0.2)
plt.xlabel('t (µs)')
plt.legend(['data'])

plt.subplot(2,2,2)
plt.plot(r,P,'k',linewidth=1.5)
plt.xlabel('r (nm)')
plt.ylabel('P (nm⁻¹)')
plt.legend(['Ground truth'])
plt.grid(alpha=0.2)

modelnames = [model.__name__ for model in models]

plt.subplot(2,2,3)
plt.bar(modelnames,daic,color='b',alpha=0.5)
plt.ylabel('$\Delta$AIC')
plt.grid(alpha=0.2)
plt.xticks(rotation=45)

# Plot the results
plt.subplot(2,2,4)
plt.bar(modelnames,weights,color='b',alpha=0.5)
plt.ylabel('Akaike Weights (%)')
plt.xticks(rotation=45)
plt.grid(alpha=0.2)

# %% [markdown]
# Typically there is not a single optimal model unless the noise level is very
# low. Usually several models have similar probabilities and should therefore be presented together. 


# %%
