# %% [markdown]
"""
Fitting a mixed distance-distribution model
===========================================

Basic manipulation of parametric models and creating mixed models 
for fitting distance distributions.
"""

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %% [markdown]
# Simulate the data
# -----------------
#
# Let's start by creating a simple dipolar evolution function (i.e. no background 
# and full modulation depth) corresponding to a simple 4-pulse DEER signal.

#Axis definition
t = np.linspace(-0.5,4,350)
r = np.linspace(2,6,200)

# Distribution parameters
rmean = 4.5
sigma = 0.2
chain = 4.3
pers = 10
amp = 0.35

# Generate distribution
P = dl.dd_gauss(r,[rmean, sigma])
P = amp*P + (1 - amp)*dl.dd_wormchain(r,[chain, pers])
# Normalize distribution
P = P/sum(P)/np.mean(np.diff(r))
# Generate dipolar evolution function
K = dl.dipolarkernel(t,r)
V = K @ P + dl.whitegaussnoise(t,0.02,seed=0)

# %%
# Generating a mixed parametric model
# -----------------------------------
#
# Let's say our intuiton (which, since we know the ground truth, is exact) on 
# the sample indicates that our distribution is a llinear combination of a Gaussian 
# distribution and a worm-like chain model. While DeerLab provides built-in 
# parametric models for both models, we require a combination of both. 
#
# For such cases we can use the ``mixmodels`` function to create a custom mixed 
# parametric model. It's syntax is rather simple, we just have to pass the desired 
# parametric models as lambda functions. 

#Mix the models into new one
gausswlc = dl.mixmodels(dl.dd_gauss,dl.dd_wormchain)

# %% [markdown]
# Our new model ``gausswlc`` will now describe our sought linear combination of 
# both parametric models. We can check the state of the model by retrieving its 
# information

#Get information on the mixed model
info = gausswlc()

# %% [markdown]
# We can see that the ``mixmodels`` function has introduced an ampitude parameters 
# as the first parameter of the model. This parameters weights the contribution 
# of each individual parametric model. We see also that this is followed by the 
# parameters of the Gaussian model and finally with the parameters of the worm-
# like chain model.
#
# Our model is ready, and since it was generated from built-in models we do 
# not need to specify any parameters initial values or boundary constraints. These 
# can, however, by re-defined if the built-in defaults are not appropiate (see 
# other examples). 
# 
# Since we are dealing with a distance-domain model we require a dipolar kernel 
# to transform our model into time-domain. Remember that our signal in this example 
# is a dipolar evolution function, therefore we do not require anything else than 
# a very basic dipolar kernel.

# Generate the dipolar evolution function kernel
K = dl.dipolarkernel(t,r)

# Fit the model to the data
Vmodel = lambda par: K @ gausswlc(r,par)
info = gausswlc()
par0 = info['Start'] # built-in start values
lb = info['Lower'] # built-in lower bounds
ub = info['Upper'] # built-in upper bounds
fit = dl.fitparamodel(V,Vmodel,par0,lb,ub,multistart=10)
fitpar = fit.param
# %% [markdown]
# From the fitted parameter set ``fitpar`` we can now generate our fitted distance 
# distribution and the corresponding time-domain fit.

# Calculate the fitted model
Pfit = gausswlc(r,fitpar)
Vfit = Vmodel(fitpar)

# %% [markdown]
# Since we know both the ground truth for the distance distribution and the 
# dipolar signal, let's see how our fit turned out.

# Plot results
plt.subplot(2,1,1)
plt.plot(t,V,'k.',t,Vfit,'r',linewidth=1.5)
plt.xlabel('t (µs)')
plt.ylabel('V')
plt.legend(['data','fit'])

plt.subplot(2,1,2)
plt.plot(r,P,'k',r,Pfit,'r',linewidth=1.5)
plt.xlabel('r (nm)')
plt.ylabel('P (nm⁻¹)')
plt.legend(['truth','fit'])
