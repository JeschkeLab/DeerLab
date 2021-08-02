# %% [markdown] 
"""
Comparing confidence intervals for regularization results
=========================================================

A simpe example of uncertainty estimation for Tikhonov regularization 
results. The example will cover the use of confidence intervals
obtained from curvature matrices and boostrap analysis.
""" 

import numpy as np 
import matplotlib.pyplot as plt
import deerlab as dl

# %% [markdown] 
# Simulate the data
# ------------------
#
# Let's start by generating some data.

t = np.linspace(-0.4,3.5,160)           # time axis, µs
r = np.linspace(2,5,160)                # distance axis, nm
P = dl.dd_gauss2(r,[3, 0.1, 0.6, 3.5, 0.2, 0.4]) # model distribution
lam = 0.32                              # modulatio depth
B = dl.bg_strexp(t,[0.04,1])        # background decay
K = dl.dipolarkernel(t,r,mod=lam,bg=B)         # dipolar kernel matrix
V = K@P + dl.whitegaussnoise(t,0.01)    # signal with added noise

# %% [markdown] 
# For the sake of simplicity, in this examples we will assume that we know the 
# background exactly. Our first step is to generate the proper dipolar kernel.
#
# Covariance-based confidence intervals
# -------------------------------------
#
# Fit a Tikhonov model to the data, using AIC to select the regularization parameter
fit = dl.rlls(V,K,regparam='aic')

Pfit = fit.param        # fitted distribution
Vfit = fit.model        # fitted DEER trace
normfactor = np.trapz(Pfit,r)
Pfit = Pfit/normfactor

# curvature matrix confidence intervals for distribution
Pci95_cm = fit.paramUncert.ci(95)/normfactor
Pci50_cm = fit.paramUncert.ci(50)/normfactor

# %% [markdown]
# Bootstrapped confidence intervals
# ---------------------------------
#
# Now we calculate the bootstrap confidence intervals. For this, we
# need to define a function that takes a signal as input and returns
# the outputs of interest (``Pfit`` in our example).

def mybootfcn(V):
    fit = dl.rlls(V,K,regparam='aic')
    Pfit = fit.param
    Pfit = Pfit/np.trapz(Pfit,r)
    return Pfit


# Launch bootstrapping
Nsamples = 100
booci = dl.bootan(mybootfcn,V,Vfit,Nsamples)
Pci95_bs = booci.ci(95)
Pci50_bs = booci.ci(50)

# %% [markdown]
# By plotting the results, one can see that the bootstrapped confidence intervals 
# are narrower in comparison to the ones obtained via the curvature
# matrices. This is because bootstrapping takes the nonnegativity constraint of P into
# account, whereas the curvature matrix CIs do not. 

fig, ax = plt.subplots(2,1,sharey=True)
ax[0].plot(r,Pfit,'tab:red',linewidth=1)
ax[0].fill_between(r,Pci50_cm[:,0],Pci50_cm[:,1],color='tab:red',linestyle='None',alpha=0.45)
ax[0].fill_between(r,Pci95_cm[:,0],Pci95_cm[:,1],color='tab:red',linestyle='None',alpha=0.25)

ax[1].plot(r,Pfit,'tab:blue',linewidth=1)
ax[1].fill_between(r,Pci50_bs[:,0],Pci50_bs[:,1],color='tab:blue',linestyle='None',alpha=0.45)
ax[1].fill_between(r,Pci95_bs[:,0],Pci95_bs[:,1],color='tab:blue',linestyle='None',alpha=0.25)

ax[0].grid(alpha=0.5)
ax[0].set_xlabel('r (nm)')
ax[0].set_ylabel('P (nm⁻¹)')
ax[0].set_title('Curvature Matrix CI')
ax[0].legend(['Truth','Fit','50%-CI','95%-CI'])

ax[1].grid(alpha=0.5)
ax[1].set_xlabel('r (nm)')
ax[1].set_ylabel('P (nm⁻¹)')
ax[1].set_title('Bootstrapped CI')
ax[1].legend(['Truth','Fit','50%-CI','95%-CI'])
plt.tight_layout()
plt.show()
# %%
