# %%
"""
Comparing confidence intervals for regularization results
=========================================================

A simple example of uncertainty estimation of results obtained from dipolar signals.
This example runs the analysis of a 4-pulse DEER signal and compares the uncertainty of the
distance distribution between the covariance (curvature matrix) and bootstrap methods.

By plotting the results, one can see that the bootstrapped confidence intervals 
are narrower in comparison to the ones obtained via the curvature
matrices. This is because bootstrapping takes the nonnegativity constraint of P(r) into
account, whereas the curvature matrix CIs do not. 
""" 

import numpy as np 
import matplotlib.pyplot as plt
import deerlab as dl
# Use the seaborn style for nicer plots
from seaborn import set_theme
set_theme()

# %% 

# Load the experimental data
t,Vexp = np.load('../data/example_data_#1.npy')

# Distance vector
r = np.linspace(1,7,80)

# Construct the 4-pulse DEER dipolar model
Vmodel = dl.dipolarmodel(t,r)

# Fit the model to the data using covariane-based uncertainty
fit_cm = dl.fit(Vmodel,Vexp)

# Fit the model to the data using bootstrapped uncertainty
fit_bs = dl.fit(Vmodel,Vexp,bootstrap=10)

# Compute the covariance-based uncertainty bands of the distance distribution
Pci50_cm = fit_cm.PUncert.ci(50)
Pci95_cm = fit_cm.PUncert.ci(95)

# Compute the bootstrapped uncertainty bands of the distance distribution
Pci50_bs = fit_bs.PUncert.ci(50)
Pci95_bs = fit_bs.PUncert.ci(95)

#%%

# Plot the results
fig, ax = plt.subplots(1,2,sharey=True)
ax[0].plot(r,fit_cm.P,'tab:red',linewidth=1)
ax[0].fill_between(r,Pci50_cm[:,0],Pci50_cm[:,1],color='tab:red',linestyle='None',alpha=0.45)
ax[0].fill_between(r,Pci95_cm[:,0],Pci95_cm[:,1],color='tab:red',linestyle='None',alpha=0.25)

ax[1].plot(r,fit_bs.P,'tab:blue',linewidth=1)
ax[1].fill_between(r,Pci50_bs[:,0],Pci50_bs[:,1],color='tab:blue',linestyle='None',alpha=0.45)
ax[1].fill_between(r,Pci95_bs[:,0],Pci95_bs[:,1],color='tab:blue',linestyle='None',alpha=0.25)

ax[0].set_xlabel('Distance $r$ (nm)')
ax[0].set_ylabel('$P(r)$ (nm$^{-1}$)')
ax[0].set_title('Curvature Matrix CI')
ax[0].legend(['Median','50%-CI','95%-CI'],frameon=False,loc='best')

ax[1].set_xlabel('Distance $r$ (nm)')
ax[1].set_title('Bootstrapped CI')
ax[1].legend(['Median','50%-CI','95%-CI'],frameon=False,loc='best')
plt.tight_layout()
plt.show()
# %%
