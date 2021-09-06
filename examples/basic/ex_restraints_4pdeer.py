# %%
""" 
Distance restraints from 4-pulse DEER data, non-parametric distribution
-----------------------------------------------------------------------

How to fit a simple 4-pulse DEER signal and derive distance restraints from 
the fitted non-parametric distance distribution.

Once we have a fit of the distance distribution we can obtain distance
restraints in the form of different statistical descriptors such as the mean distance
and standard deviation of distances. While we could calculate this manually, DeerLab 
provides a convenient function ``diststats`` which will automatically compute
these for you and even propagate the uncertainty in the distributions to those 
values to get confidence intervals on the restraints. 
""" 

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl
# Use the seaborn style for nicer plots
from seaborn import set_theme
set_theme()

#%% 

# Load experimental data
t,Vexp = np.load('../data/example_4pdeer_#1.npy')

# Distance vector
r = np.linspace(2,5,100) # nm

# Construct dipolar model
Vmodel = dl.dipolarmodel(t,r)

# Fit the model to the data
fit = dl.fit(Vmodel,Vexp)

# Get printed summary of all statistical descriptors available with confidence intervals
estimators,uq = dl.diststats(r,fit.P,fit.PUncert,verbose=True)

# Get the mean distance
rmean = estimators['mean']
rmean_ci = uq['mean'].ci(95)
# Get the standard deviation of distances
r_std = estimators['std']
r_std_ci = uq['std'].ci(95)

# Print out the results
print(f'Mean distance: {rmean:.3f} ({rmean_ci[0]:.3f}-{rmean_ci[1]:.3f}) nm')
print(f'Standard deviation: {r_std:.3f} ({r_std_ci[0]:.3f}-{r_std_ci[1]:.3f}) nm')

#%%

# For display, you can plot the mean distance with its confidence intervals without further calculations.

# Plot distribution and confidence bands
Pci95 = fit.PUncert.ci(95)/np.trapz(fit.P,r)
Pci50 = fit.PUncert.ci(50)/np.trapz(fit.P,r)
plt.plot(r,fit.P/np.trapz(fit.P,r),linewidth=2,label='Distance distribution fit')
plt.fill_between(r,Pci95[:,0],Pci95[:,1],color='tab:blue',alpha=0.3)
plt.fill_between(r,Pci50[:,0],Pci50[:,1],color='tab:blue',alpha=0.4)

# Plot mean distance and confidence interval
plt.vlines(rmean,0,max(Pci95[:,1]),color='tab:red',linestyles='dotted',linewidth=3,label='Mean distance')
plt.fill_between(rmean_ci,0,max(Pci95[:,1]),color='tab:red',alpha=0.3,linewidth=0)

plt.legend(frameon=False,loc='best')

plt.ylim([0,max(Pci95[:,1])])
plt.xlabel('Distance $r$ (nm)')
plt.ylabel('P(r) (nm$^{-1}$)')
plt.autoscale(enable=True, axis='both', tight=True)
plt.tight_layout() 
plt.show()
# %%
