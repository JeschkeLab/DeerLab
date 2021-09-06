# %% [markdown]
"""
Fitting Gaussians to a non-parametric distance distribution fit
============================================================================

This example shows how to fit multi-Gaussian model to a non-parametric distance
distribution calculated from Tikhonov regularization.  
""" 
# %%

# Import the required libraries
import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl
# Use the seaborn style for nicer plots
from seaborn import set_theme
set_theme()

# %%
    
# Load the experimental dataset
t,V = np.load('../data/example_data_#1.npy')

# Construct the dipolar signal model
r = np.linspace(1,7,100)
Vmodel = dl.dipolarmodel(t,r)

# Fit the model to the data
fit = dl.fit(Vmodel,V)
fit.plot(axis=t)
plt.ylabel('V(t)')
plt.xlabel('Time $t$ (μs)')
plt.show()

# From the fit results, extract the distribution and the covariance matrix
Pfit = fit.P
Pci95 = fit.PUncert.ci(95)

# Select a bimodal Gaussian model for the distance distribution
Pmodel = dl.dd_gauss2

# Fit the Gaussian model to the non-parametric distance distribution
fit = dl.fit(Pmodel,Pfit,r)

# Extract the fit results
PGauss = fit.model
PGauss_ci95 = fit.modelUncert.ci(95)

# Print the parameters nicely
print(f'Gaussian components with (95%-confidence intervals):')
print(f'       mean1 = {fit.mean1:2.2f} ({fit.mean1Uncert.ci(95)[0]:2.2f}-{fit.mean1Uncert.ci(95)[1]:2.2f}) nm')
print(f'       mean2 = {fit.mean2:2.2f} ({fit.mean2Uncert.ci(95)[0]:2.2f}-{fit.mean2Uncert.ci(95)[1]:2.2f}) nm')
print(f'      width1 = {fit.width1:2.2f} ({fit.width1Uncert.ci(95)[0]:2.2f}-{fit.width1Uncert.ci(95)[1]:2.2f}) nm')
print(f'      width2 = {fit.width2:2.2f} ({fit.width2Uncert.ci(95)[0]:2.2f}-{fit.width2Uncert.ci(95)[1]:2.2f}) nm')
print(f'  amplitude1 = {fit.amp1:2.2f} ({fit.amp1Uncert.ci(95)[0]:2.2f}-{fit.amp1Uncert.ci(95)[1]:2.2f})')
print(f'  amplitude2 = {fit.amp2:2.2f} ({fit.amp2Uncert.ci(95)[0]:2.2f}-{fit.amp2Uncert.ci(95)[1]:2.2f})')

# %%

# sphinx_gallery_thumbnail_number = 2

# Plot the fitted constraints model on top of the non-parametric case
plt.plot(r,Pfit,linewidth=1.5,label='Non-param. fit')
plt.fill_between(r,Pci95[:,0],Pci95[:,1],alpha=0.4,linewidth=0)
plt.plot(r,PGauss,linewidth=1.5,label='2-Gauss fit to non-param. fit',color='green')
plt.fill_between(r,PGauss_ci95[:,0],PGauss_ci95[:,1],alpha=0.2,linewidth=0,color='green')
# Formatting settings 
plt.xlabel('Distance (nm)')
plt.ylabel('P (nm$^{-1}$)')
plt.autoscale(enable=True, axis='both', tight=True)
plt.legend(loc='best',frameon=False)
plt.tight_layout()
plt.show()

# %%
