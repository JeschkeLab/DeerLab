# %% [markdown]
""" 
Basic analysis of a 4-pulse DEER signal with multiple dipolar pathays
-------------------------------------------------------------------------

Fit a simple 4-pulse DEER signal with a model with a non-parametric
distribution and a homogeneous background, using Tikhonov regularization.
""" 

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl


# %%

# Load the experimental data
t,Vexp = np.load('../data/example_data_#3.npy')

# Distance vector
r = np.linspace(2,5,100) # nm

# Construct the model
Vmodel = dl.dipolarmodel(t,r,npathways=3)

# Adjust the boundaries for the refocusing times
Vmodel.reftime1.set(par0=0.5, lb=0.0, ub=1.0) # Main pathway contribution
Vmodel.reftime2.set(par0=0.0, lb=0.0, ub=0.2) # Pathway refocusing at the start of the signal
Vmodel.reftime3.set(par0=4.5, lb=4.0, ub=5.0) # Pathway refocusing at the end of the signal

# Fit the model to the data
fit = dl.fit(Vmodel,Vexp)

#%%

# Extract fitted dipolar signal
Vfit = fit.model
Vci = fit.modelUncert.ci(95)

# Extract fitted distance distribution
Pfit = fit.P
scale = np.trapz(Pfit,r)
Pci95 = fit.PUncert.ci(95)/scale
Pci50 = fit.PUncert.ci(50)/scale
Pfit =  Pfit/scale


plt.figure(figsize=[6,7])
plt.subplot(211)
# Plot experimental data
plt.plot(t,Vexp,'.',color='grey',label='Data')
# Plot the fitted signal 
plt.plot(t,Vfit,linewidth=3,label='Fit')
plt.fill_between(t,Vci[:,0],Vci[:,1],alpha=0.3)
plt.legend(frameon=False,loc='best')
plt.xlabel('Time $t$ (μs)')
plt.ylabel('$V(t)$ (arb.u.)')
# Plot the distance distribution
plt.subplot(212)
plt.plot(r,Pfit,linewidth=3,label='Fit')
plt.fill_between(r,Pci95[:,0],Pci95[:,1],alpha=0.3,color='tab:blue',label='95%-Conf. Inter.',linewidth=0)
plt.fill_between(r,Pci50[:,0],Pci50[:,1],alpha=0.5,color='tab:blue',label='50%-Conf. Inter.',linewidth=0)
plt.legend(frameon=False,loc='best')
plt.autoscale(enable=True, axis='both', tight=True)
plt.xlabel('Distance $r$ (nm)')
plt.ylabel('$P(r)$ (nm$^{-1}$)')
plt.tight_layout()
plt.show()

# %%
