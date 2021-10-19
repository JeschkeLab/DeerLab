# %% 
""" 
Global fitting of multiple different DEER signals
-------------------------------------------------------------------------------

How to fit multiple signals from different DEER experiments to a model with a non-parametric
distribution and a homogeneous background, using Tikhonov regularization.
""" 

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl
# Use the seaborn style for nicer plots
from seaborn import set_theme
set_theme()
#%%

# Load the experimental 4-pulse and 5-pulse DEER datasets
t4p, V4p = np.load('../data/example_4pdeer_#2.npy')
t5p, V5p = np.load('../data/example_5pdeer_#2.npy')

# Since they have different scales, normalize both datasets
V4p = V4p/np.max(V4p)
V5p = V5p/np.max(V5p)

# Run fit
r = np.linspace(2,4.5,100)

# Construct the individual dipolar signal models
V4pmodel = dl.dipolarmodel(t4p,r,npathways=1)
V5pmodel = dl.dipolarmodel(t5p,r,npathways=2)
V5pmodel.reftime2.set(lb=3,ub=3.5,par0=3.2)

# Make the joint model with the distribution as a global parameters
globalmodel = dl.expand(V4pmodel,V5pmodel)  
globalmodel = dl.link(globalmodel, P = ['P_1','P_2'])

# Fit the model to the data (with fixed regularization parameter)
fit = dl.fit(globalmodel,[V4p,V5p],regparam=0.5)

# %%

plt.figure(figsize=[10,7])

# Extract fitted distance distribution
Pfit = fit.P
scale = np.trapz(Pfit,r)
Pci95 = fit.PUncert.ci(95)/scale
Pci50 = fit.PUncert.ci(50)/scale
Pfit =  Pfit/scale
for n,(t,V) in enumerate(zip([t4p,t5p],[V4p,V5p])):

    # Extract fitted dipolar signal
    Vfit = fit.model[n]
    Vci = fit.modelUncert[n].ci(95)

    plt.subplot(2,2,n+1)
    # Plot experimental data
    plt.plot(t,V,'.',color='grey',label='Data')
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
