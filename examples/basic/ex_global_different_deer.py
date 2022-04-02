# %% 
""" 
Global fitting of different multi-pulse DEER signals
-------------------------------------------------------------------------------

How to fit multiple signals from different multi-pulse DEER experiments to a single
global distance distribution.
""" 

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

#%%

# Load the experimental 4-pulse and 5-pulse DEER datasets
t4p, V4p = np.load('D:\lufa\projects\DeerLab\DeerLab\examples/data/example_4pdeer_#1.npy')
t5p, V5p = np.load('D:\lufa\projects\DeerLab\DeerLab\examples/data/example_5pdeer_#2.npy')

# Since they have different scales, normalize both datasets
V4p = V4p/np.max(V4p)
V5p = V5p/np.max(V5p)

# Run fit
r = np.arange(2,5,0.05)

# Construct the individual dipolar signal models:
# 4-pulse DEER model
ex4pdeer = dl.ex_4pdeer(tau1=0.1,tau2=3.9,pathways=[1,2,3])
V4pmodel = dl.dipolarmodel(t4p,r,experiment=ex4pdeer)
# 5-pulse DEER model
ex5pdeer = dl.ex_rev5pdeer(tau1=2.7,tau2=3.3,tau3=0.1,pathways=[1,2])
V5pmodel = dl.dipolarmodel(t5p,r,experiment=ex5pdeer)

# Make the joint model with the distribution as a global parameters
globalmodel = dl.merge(V4pmodel,V5pmodel)  
globalmodel = dl.link(globalmodel, P = ['P_1','P_2'])

# Fit the model to the data (with fixed regularization parameter)
compactness = dl.dipolarpenalty(Pmodel=None,r=r,type='compactness')
fit = dl.fit(globalmodel,[V4p,V5p],regparam=0.1, penalties=compactness)

# %%

plt.figure(figsize=[10,7])
violet = '#4550e6'

# Extract fitted distance distribution
Pfit = fit.P
Pci95 = fit.PUncert.ci(95)
Pci50 = fit.PUncert.ci(50)

for n,(t,V) in enumerate(zip([t4p,t5p],[V4p,V5p])):

    # Extract fitted dipolar signal
    Vfit = fit.model[n]
    Vci = fit.modelUncert[n].ci(95)

    plt.subplot(2,2,n+1)
    # Plot experimental data
    plt.plot(t,V,'.',color='grey',label='Data')
    # Plot the fitted signal 
    plt.plot(t,Vfit,linewidth=3,color=violet,label='Fit')
    plt.fill_between(t,Vci[:,0],Vci[:,1],color=violet,alpha=0.3)
    plt.legend(frameon=False,loc='best')
    plt.xlabel('Time $t$ (μs)')
    plt.ylabel('$V(t)$ (arb.u.)')

# Plot the distance distribution
plt.subplot(212)
plt.plot(r,Pfit,linewidth=3,color=violet,label='Fit')
plt.fill_between(r,Pci95[:,0],Pci95[:,1],alpha=0.3,color=violet,label='95%-Conf. Inter.',linewidth=0)
plt.fill_between(r,Pci50[:,0],Pci50[:,1],alpha=0.5,color=violet,label='50%-Conf. Inter.',linewidth=0)
plt.legend(frameon=False,loc='best')
plt.autoscale(enable=True, axis='both', tight=True)
plt.xlabel('Distance $r$ (nm)')
plt.ylabel('$P(r)$ (nm$^{-1}$)')
plt.tight_layout()
plt.show()

# %%
