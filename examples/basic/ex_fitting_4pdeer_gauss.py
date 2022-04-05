# %% [markdown]
""" 
Basic analysis of a 4-pulse DEER signal with a bimodal Gaussian model
-------------------------------------------------------------------------

Fit a simple 4-pulse DEER signal with a model with a bimodal Gaussian parametric 
distribution and a homogeneous background.
""" 

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl


# %%

# Load the experimental data
t,Vexp = np.load('../data/example_4pdeer_#2.npy')

# Distance vector
r = np.arange(1.5,6,0.01) # nm

# Construct the model
Pmodel= dl.dd_gauss2
Vmodel = dl.dipolarmodel(t,r,Pmodel)

# Fit the model to the data
results = dl.fit(Vmodel,Vexp,reg=False)

# Print results summary
print(results)

#%%

# Extract fitted dipolar signal
Vfit = results.model
Vci = results.modelUncert.ci(95)

# Extract fitted distance distribution
Pfit = results.evaluate(Pmodel,r)
scale = np.trapz(Pfit,r)
Puncert = results.propagate(Pmodel,r,lb=np.zeros_like(r))
Pfit = Pfit/scale
Pci95 = Puncert.ci(95)/scale
Pci50 = Puncert.ci(50)/scale

# Extract the unmodulated contribution
Bfcn = lambda mod,conc: scale*(1-mod)*dl.bg_hom3d(t,conc,mod)
Bfit = Bfcn(results.mod,results.conc)
Bci = results.propagate(Bfcn).ci(95)

plt.figure(figsize=[6,7])
violet = '#4550e6'
plt.subplot(211)
# Plot experimental and fitted data
plt.plot(t,Vexp,'.',color='grey',label='Data')
plt.plot(t,Vfit,linewidth=3,color=violet,label='Fit')
plt.fill_between(t,Vci[:,0],Vci[:,1],color=violet,alpha=0.3)
plt.plot(t,Bfit,'--',linewidth=3,color=violet,label='Unmodulated contribution')
plt.fill_between(t,Bci[:,0],Bci[:,1],color=violet,alpha=0.3)
plt.legend(frameon=False,loc='best')
plt.xlabel('Time $t$ (μs)')
plt.ylabel('$V(t)$ (arb.u.)')
# Plot the distance distribution
plt.subplot(212)
plt.plot(r,Pfit,color=violet,linewidth=3,label='Fit')
plt.fill_between(r,Pci95[:,0],Pci95[:,1],alpha=0.3,color=violet,label='95%-Conf. Inter.',linewidth=0)
plt.fill_between(r,Pci50[:,0],Pci50[:,1],alpha=0.5,color=violet,label='50%-Conf. Inter.',linewidth=0)
plt.legend(frameon=False,loc='best')
plt.autoscale(enable=True, axis='both', tight=True)
plt.xlabel('Distance $r$ (nm)')
plt.ylabel('$P(r)$ (nm$^{-1}$)')
plt.tight_layout()
plt.show()

# %%