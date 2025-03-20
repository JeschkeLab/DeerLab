# %% [markdown]
""" 
Basic analysis of a RIDME signal
-------------------------------------------------------------------------

Fit a simple RIDME signal with a model with a non-parametric
distribution and a stretched exponetial background, using Tikhonov regularization.
""" 

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl


# %%

# File location
path = '../data/'
file = 'example_ridme_1.DTA'

# Experimental parameters
tau1 = 0.4      # First inter-pulse delay, μs
tau2 = 4.2      # Second inter-pulse delay, μs
tmin = 0.28  # Start time, μs

# Load the experimental data
t,Vexp = dl.deerload(path + file)

# Pre-processing
Vexp = dl.correctphase(Vexp) # Phase correction
Vexp = Vexp/np.max(Vexp)     # Rescaling (aesthetic)
t = t - t[0]             # Account for zerotime
t = t + tmin             

# Distance vector
r = np.linspace(1.5,6,50) # nm

# Construct the model
experimentInfo = dl.ex_ridme(tau1,tau2, pathways=[1])
Vmodel = dl.dipolarmodel(t,r,Bmodel=dl.bg_strexp, experiment =experimentInfo)

# Fit the model to the data
results = dl.fit(Vmodel,Vexp)

# Print results summary
print(results)

#%%

# Extract fitted dipolar signal
Vfit = results.model

# Extract fitted distance distribution
Pfit = results.P
Pci95 = results.PUncert.ci(95)
Pci50 = results.PUncert.ci(50)

# Extract the unmodulated contribution
Bfcn = dl.dipolarbackgroundmodel(experimentInfo,dl.bg_strexp)
Bfit = results.P_scale*results.evaluate(Bfcn,t)
Bci = results.P_scale*results.propagate(Bfcn,t).ci(95)

plt.figure(figsize=[6,7])
violet = '#4550e6'
plt.subplot(211)
# Plot experimental and fitted data
plt.plot(t,Vexp,'.',color='grey',label='Data')
plt.plot(t,Vfit,linewidth=3,color=violet,label='Fit')
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
