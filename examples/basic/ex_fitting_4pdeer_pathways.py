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
r = np.arange(2,5,0.025) # nm

# Construct the model
Vmodel = dl.dipolarmodel(t,r,npathways=3)

# Adjust the boundaries for the refocusing times
Vmodel.reftime1.set(par0=0.5, lb=0.0, ub=1.0) # Main pathway contribution
Vmodel.reftime2.set(par0=0.0, lb=0.0, ub=0.2) # Pathway refocusing at the start of the signal
Vmodel.reftime3.set(par0=4.5, lb=4.0, ub=5.0) # Pathway refocusing at the end of the signal

# Fit the model to the data
results = dl.fit(Vmodel,Vexp)

#%%

# Extract fitted dipolar signal
Vfit = results.model
Vci = results.modelUncert.ci(95)

# Extract fitted distance distribution
Pfit = results.P
Pci95 = results.PUncert.ci(95)
Pci50 = results.PUncert.ci(50)


plt.figure(figsize=[6,9])
violet = '#4550e6'
plt.subplot(311)
# Plot experimental data
plt.plot(t,Vexp,'.',color='grey',label='Data')
# Plot the fitted signal 
plt.plot(t,Vfit,linewidth=3,color=violet,label='Fit')
plt.fill_between(t,Vci[:,0],Vci[:,1],color=violet,alpha=0.3)
plt.legend(frameon=False,loc='best')
plt.xlabel('Time $t$ (μs)')
plt.ylabel('$V(t)$ (arb.u.)')

plt.subplot(312)
lams = [results.lam1, results.lam2, results.lam3]
reftimes = [results.reftime1, results.reftime2, results.reftime3]
for n,(lam,reftime) in enumerate(zip(lams,reftimes)):
    Vpath = lam*dl.dipolarkernel(t-reftime,r)@Pfit
    plt.plot(t,Vpath,linewidth=3,label=f'Pathway #{n+1}')
plt.legend(frameon=False,loc='best')
plt.xlabel('Time $t$ (μs)')
plt.ylabel('$V(t)$ (arb.u.)')

# Plot the distance distribution
plt.subplot(313)
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
