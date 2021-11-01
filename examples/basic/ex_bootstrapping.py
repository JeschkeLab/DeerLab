# %%
""" 
Bootstrapped confidence intervals in routine analysis
-------------------------------------------------------------------

How to obtain bootstrapped confidence intervals for simple routine operations.

Unless specified otherwise, the function ``fitmodel`` will return asymptotic confidence intervals based on the covariance matrix 
of the objective function used to fit the data. These are quick to calculate and therefore very comfortable for quick estimates of
the uncertainty during routine analysis or testing. 

However, for publication-level analysis, these confidence intervals might be inaccurate. It is strongly recommended to employ bootstrapped 
confidence intervals to get accurate estimates of the uncertainty. 
Conviniently, ``fit`` integrates bootstrapping to make it accessible via the keyword argument ``bootstrap`` which specifies the
number of samples to analyze to estimate the uncertainty. The larger this number, the more accurate 
the confidence intervals but the longer the analysis will be. The standard for publication is typically 1000 samples. 
""" 

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl


# %%

# Load the experimental data
t,Vexp = np.load('../data/example_4pdeer_#1.npy')

# Distance vector
r = np.linspace(2,5,100) # nm

# Construct the model
Vmodel = dl.dipolarmodel(t,r)

# Fit the model to the data
fit = dl.fit(Vmodel,Vexp,bootstrap=20)

# In this example, just for the sake of time, we will just use 20 bootstrap samples.  

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

# Extract the unmodulated contribution
Bfcn = lambda mod,conc: scale*(1-mod)*dl.bg_hom3d(t,conc,mod)
Bfit = Bfcn(fit.mod,fit.conc)
Bci = fit.propagate(Bfcn).ci(95)

plt.figure(figsize=[6,7])
violet = '#4550e6'
plt.subplot(211)
# Plot experimental data
plt.plot(t,Vexp,'.',color='grey',label='Data')
# Plot the fitted signal 
plt.plot(t,Vfit,linewidth=3,label='Bootstrap median',color=violet)
plt.fill_between(t,Vci[:,0],Vci[:,1],alpha=0.3,color=violet)
plt.plot(t,Bfit,'--',linewidth=3,color=violet,label='Unmodulated contribution')
plt.fill_between(t,Bci[:,0],Bci[:,1],alpha=0.3,color=violet)
plt.legend(frameon=False,loc='best')
plt.xlabel('Time $t$ (μs)')
plt.ylabel('$V(t)$ (arb.u.)')
# Plot the distance distribution
plt.subplot(212)
plt.plot(r,Pfit,linewidth=3,label='Bootstrap median',color=violet)
plt.fill_between(r,Pci95[:,0],Pci95[:,1],alpha=0.3,color=violet,label='95%-Conf. Inter.',linewidth=0)
plt.fill_between(r,Pci50[:,0],Pci50[:,1],alpha=0.5,color=violet,label='50%-Conf. Inter.',linewidth=0)
plt.legend(frameon=False,loc='best')
plt.autoscale(enable=True, axis='both', tight=True)
plt.xlabel('Distance $r$ (nm)')
plt.ylabel('$P(r)$ (nm$^{-1}$)')
plt.tight_layout()
plt.show()

# %%
