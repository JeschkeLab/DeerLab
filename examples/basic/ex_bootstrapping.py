# %%
""" 
Bootstrapped confidence intervals in routine analysis
-------------------------------------------------------------------

How to obtain bootstrapped confidence intervals for simple routine operations.

Unless specified otherwise, the function ``fit`` will return asymptotic confidence intervals based on the covariance matrix 
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

# File location
path = '../data/'
file = 'example_4pdeer_1.DTA'

# Experimental parameters
tau1 = 0.3      # First inter-pulse delay, μs
tau2 = 4.0      # Second inter-pulse delay, μs
deadtime = 0.1  # Acquisition deadtime, μs

# Load the experimental data
t,Vexp = dl.deerload(path + file)

# Pre-processing
Vexp = dl.correctphase(Vexp) # Phase correction
Vexp = Vexp/np.max(Vexp)     # Rescaling (aesthetic)
t = t + deadtime             # Account for deadtime

# Distance vector
r = np.linspace(2,5,100) # nm

# Construct the model
Vmodel = dl.dipolarmodel(t,r, experiment = dl.ex_4pdeer(tau1,tau2, pathways=[1]))

# Fit the model to the data
results = dl.fit(Vmodel,Vexp,bootstrap=20)

# In this example, just for the sake of time, we will just use 20 bootstrap samples.  

# Print results summary
print(results)

#%%

# Extract fitted dipolar signal
Vfit = results.model
Vci = results.modelUncert.ci(95)

# Extract fitted distance distribution
Pfit = results.P
Pci95 = results.PUncert.ci(95)
Pci50 = results.PUncert.ci(50)

# Extract the unmodulated contribution
Bfcn = lambda mod,conc: results.P_scale*(1-mod)*dl.bg_hom3d(t,conc,mod)
Bfit = Bfcn(results.mod,results.conc)
Bci = results.propagate(Bfcn).ci(95)

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
