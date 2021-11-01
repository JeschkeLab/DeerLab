# %% [markdown]
"""
Basic fitting of a 5-pulse DEER signal
====================================================================

This example shows how to model and fit a 5-pulse DEER signal, including
the typically present additional dipolar pathway.  

Now, the simple 5pDEER models contain 3 additional parameters compared to 4pDEER (due
to the additional dipolar pathway present in the signal). However, the
refocusing time of the second dipolar pathway is very easy to constrain
and strongly helps stabilizing the fit. 
 
This pathway can even estimated visually from the signal or estimated from the 
pulse sequence timings. Thus, we can strongly constraint this parameters while leaving the
pathway amplitudes pretty much unconstrained.
""" 
# %%
# Import required packages
import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl


# %% 

# Load the experimental data
t,Vexp = np.load('../data/example_5pdeer_#1.npy')

# Distance vector
r = np.linspace(2,5,200) # nm

# Construct dipolar model with two dipolar pathways
Vmodel = dl.dipolarmodel(t,r,npathways=2)

# The refocusing time of the second pathway can be well estimated by visual inspection
Vmodel.reftime2.set(lb=3, ub=4, par0=3.5)

# Fit the model to the data
fit = dl.fit(Vmodel,Vexp)

# %%

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
Bfcn = lambda lam1,lam2,reftime1,reftime2,conc: scale*(1-lam1-lam2)*dl.bg_hom3d(t-reftime1,conc,lam1)*dl.bg_hom3d(t-reftime2,conc,lam2)
Bfit = Bfcn(fit.lam1,fit.lam2,fit.reftime1,fit.reftime2,fit.conc)
Bci = fit.propagate(Bfcn).ci(95)

plt.figure(figsize=[6,7])
plt.subplot(211)
# Plot experimental data
plt.plot(t,Vexp,'.',color='grey',label='Data')
# Plot the fitted signal 
plt.plot(t,Vfit,linewidth=3,label='Fit')
plt.fill_between(t,Vci[:,0],Vci[:,1],alpha=0.3)
plt.plot(t,Bfit,'--',linewidth=3,label='Unmodulated contribution')
plt.fill_between(t,Bci[:,0],Bci[:,1],alpha=0.3)
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
