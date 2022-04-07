# %% [markdown]
"""
Basic fitting of a 5-pulse DEER signal
====================================================================

This example shows how to model and fit a 5-pulse DEER signal, including
the typically present additional dipolar pathway.  

Now, the simple 5pDEER models contain 3 additional parameters compared to 4pDEER (due
to the additional dipolar pathway present in the signal). However, the
refocusing time of the second dipolar pathway is very easy to constrain
and strongly helps stabilizing the results. 
 
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

# File location
path = dl.__path__[0] + '/../examples/data/'
file = 'example_5pdeer_#1.DTA'

# Experimental parameters (reversed 5pDEER)
tau1 = 3.9               # First inter-pulse delay, μs
tau2 = 3.7               # Second inter-pulse delay, μs
tau3 = 0.5               # Third inter-pulse delay, μs
deadtime = 0.3           # Acquisition deadtime, μs

# Load the experimental data
t,Vexp = dl.deerload(path + file)
Vexp = dl.correctphase(Vexp)    # Phase correction
Vexp = Vexp/np.max(Vexp)         # Rescaling (aesthetic)
t = t + deadtime       # Account for deadtime

# Distance vector
r = np.arange(3,5,0.025) # nm

# Construct dipolar model with two dipolar pathways
experimentInfo = dl.ex_rev5pdeer(tau1, tau2, tau3, pathways=[1,2])
Vmodel = dl.dipolarmodel(t,r, experiment=experimentInfo)

# Fit the model to the data
results = dl.fit(Vmodel,Vexp)

# Print results summary
print(results)

# %%

# Extract fitted dipolar signal
Vfit = results.model
Vci = results.modelUncert.ci(95)

# Extract fitted distance distribution
Pfit = results.P
Pci95 = results.PUncert.ci(95)
Pci50 = results.PUncert.ci(50)
Pfit =  Pfit

# Extract the unmodulated contribution
Bfcn = lambda lam1,lam2,reftime1,reftime2,conc: results.P_scale*(1-lam1-lam2)*dl.bg_hom3d(t-reftime1,conc,lam1)*dl.bg_hom3d(t-reftime2,conc,lam2)
Bfit = Bfcn(results.lam1,results.lam2,results.reftime1,results.reftime2,results.conc)
Bci = results.propagate(Bfcn).ci(95)

plt.figure(figsize=[6,7])
violet = '#4550e6'
plt.subplot(211)
# Plot experimental data
plt.plot(t,Vexp,'.',color='grey',label='Data')
# Plot the fitted signal 
plt.plot(t,Vfit,linewidth=3,color=violet,label='Fit')
plt.fill_between(t,Vci[:,0],Vci[:,1],color=violet,alpha=0.3)
plt.plot(t,Bfit,'--',linewidth=3,color=violet,alpha=0.5,label='Unmodulated contribution')
plt.fill_between(t,Bci[:,0],Bci[:,1],color=violet,alpha=0.1)
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
