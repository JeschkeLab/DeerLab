# %% [markdown]
""" 
Analysis of a 5-pulse DEER signal with multiple dipolar pathways
-------------------------------------------------------------------------

Fit a 5-pulse DEER with multiple dipolar pathways and display the individual pathway contributions. 
""" 

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl


# %%

# File location
path = '../data/'
file = 'example_5pdeer_3.DTA'

# Experimental parameters (reversed 5pDEER)
tau1 = 3.7              # First inter-pulse delay, μs
tau2 = 3.5              # Second inter-pulse delay, μs
tau3 = 0.3              # Third inter-pulse delay, μs
tmin = 0.1              # Start time, μs

# Load the experimental data
t,Vexp = dl.deerload(path + file)
Vexp = dl.correctphase(Vexp)   # Phase correction
Vexp = Vexp/np.max(Vexp)       # Rescaling (aesthetic)
t = t - t[0]                     # Account for zerotime
t = t + tmin   
 
# Distance vector
r = np.arange(2.5,5.5,0.025) # nm

# Construct the model
experimentInfo = dl.ex_rev5pdeer(tau1,tau2,tau3, pathways=[1,2,3,5])
Vmodel = dl.dipolarmodel(t,r,experiment=experimentInfo)

# Fit the model to the data
results = dl.fit(Vmodel,Vexp)

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

plt.figure(figsize=[8,6])
violet = '#4550e6'
green = '#3cb4c6' 
red = '#f84862'
plt.subplot(221)
# Plot experimental data
plt.plot(t,Vexp,'.',color='grey',label='Data')
# Plot the fitted signal 
plt.plot(t,Vfit,linewidth=3,color=violet,label='Fit')
plt.legend(frameon=False,loc='best')
plt.xlabel('Time $t$ (μs)')
plt.ylabel('$V(t)$ (arb.u.)')

plt.subplot(222)
labels = [1,2,3,5]
lams = [results.lam1, results.lam2, results.lam3, results.lam5]
reftimes = [results.reftime1, results.reftime2, results.reftime3, results.reftime5]
colors= ['tab:blue','tab:orange', red, green] 
Vinter = results.P_scale*results.evaluate(dl.dipolarbackgroundmodel(experimentInfo),t)
for (lam,reftime,color,label) in zip(lams,reftimes,colors,labels):
    Vpath = (1-np.sum(lams) + lam*dl.dipolarkernel(t-reftime,r)@Pfit)*Vinter
    plt.plot(t,Vpath,linewidth=3,label=f'Pathway #{label}',color=color)
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
