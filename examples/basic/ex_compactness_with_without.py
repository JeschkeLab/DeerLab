# %% [markdown]
""" 
Comparing the analysis with and without compactness criterion
-------------------------------------------------------------------------

This examples shows how to quickly compare the effects of the compactness criterion
on the results of the analysis of a 4-pulse DEER signal. 
""" 

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl


# %%

# File location
path = dl.__path__[0] + '/../examples/data/'
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

# Artifically truncate the signal 
Vexp = Vexp[t<=3]
t = t[t<=3]

# Distance vector
r = np.arange(1.5,7,0.05) # nm

# Construct the model
Vmodel = dl.dipolarmodel(t,r, experiment=dl.ex_4pdeer(tau1,tau2, pathways=[1]))
compactness = dl.dipolarpenalty(Pmodel=None,r=r,type='compactness')

# Fit the model to the data with compactness criterion
results_with = dl.fit(Vmodel,Vexp,penalties=compactness)
print(results_with)

# Fit the model to the data without compactness criterion
results_without = dl.fit(Vmodel,Vexp)
print(results_without)

#%%

# Define colors 
green = '#3cb4c6' 
red = '#f84862'
colors = [green,red]
plt.figure(figsize=[10,7])
for n,results in enumerate([results_with, results_without]):

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

    plt.subplot(2,2,n+1)
    # Plot experimental and fitted data
    plt.plot(t,Vexp,'.',color='grey',label='Data')
    plt.plot(t,Vfit,linewidth=3,color=colors[n],label='Fit')
    plt.fill_between(t,Vci[:,0],Vci[:,1],color=colors[n],alpha=0.3)
    plt.plot(t,Bfit,'--',linewidth=3,color=colors[n],label='Unmodulated contribution')
    plt.fill_between(t,Bci[:,0],Bci[:,1],color=colors[n],alpha=0.3)
    plt.legend(frameon=False,loc='best')
    plt.xlabel('Time $t$ (μs)')
    plt.ylabel('$V(t)$ (arb.u.)')
    plt.title('With compactness criterion' if n==0 else 'Without compactness criterion')
    # Plot the distance distribution
    plt.subplot(2,2,n+3)
    plt.plot(r,Pfit,color=colors[n],linewidth=3,label='Fit')
    plt.fill_between(r,Pci95[:,0],Pci95[:,1],alpha=0.3,color=colors[n],label='95%-Conf. Inter.',linewidth=0)
    plt.fill_between(r,Pci50[:,0],Pci50[:,1],alpha=0.5,color=colors[n],label='50%-Conf. Inter.',linewidth=0)
    plt.legend(frameon=False,loc='best')
    plt.autoscale(enable=True, axis='both', tight=True)
    plt.xlabel('Distance $r$ (nm)')
    plt.ylabel('$P(r)$ (nm$^{-1}$)')
plt.tight_layout()
plt.show()

# %%
