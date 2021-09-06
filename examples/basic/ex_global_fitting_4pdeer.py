# %% [markdown]
""" 
Global fitting of multiple 4-pulse DEER signals, non-parametric distribution
-----------------------------------------------------------------------------

How to fit multiple 4-pulse DEER signals to a model with a non-parametric
distribution and a homogeneous background.
""" 
import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl
# Use the seaborn style for nicer plots
from seaborn import set_theme
set_theme()

# %%

# Load the experimental data
t1,Vexp1 = np.load('../data/example_4pdeer_#2.npy')
t2,Vexp2 = np.load('../data/example_4pdeer_#3.npy')
t3,Vexp3 = np.load('../data/example_4pdeer_#4.npy')

# Put the datasets into lists
ts = [t1,t2,t3]
Vs = [Vexp1,Vexp2,Vexp3]

# Distance vector
r = np.linspace(1.5,5.5,150) # nm

# Construct the dipolar models for the individual signals
V1model = dl.dipolarmodel(ts[0],r)
V2model = dl.dipolarmodel(ts[1],r)
V3model = dl.dipolarmodel(ts[2],r)

# Make the global model by joining the individual models
globalmodel = dl.expand(V1model,V2model,V3model)

# Link the distance distribution into a global parameter 
globalmodel = dl.link(globalmodel,
        P=[globalmodel.P_1,globalmodel.P_2,globalmodel.P_3])

# Fit the model to the data
fit = dl.fit(globalmodel,Vs)

# %%

plt.figure(figsize=[10,7])
for n in range(len(fit.model)):

    # Extract fitted dipolar signal
    Vfit = fit.model[n]
    Vci = fit.modelUncert[n].ci(95)

    # Extract fitted distance distribution
    Pfit = fit.P
    scale = np.trapz(Pfit,r)
    Pci95 = fit.PUncert.ci(95)/scale
    Pci50 = fit.PUncert.ci(50)/scale
    Pfit =  Pfit/scale

    plt.subplot(2,3,n+1)
    # Plot experimental data
    plt.plot(ts[n],Vs[n],'.',color='grey',label='Data')
    # Plot the fitted signal 
    plt.plot(ts[n],Vfit,linewidth=3,label='Fit')
    plt.fill_between(ts[n],Vci[:,0],Vci[:,1],alpha=0.3)
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
