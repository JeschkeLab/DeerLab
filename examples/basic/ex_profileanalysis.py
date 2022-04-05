# %%
""" 
Profile analysis in routine analysis
-------------------------------------------------------------------

How to obtain objective function profiles for the non-linear parameters of a model. 
""" 

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %%

# Load the experimental data
t,Vexp = np.load('../data/example_4pdeer_#1.npy')

# Distance vector
r = np.arange(2,5,0.05) # nm

# Construct the model
Vmodel = dl.dipolarmodel(t,r)

# Fit the model to the data
results = dl.fit(Vmodel,Vexp)

# Print results summary
print(results)

# Compute uncertainty with the likelihood profile method for the spin concentration and modulation depth parameters
profile_uq = dl.profile_analysis(Vmodel,Vexp,samples=20, parameters=['conc','mod']) 

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

plt.figure(figsize=[9,7])
violet = '#4550e6'
plt.subplot(221)
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
plt.subplot(222)
plt.plot(r,Pfit,linewidth=3,label='Bootstrap median',color=violet)
plt.fill_between(r,Pci95[:,0],Pci95[:,1],alpha=0.3,color=violet,label='95%-Conf. Inter.',linewidth=0)
plt.fill_between(r,Pci50[:,0],Pci50[:,1],alpha=0.5,color=violet,label='50%-Conf. Inter.',linewidth=0)
plt.legend(frameon=False,loc='best')
plt.autoscale(enable=True, axis='both', tight=True)
plt.xlabel('Distance $r$ (nm)')
plt.ylabel('$P(r)$ (nm$^{-1}$)')

for n,param in enumerate(['conc','mod']):
    plt.subplot(2,2,2+n+1)
    profile = profile_uq[param].profile
    threshold = profile_uq[param].threshold(0.95)
    plt.plot(profile['x'],profile['y'] - threshold,'-',color=violet, linewidth=3, label='Profile')
    plt.hlines(0,min(profile['x']),max(profile['x']), linewidth=3, linestyles='--',color='grey',alpha=0.6, label='Threshold')
    plt.autoscale(True,'both',tight=True)
    plt.ylim([7*np.min(profile['y'] - threshold),1e8])  
    plt.xlabel(f'{getattr(Vmodel,param).description} ({getattr(Vmodel,param).unit})')
    plt.ylabel('Profile objective function')
    plt.legend(frameon=False,loc='best')

plt.tight_layout()
plt.show()

# %%
