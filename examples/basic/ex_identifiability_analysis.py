# %%
""" 
Identifiability analysis of a 4-pulse DEER signal via the profile analysis method
----------------------------------------------------------------------------------

How to use the profile analysis method to asses whether the 4-pulse DEER model parameters 
are identifiable (unique) for a given experimental signal.  
""" 

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl
from deerlab.constants import D
green = '#3cb4c6' 
red = '#f84862'

# %%

# Load the experimental data
t,Vexp = np.load('../data/example_4pdeer_#1.npy')

Vexp = Vexp/np.max(Vexp)

# Truncate the signal
Vexp_truncated = Vexp[t<=2]
t_truncated = t[t<=2]

# Distance vector
r = np.arange(2,7,0.05) # nm

# Construct the model
Vmodel = dl.dipolarmodel(t,r)
Vmodel_truncated = dl.dipolarmodel(t_truncated,r)

# Compute uncertainty with the likelihood profile method for the spin concentration and modulation depth parameters
grids = {
    'conc': np.linspace(10,600,15),
    'mod': np.linspace(0.4,0.7,15),
}

profile_long = dl.profile_analysis(Vmodel,Vexp,samples=10, parameters=['conc','mod'], grids=grids) 
profile_short = dl.profile_analysis(Vmodel_truncated,Vexp_truncated,samples=10, parameters=['conc','mod'], grids=grids) 

#%%

plt.figure(figsize=[8,4])
for n,param in enumerate(['conc','mod']):
    plt.subplot(1,2,n+1)
    for profile_uq,color in zip([profile_long,profile_short],[green,red]):
        profile = profile_uq[param].profile
        threshold = profile_uq[param].threshold(0.95)
        plt.plot(profile['x'],profile['y'] - threshold,'-',color=color, linewidth=3)
        plt.hlines(0,min(profile['x']),max(profile['x']), linewidth=3, linestyles='--',color='grey',alpha=0.6)
    plt.autoscale(True,'both',tight=True)
    plt.ylim([1.5*np.min(profile['y'] - threshold),0.01])  
    plt.xlabel(f'{getattr(Vmodel,param).description} ({getattr(Vmodel,param).unit})')
    plt.ylabel('Profile objective function')
    plt.legend(['Profile (long trace)','Profile (short trace)','Threshold'],frameon=False,loc='best')

plt.tight_layout()
plt.show()

# %%
