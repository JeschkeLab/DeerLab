# %% [markdown]
""" 
Distance restraints from 4-pulse DEER data, non-parametric distribution
-----------------------------------------------------------------------

How to fit a simple 4-pulse DEER signal and derive distance restraints from 
the fitted non-parametric distance distribution.
""" 

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %% [markdown]
# Uncomment and use the following lines if you have experimental data:
#
# t,Vexp = dl.deerload('my\path\4pdeer_data.DTA')
# Vexp = dl.correctphase(Vexp)
# t = dl.correctzerotime(Vexp,t)
#

# %% [markdown]#
# In this example we will use simulated data instead:

#%% 

# Define a function that generates synthetic data
def generatedata():
    t = np.linspace(-0.1,4,250)        # time axis, µs
    r = np.linspace(1,6,200)           # distance axis, nm
    param = [3, 0.1, 0.2, 3.5, 0.1, 0.65, 3.8, 0.05, 0.15] # parameters for three-Gaussian model
    P = dl.dd_gauss3(r,param)          # model distance distribution
    lam = 0.5                          # modulation depth
    B = dl.bg_hom3d(t,300,lam)         # background decay
    K = dl.dipolarkernel(t,r,mod=lam,bg=B)    # kernel matrix
    Vexp = K@P + dl.whitegaussnoise(t,0.03,seed=0)
    return t, Vexp

t, Vexp = generatedata()

# %%

# Run fit
r = np.linspace(1,6,200)
fit = dl.fitmodel(Vexp,t,r,'P',dl.bg_hom3d,dl.ex_4pdeer)
fit.plot();

# %% [markdown]
# Now that we have a fit of the distance distribution we can obtain distance restraints in the form of different statistical 
# descriptors such as the mean distance and standard deviation of distances. 
# While we could calculate this manually, DeerLab provides a convenient function ``diststats`` which will automatically compute
# these for you and even propagate the uncertainty in the distributions to those values to get confidence intervals on the restraints. 
#
# If we specify the ``verbose`` option, we can get a complete summary of all statistical descriptors of the fitted distributions 
# including 95% confidence intervals.

#%% 

# Get printed summary of all statistical descriptors available with confidence intervals
estimators,uq = dl.diststats(r,fit.P,fit.Puncert,verbose=True)

# %% [markdown]
# However, if you are just interested in specific quantities to use as restraints, you can extract them from the returned dictionary.
# For example, let's just get the mean distance and standard deviation of the fitted distribution.

#%% 

# Mean distance
rmean = estimators['mean']
rmean_ci = uq['mean'].ci(95)
# Standard deviation
r_std = estimators['std']
r_std_ci = uq['std'].ci(95)

# Print out the results
print(f'Mean distance: {rmean:.3f} ({rmean_ci[0]:.3f}-{rmean_ci[1]:.3f}) nm')
print(f'Standard deviation: {r_std:.3f} ({r_std_ci[0]:.3f}-{r_std_ci[1]:.3f}) nm')

# %% [markdown]
# For display, you can plot the mean distance with its confidence intervals without further calculations.

# %%

# sphinx_gallery_thumbnail_number = 2

# Plot distribution and confidence bands
Pci95 = fit.Puncert.ci(95)
Pci50 = fit.Puncert.ci(50)
plt.plot(r,fit.P,linewidth=2,label='Distance distribution fit')
plt.fill_between(r,Pci95[:,0],Pci95[:,1],color='tab:blue',alpha=0.1)
plt.fill_between(r,Pci50[:,0],Pci50[:,1],color='tab:blue',alpha=0.3)

# Plot mean distance and confidence interval
plt.vlines(rmean,0,max(Pci95[:,1]),color='tab:red',linestyles='dotted',linewidth=3,label='Mean distance')
plt.fill_between(rmean_ci,0,max(Pci95[:,1]),color='tab:red',alpha=0.3)

plt.legend()

plt.ylim([0,max(Pci95[:,1])])
plt.xlabel('r (nm)')
plt.ylabel('P (nm⁻¹)')

plt.tight_layout() 
plt.show()
# %%
