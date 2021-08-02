# %% [markdown]
"""
Fitting Gaussians to a non-parametric distance distribution
============================================================================

This example shows how to fit Gaussians to a non-parametric distance
distribution obtained via Tikhonov regularization and how to calculate
the corresponding uncertainty.
""" 
# %%

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %%

# Simulate a 4pDEER signal

def simulatedata():
    r = np.linspace(1,7,200)   # distance axis, nm
    P = dl.dd_gauss3(r,[4.5, 0.35, 0.4, 3, 0.25, 0.3, 4, 0.4, 0.5])  # distance distribution
    lam = 0.3                  # modulation depth
    conc = 80                  # spin concentration, µM
    t = np.linspace(0,5,250)   # time axis, µs
    B = dl.bg_hom3d(t,conc,lam)                   # background
    K = dl.dipolarkernel(t,r,mod=lam,bg=B)        # kernel matrix
    V = K@P + dl.whitegaussnoise(t,0.01,seed=0)   # DEER trace, with added noise
    return t, V
    
t, V = simulatedata()

# %% [markdown]
# Fit the dipolar signal
#----------------------
# First, we fit the non-parametric distance distribution using ``fitmodel()``
# %%
r = np.linspace(1,7,200)
fit = dl.fitmodel(V,t,r,'P',dl.bg_exp,dl.ex_4pdeer)
fit.plot()
plt.show() 

# %% [markdown]
# Fit Gaussians to the distance distribution
# ------------------------------------------
# Next, we fit a multi-Gauss distribution to the fitted non-parametric
# distribution. We can do this by using the ``nlls()`` function (in
# this example, fitting a two-Gauss model). 
#
# However, in order to get the correct uncertainty quantification, we need
# to specify the covariance matrix of the fitted distribution.
# ``nlls()`` can then use that information to propagate the error in
# ``Pfit`` to the Gauss constraints that we then fit.

# %%

# From the fit results, extract the distribution and the covariance matrix
Pfit = fit.P
Pfit_uq = fit.Puncert

# %% [markdown]
#Fit a 2-Gauss model to the fitted parameter-free distribution:
#
#    - ``parfit```: will contain the Gaussian constraints
#    - ``PGauss```: the corresponding distribution
#    - ``paruq```: the uncertainty quantification of our constraints

# %%
Pmodel = lambda par: dl.dd_gauss2(r,par)

# Get information on the model
par0 = dl.dd_gauss2.start
lb = dl.dd_gauss2.lower
ub = dl.dd_gauss2.upper

# Fit the Gaussians
fit = dl.nlls(Pfit,Pmodel,par0,lb,ub,fitscale=False)

# Extract the fit results
parfit = fit.param
paruq = fit.paramUncert
PGauss = dl.dd_gauss2(r,parfit)

# Extract the 95%-confidence intervals...
par95 = paruq.ci(95)
# ... and print the results
print('\nGaussian components:')
for i in range(len(parfit)):
    print(f'  parfit[{i}] = {parfit[i]:2.2f} {dl.dd_gauss2.parameters[i]}')

# %%
# sphinx_gallery_thumbnail_number = 2

# Plot the fitted constraints model on top of the non-parametric case
plt.plot(r,Pfit,'r',linewidth=1.5,label='non-param. fit')
plt.fill_between(r,Pfit_uq.ci(95)[:,0], Pfit_uq.ci(95)[:,1],facecolor='r',linestyle='None',alpha=0.2,label=r'95% confidence intervals')

plt.plot(r,PGauss,'b',linewidth=1.5,label='2-Gauss fit to nonparam. fit')

plt.xlabel('Distance (nm)')
plt.ylabel('P (nm⁻¹)')
plt.tight_layout()
plt.grid(alpha=0.3)
plt.legend()
plt.show()

# %%
