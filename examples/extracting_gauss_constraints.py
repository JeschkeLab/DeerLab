# %% [markdown]
"""
Fitting Gaussians to a non-parametric distance distribution
=======================================================================

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

t = np.linspace(0,5,250)   # time axis, µs
r = np.linspace(1,7,200)   # distance axis, nm
P = dl.dd_gauss3(r,[4.5, 0.35, 0.4, 3, 0.25, 0.3, 4, 0.4, 0.5])  # distance distribution
lam = 0.3                  # modulation depth
conc = 80                  # spin concentration, µM

B = dl.bg_hom3d(t,conc,lam)                   # background
K = dl.dipolarkernel(t,r,mod=lam,bg=B)               # kernel matrix
V = K@P + dl.whitegaussnoise(t,0.01,seed=0)   # DEER trace, with added noise

# %% [markdown]
# Fit the dipolar signal
# ----------------------
# First, we fit the non-parametric distance distribution using ``fitmodel()``
# %%
fit = dl.fitmodel(V,t,r,'P',dl.bg_exp,dl.ex_4pdeer)
fit.plot()
plt.show() 

# %% [markdown]
# Fit Gaussians to the distance distribution
# ------------------------------------------
# Next, we fit a multi-Gauss distribution to the fitted non-parametric
# distribution. We can do this by using the ``fitparamodel()`` function (in
# this example, fitting a two-Gauss model). 
#
# However, in order to get the correct uncertainty quantification, we need
# to specify the covariance matrix of the fitted distribution.
# ``fitparamodel()`` can then use that information to propagate the error in
# ``Pfit`` to the Gauss constraints that we then fit.

# %%

# From the fit results, extract the distribution and the covariance matrix
Pfit = fit.P
Pfit_uq = fit.Puncert
Pfit_covmat = Pfit_uq.covmat

# %% [markdown]
#Fit a 2-Gauss model to the fitted parameter-free distribution:
#
#    - ``parfit```: will contain the Gaussian constraints
#    - ``PGauss```: the corresponding distribution
#    - ``paruq```: the uncertainty quantification of our constraints

# %%
Pmodel = lambda p: dl.dd_gauss2(r,p)
# Get information on the model
info = dl.dd_gauss2()
par0 = info['Start']
lb = info['Lower']
ub = info['Upper']

# Fit the Gaussians
fit = dl.fitparamodel(Pfit,Pmodel,par0,lb,ub,covmatrix=Pfit_covmat,fitscale=False)

# Extract the fit results
parfit = fit.param
paruq = fit.paramUncert
PGauss = dl.dd_gauss2(r,parfit)

# Extract the 95%-confidence intervals...
par95 = paruq.ci(95)
# ... and print the results
print('\nGaussian components:')
info = dl.dd_gauss2()
for i in range(len(parfit)):
    print(f'  parfit[{i}] = {parfit[i]:2.2f} ({par95[i,0]:2.2f}, {par95[i,1]:2.2f}) {info["Parameters"][i]}')

# Now propagate the error of the Gaussian parameters to the distribution
lb = np.zeros_like(r) # non-negativity constraint
PGauss_uq = paruq.propagate(lambda par: dl.dd_gauss2(r,par),lb)
PGauss95 = PGauss_uq.ci(95)

# %%

# Plot the fitted constraints model on top of the non-parametric case
plt.plot(r,Pfit,'r',linewidth=1.5,label='non-param. fit')
plt.fill_between(r,Pfit_uq.ci(95)[:,0], Pfit_uq.ci(95)[:,1],facecolor='r',linestyle='None',alpha=0.2,label='95% CI')

plt.plot(r,PGauss,'b',linewidth=1.5,label='2-Gauss fit to nonparam. fit')
plt.fill_between(r,PGauss95[:,0], PGauss95[:,1],facecolor='b',linestyle='None',alpha=0.2,label='95% CI')

plt.xlabel('Distance (nm)')
plt.ylabel('P (nm⁻¹)')
plt.tight_layout()
plt.grid(alpha=0.3)
plt.legend()
plt.show()

# %%
