# %% [markdown]
"""
Analyzing the Pake pattern of a dipolar signal
============================================================================

A very basic example for displaying the dipolar spectrum (Pake pattern) of a given dipolar signal.
""" 
# %%
import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %% [markdown]
# Generate a dipolar signal
# -------------------------
# Let's start by simulating a dipolar signal with some background and noise.

# %%
# Prepare components
t = np.linspace(0,5,400)   # µs
r = np.linspace(2,5,100)   # nm
P = dl.dd_gauss2(r,[3.5, 0.1, 0.2, 4, 0.05, 0.8])
B = dl.bg_exp(t,0.2)
lam = 0.3
K = dl.dipolarkernel(t,r,mod=lam,bg=B)
V = K@P + dl.whitegaussnoise(t,0.005,seed=0)

# Plot
plt.plot(t,V,'k.')
plt.grid(alpha=0.3)
plt.xlabel('Time (µs)')
plt.ylabel('V')
plt.tight_layout()
plt.show()
# %% [markdown]
# Prepare the signal
# ------------------
# Since experimental dipolar signals contain the background, this must be fitted 
# removed prior to Fourier transform.
# 
# First we proceed to fit the background function using some time-domain parametric 
# model. In this example we will use an exponential function (``bg_exp``). 
# Using the ``fitparamodel`` function we obtain the fitted background as well as 
# the fitted modulation depth.

# %%

tstart = 3 # Time to start fitting background, in µs
mask = t>tstart
# Model for the background component (1-lambda)*B
def Bmodel(par):
    lam,kappa = par 
    B = (1 - lam)*dl.bg_exp(t[mask],kappa)
    return B

# Fit the background function
fit = dl.fitmodel(V,t,r,'P',dl.bg_exp,dl.ex_4pdeer)
Bfit = fit.B
lam = fit.exparam
kappa = fit.bgparam

# %% [markdown]
# Now we can use these fitted parameters to isolate the dipolar evolution function 
# from the primary data. Removal of the background via division leads to a noise 
# increase at later times and thus to an approximation ``Vcorr`` of the real dipolar 
# evolution function.

# %%

# "Correct" for the background and modulation depth
Vcorr = (V/Bfit - (1 - lam))/lam

plt.plot(t,Vcorr,'k.')
plt.grid(alpha=0.3)
plt.xlabel('Time (µs)')
plt.ylabel('V')
plt.tight_layout()
plt.show()
# %% [markdown]
# Computing the dipolar spectrum
# --------------------------------
#
# Now that the signal has the appropiate structure for Fourier transform it, 
# we can call the ``fftspec`` function to obtain the dipolar spectrum.

# %%

# Compute spectrum
nu,pake = dl.fftspec(Vcorr,t,apodization=False)
 
 # %% [markdown]
# In order to avoid truncation ripples in the Fourier spectrum and at the same 
# time to compensate for the increase of noise, we recommend the use of apodization 
# using the appropiate option in ``fftspec``.

# %%
# Compute spectrum with apodization
nuapo,pakeapo = dl.fftspec(Vcorr,t,apodization=False,mode='real')

# Plot results
plt.plot(nu,pake,'k',nuapo,pakeapo,'b',linewidth=1.5)
plt.grid(alpha=0.3)
plt.xlim([-10, 10])
plt.xlabel('Frequency (MHz)')
plt.ylabel('Intensity (arb.u.)')
plt.legend(['Raw','Apodized'])
plt.tight_layout()
plt.show()
# %% [markdown]
# We do not need to worry about the zero-filling since ``fftspec`` takes care 
# of setting it to twice the amount of points in the signal, to preserve all information. 
# Adding more points will artificially increase the resolution of the Pake pattern.
# The improvement will only be visual as no further information can be gained 
# from additional zero-filling.


# %%
