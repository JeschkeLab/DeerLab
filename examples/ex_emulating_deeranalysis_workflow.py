# %% [markdown]
"""
Emulating the DeerAnalysis workflow
===================================

This example shows how to reproduce the type of workflow implemented in
DeerAnalysis, using DeerLab functions. This kind of analysis workflow is 
outdated and not recommended for routine or accurate data analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %% [markdown]
# Generating a dataset
#---------------------
#
# For this example we will simulate a simple 4pDEER signal

# Parameters
t = np.linspace(-0.1,3,250) # µs
rtrue = np.linspace(1,7,200) # nm
Ptrue = dl.dd_gauss3(rtrue,[4.5, 0.35, 0.4, 3, 0.25, 0.3, 4, 0.4, 0.5])
lam = 0.3
conc = 180 # µM

# Simulate an experimental signal with some scale and phase
Bmodel = lambda t, lam: dl.bg_hom3d(t,conc,lam)
K = dl.dipolarkernel(t,rtrue,mod=lam,bg=Bmodel)
V = K@Ptrue*np.exp(1j*np.pi/16) # add a phase shift 
rnoise = dl.whitegaussnoise(t,0.01,seed=1) # real-component noise 
inoise = 1j*dl.whitegaussnoise(t,0.01,seed=2) # imaginary-component noise 
V = V + rnoise + inoise # complex-valued noisy signal
V = V*3e6 # add an arbitrary amplitude scale

plt.plot(t,V.real,'.',t,V.imag,'.'),
plt.xlabel('t (µs)')
plt.ylabel('V')
plt.grid(alpha=0.3)
plt.legend(['real','imag'])
plt.tight_layout()
plt.show()
# %% [markdown]
# DeerAnalysis workflow
# ---------------------
#

# Pre-processing
V = dl.correctphase(V)
t = dl.correctzerotime(V,t)
V = V/max(V)

# Distance axis estimation
r = dl.time2dist(t)

# Background fit
tstart = 1.0 # background fit start, in µs
mask = t>tstart
def Bmodel(par):
    lam,kappa,d = par # unpack parameters
    B = (1-lam)*dl.bg_strexp(t[mask],[kappa,d])
    return B

#       lam     k   d
par0 = [0.5,   0.5, 3]
lb   = [0.1,    0,  1]
ub   = [1,      5,  6]
fit = dl.fitparamodel(V[mask],Bmodel,par0,lb,ub,fitscale=False)

lamfit,kappa,d = fit.param
Bfit = dl.bg_strexp(t,[kappa,d])

# Background "correction" by division
Vcorr = (V/Bfit - 1 + lamfit)/lamfit

# Tikhonov regularization using the L-curve criterion
K = dl.dipolarkernel(t,r)
fit = dl.fitregmodel(Vcorr,K,r,'tikhonov','lr',)
Pfit = fit.P

# %% [markdown]
# Plots
# -----

plt.subplot(311)
plt.plot(t,V,'k.',t,(1-lamfit)*Bfit,'r',linewidth=1.5)
plt.xlabel('t (µs)')
plt.ylabel('V')
plt.legend(['data','$(1-\lambda)B_{fit}$'])

plt.subplot(312)
plt.plot(t,Vcorr,'k.',t,K@Pfit,'r',linewidth=1.5)
plt.xlabel('t (µs)')
plt.ylabel('V')
plt.legend(['corrected data','fit'])

plt.subplot(313)
plt.plot(rtrue,Ptrue,'k',r,Pfit,'r',linewidth=1.5)
plt.xlabel('r (nm)')
plt.ylabel('P (nm⁻¹)')
plt.legend(['truth','fit'])
plt.tight_layout()
plt.show()

# %%
