# %%
"""
Simulating a 5-pulse DEER signal
============================================================

An example on how to simulate a 5-pulse DEER dipolar signal containing the
two main dipolar pathway contributions, using a Gaussian distance distribution.
"""

# Import the required libraries
import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %%

# Simulation parameters
tau1, tau2, tau3 = 3.5, 4.2, 0.3  # Inter-pulse delays, μs
Δt = 0.008              # Time increment, μs 
tmin = 0.1              # Start time, μs

rmean = 4.0             # Mean distance, nm
rstd = 0.4              # Distance standard deviation, nm
rmin, rmax = 1.5, 6     # Range of the distance axis, nm
Δr = 0.05               # Distance resolution, nm

conc = 150              # Spin concentration, μM
lam1 = 0.30             # Amplitude of dipolar pathway refocusing at t=tau3
lam2 = 0.15             # Amplitude of dipolar pathway refocusing at t=tau2
V0 = 1                  # Overall echo amplitude

# Time vector
tmax = tau1+tau2+tau3
t = np.arange(tmin, tmax, Δt)
# Distance vector 
r = np.arange(rmin, rmax, Δr)

experiment = dl.ex_rev5pdeer(tau1, tau2, tau3, pathways=[1,2])
reftime1, reftime2 = experiment.reftimes(tau1, tau2, tau3)

# Construct the dipolar signal model
Vmodel = dl.dipolarmodel(t, r, Pmodel=dl.dd_gauss, experiment=experiment) 

# Simulate the signal with orientation selection
Vsim = Vmodel(mean=rmean, std=rstd, conc=conc, scale=V0, lam1=lam1, lam2=lam2, reftime1=reftime1, reftime2=reftime2)

# Scaled background (for plotting)
Vinter = V0*(1-lam1-lam2)*dl.bg_hom3d(t-reftime1,conc,lam1)*dl.bg_hom3d(t-reftime2,conc,lam2)

# Plot the simulated signal
violet = '#4550e6'
plt.figure(figsize=[4,3])
plt.plot(t, Vsim, lw=2, label='V(t)', color=violet)
plt.plot(t, Vinter, '--', color=violet, lw=2, label='(1-λ)$V_{inter}$')
plt.legend()
plt.xlabel('Time (μs)')
plt.ylabel('V(t)')
plt.tight_layout()
plt.show()

# %%
