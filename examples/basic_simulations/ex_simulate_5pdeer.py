# %%
"""
Simulating a two-pathway 5-pulse DEER signal
============================================================

An example on how to simulate a 5-pulse DEER dipolar signal.  

Specifically, we simulate a single-pathway 5-pulse DEER dipolar signal arising from 
a Gaussian distance distribution
"""

# Import the required libraries
import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl
violet = '#4550e6'

# %%

# Simulation parameters
tau1 = 3.5          # 1st experimental time delay, μs
tau2 = 4.2          # 2nd experimental time delay, μs
tau3 = 0.3          # 3rd experimental time delay, μs
conc = 150           # Spin concentration, μs
lam1 = 0.30         # Amplitude of dipolar pathway refocusing at t=tau3, μs    
lam2 = 0.15         # Amplitude of dipolar pathway refocusing at t=tau2, μs    
rmean = 4.0         # Mean distance, nm
rstd = 0.4          # Distance standard deviation, nm
rmin,rmax = 1.5,6   # Range of the distance axis, nm
Δr = 0.05           # Distance resolution, nm
Δt = 0.008          # Time resolution, μs 
deadtime = 0.1      # Acquisition deadtime, μs

# Experimental time vector
t = np.arange(deadtime,tau1+tau2+tau3,Δt)
# Distance vector 
r = np.arange(rmin,rmax,Δr)

experiment = dl.ex_rev5pdeer(tau1,tau2,tau3, pathways=[1,2])
reftime1, reftime2 = experiment.reftimes(tau1,tau2,tau3)

# Construct the dipolar signal model
Vmodel = dl.dipolarmodel(t,r,Pmodel=dl.dd_gauss, experiment=experiment) 
# Function for the scaled background
Vinter_fcn = lambda lam1,lam2,conc: (1-lam1-lam2)*dl.bg_hom3d(t-reftime1,conc,lam1)*dl.bg_hom3d(t-reftime2,conc,lam2)

# Simulate the signal with orientation selection
Vsim = Vmodel(mean=rmean, std=rstd, conc=conc, scale=1, lam1=lam1, lam2=lam2, reftime1=reftime1, reftime2=reftime2)

# Plot the simulated signal
plt.figure(figsize=[4,3])
plt.plot(t,Vsim,lw=2,label='V(t)',color=violet)
plt.plot(t,Vinter_fcn(lam1,lam2,conc),'--',color=violet,lw=2,label='(1-λ)$V_{inter}$')
plt.legend()
plt.xlabel('Time (μs)')
plt.ylabel('V(t)')
plt.tight_layout()
plt.show()

# %%
