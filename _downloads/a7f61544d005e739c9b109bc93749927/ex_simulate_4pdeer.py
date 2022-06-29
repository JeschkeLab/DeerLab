# %%
"""
Simulating a 4-pulse DEER signal
============================================================

An example on how to simulate a 4-pulse DEER dipolar signal.  

Specifically, we simulate a single-pathway 4-pulse DEER dipolar signal arising from 
a Gaussian distance distribution
"""

# Import the required libraries
import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %%

# Simulation parameters
tau1 = 0.5          # First inter-pulse time delay, μs
tau2 = 4.5          # Second inter-pulse time delay, μs
conc = 50           # Spin concentration, μs
lam1 = 0.40         # Amplitude of main contribution 
lam23 = 0.02        # Ampltiude of the "2+1" contributions
rmean = 4.0         # Mean distance, nm
rstd = 0.3          # Distance standard deviation, nm
Δr = 0.05           # Distance resolution, nm
rmin,rmax = 1.5,6   # Range of the distance vector, nm 
Δt = 0.008          # Time resolution, μs 
deadtime = 0.3      # Acquisition deadtime, μs

# Experimental time vector
t = np.arange(deadtime,tau1+tau2,Δt)
# Distance vector 
r = np.arange(rmin,rmax,Δr)

# Experiment model
experiment = dl.ex_4pdeer(tau1,tau2, pathways=[1,2,3])
reftime1 = experiment.reftimes[0]
reftime2 = experiment.reftimes[1]
reftime3 = experiment.reftimes[2]

# Construct the dipolar signal model
Vmodel = dl.dipolarmodel(t,r,Pmodel=dl.dd_gauss, experiment=experiment) 

# Function for the scaled background
Vinter_fcn = lambda lam1,lam23,conc: (1-lam1-2*lam23)*dl.bg_hom3d(t-reftime1,conc,lam1)*dl.bg_hom3d(t-reftime2,conc,lam23)*dl.bg_hom3d(t-reftime3,conc,lam23)

# Simulate the signal with orientation selection
Vsim = Vmodel(mean=rmean, std=rstd, conc=conc, scale=1, lam1=lam1, lam2=lam23, lam3=lam23, reftime1=reftime1, reftime2=reftime2, reftime3=reftime3)

# Plot the simulated signal
plt.figure(figsize=[4,3])
plt.plot(t,Vsim,'k',lw=2,label='V(t)')
plt.plot(t,Vinter_fcn(lam1,lam23,conc),'--',color='#f84862',lw=2,label='(1-λ)$V_{inter}$')
plt.legend()
plt.xlabel('Time (μs)')
plt.ylabel('V(t)')
plt.tight_layout()
plt.show()

# %%
