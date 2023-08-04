# %%
"""
Simulating a single-pathway 4-pulse DEER signal
============================================================

An example on how to simulate a single-pathway 4-pulse DEER dipolar signal.  

Specifically, we simulate a single-pathway 4-pulse DEER dipolar signal arising from 
a Gaussian distance distribution
"""

# Import the required libraries
import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl
violet = '#4550e6'

# %%

# Simulation parameters
reftime = 0.5       # Refocusing time/Zero-time
tmax = 3            # Trace length, μs
mod = 0.40          # Modulation depth
conc = 50           # Spin concentration, μM
rmean = 3.0         # Mean distance, nm
rstd = 0.2          # Distance standard deviation, nm
Δr = 0.05           # Distance resolution, nm
rmin,rmax = 1.5,6   # Range of the distance vector, nm 
Δt = 0.008          # Time resolution, μs 
deadtime = 0.3      # Acquisition deadtime, μs
Vamp = 1            # Overall echo amplitude

# Experimental time vector
t = np.arange(deadtime,tmax,Δt)
# Distance vector 
r = np.arange(rmin,rmax,Δr)

# Construct the dipolar signal model
Vmodel = dl.dipolarmodel(t,r,Pmodel=dl.dd_gauss) 

# Simulate the signal with orientation selection
Vsim = Vmodel(mean=rmean, std=rstd, conc=conc, scale=Vamp, mod=mod, reftime=reftime)

# Plot the simulated signal
plt.figure(figsize=[4,3])
plt.plot(t,Vsim, color=violet,lw=2,label='V(t)')
plt.plot(t,Vamp*(1-mod)*dl.bg_hom3d(t,conc,mod),'--',color=violet,lw=2,label='(1-λ)$V_{inter}$')
plt.legend()
plt.xlabel('Time (μs)')
plt.ylabel('V(t)')
plt.tight_layout()
plt.show()

# %%
