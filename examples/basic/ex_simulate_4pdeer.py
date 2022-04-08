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
reftime = 0         # Refocusing time, μs
conc = 50           # Spin concentration, μs
moddepth = 0.4      # Modulation depth
rmean = 4.0         # Mean distance, nm
rstd = 0.3          # Distance standard deviation, nm
rmin,rmax = 1.5,6   # Range of the distance axis, nm
Δr = 0.05           # Distance resolution, nm
tmin,tmax = -0.5,5  # Range of the time trace, μs 
Δt = 0.008          # Time resolution, μs 

# Experimental time vector
t = np.arange(tmin,tmax,Δt)
# Distance vector 
r = np.arange(rmin,rmax,Δr)

# Construct the dipolar signal model
Vmodel = dl.dipolarmodel(t,r,Pmodel=dl.dd_gauss) 
# Function for the scaled background
Vinter_fcn = lambda mod,conc: (1-mod)*dl.bg_hom3d(t,conc,mod)

# Simulate the signal with orientation selection
Vsim = Vmodel(mean=rmean, std=rstd, reftime=reftime, mod=moddepth, conc=conc, scale=1)

# Plot the simulated signal
plt.figure(figsize=[4,3])
plt.plot(t,Vsim,'k',lw=2,label='V(t)')
plt.plot(t,Vinter_fcn(moddepth,conc),'--',color='#f84862',lw=2,label='(1-λ)$V_{inter}$')
plt.legend()
plt.xlabel('Time (μs)')
plt.ylabel('V(t)')
plt.show()

# %%
