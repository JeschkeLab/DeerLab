# %% [markdown]
""" 
Plotting the dipolar spectrum / Pake pattern
-------------------------------------------------------------------------

The dipolar spectrum or Pake patern is the Fourier transform of the dipolar signal, and can be useful for visualising the dipolar frequencies present in the signal. We can also overlay the dipolar spectrum of the fitted model to see how well the fit captures the dipolar frequencies in the data.
""" 


import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl
# %%
# File location
path = '../data/'
file = 'example_4pdeer_1.DTA'

# Experimental parameters
tau1 = 0.3      # First inter-pulse delay, μs
tau2 = 4.0      # Second inter-pulse delay, μs
tmin = 0.1      # Start time, μs

# Load the experimental data
t,Vexp = dl.deerload(path + file)

# Pre-processing
Vexp = dl.correctphase(Vexp) # Phase correction
Vexp = Vexp/np.max(Vexp)     # Rescaling (aesthetic)
t = t - t[0]                  # Account for zerotime
t = t + tmin    
# Distance vector
r = np.arange(2.5,5,0.01) # nm

# Construct the model
experimentInfo = dl.ex_4pdeer(tau1,tau2, pathways=[1])
Vmodel = dl.dipolarmodel(t,r, experiment = experimentInfo, Bmodel=dl.bg_hom3d)

# Fit the model to the data
results = dl.fit(Vmodel,Vexp)

# Print results summary
print(results)
results.plot()

# %% [markdown]
""" 
Fourier transform
+++++++++++++++++
""" 
import numpy.fft as fft


# Fourier transform of the data and model
# Before background division
Vexp_ft = fft.fftshift(fft.fft(Vexp))
Vmodel_ft = fft.fftshift(fft.fft(results.model))
ft_axis = fft.fftshift(fft.fftfreq(len(t), d=t[1]-t[0]))

# After background division
Vexp_div = Vexp / results.bg -1
Vmodel_div = results.model / results.bg - 1
Vexp_div_ft = fft.fftshift(fft.fft(Vexp_div))
Vmodel_div_ft = fft.fftshift(fft.fft(Vmodel_div))

# Plotting
violet = '#4550e6'
fig,axs = plt.subplots(1,2,figsize=[10,6],sharex=True)
axs[0].set_title('Before background division')
axs[0].plot(ft_axis, np.abs(Vexp_ft), '.', color='grey', label='Data')
axs[0].plot(ft_axis, np.abs(Vmodel_ft), linewidth=3, color=violet, label='Fit')
axs[0].legend()
axs[0].set_ylabel('Magnitude (arb.u.)')
axs[0].set_xlabel('Dipolar frequency (MHz)')
axs[0].set_xlim(-10, 10) # Limit x-axis to focus on relevant frequencies


axs[1].set_title('After background division')
axs[1].plot(ft_axis, np.abs(Vexp_div_ft), '.', color='grey', label='Data')
axs[1].plot(ft_axis, np.abs(Vmodel_div_ft), linewidth=3, color=violet, label='Fit')
axs[1].legend()
axs[1].set_xlabel('Dipolar frequency (MHz)')
axs[1].set_ylabel('Magnitude (arb.u.)')
axs[1].set_xlim(-10, 10) # Limit x-axis to

# %%
