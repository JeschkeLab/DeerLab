# %% [markdown]
""" 
Manipulating DeerLab plots
-------------------------------------------------------------------

Basic example for manipulating the figures returned by DeerLab fit functions.
""" 

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl


# %% [markdown]
# All ``FitResult.plot()`` methods in DeerLab return a Matplotlib figure object that can be manipulated to change the default 
# plotting set in DeerLab. 
# In order to be able to change graphical elements in the plot, the figure object must be requested when calling the ``plot()``
# function. Also, since we want to change the graphics in the figure, we can disable the plotting provisionally until we have 
# edited the figure. This can be done via the ``show=False`` option.


#%% 

# Simulate some data
t = np.linspace(-0.1,4,250)        # time axis, µs
r = np.linspace(2,5,200)           # distance axis, nm
param = [3, 0.1, 0.2, 3.5, 0.1, 0.65, 3.8, 0.05, 0.15] # parameters for three-Gaussian model
P = dl.dd_gauss3(r,param)          # model distance distribution
lam = 0.5                          # modulation depth
B = dl.bg_hom3d(t,300,lam)         # background decay
K = dl.dipolarkernel(t,r,mod=lam,bg=B)    # kernel matrix
Vexp = K@P + dl.whitegaussnoise(t,0.01,seed=0)

# Run fit
fit = dl.fitmodel(Vexp,t,r,'P',dl.bg_hom3d,dl.ex_4pdeer)
fig = fit.plot(show=False)

# %% [markdown]
# If we call ``fig.get_children()`` we can see that there are three graphical elements: 
#   - A rectangle object (unimportant)
#   - Two subplot objects ``<AxesSubplot>`` one for each subplot
# We can get all of them separately from the method call.

#%%

# Get a list of the graphical elements inside the figure
_, plt1, plt2 = fig.get_children()

# %% 
# From these objects, every single graphical element can be manipulated or added to the existing plots. Here is an example of
# manipulation of the graphical settings:

# %% 

plt1.set_ylabel('Dipolar signal')               # Edit the y-axis of the top subplot
plt1.set_xlabel('t [μs]')                       # Edit the x-axis of the top subplot
plt1.grid(False)                                # Remove the grid lines of the top subplot

plt2.set_ylabel('Distance distribution [nm⁻¹]') # Edit the y-axis of the bottom subplot
plt2.set_xlabel('r [nm]')                       # Edit the x-axis of the bottom subplot
plt2.grid(False)                                # Remove the grid lines of the bottom subplot

plt1.get_children()[4].set_color('tab:red')     # Change the color of the data points
plt1.get_children()[5].set_color('k')        # Change the color of the fitted signal
plt1.get_children()[6].set_color('k')        # Change the color of the unmodulated contribution

# Change the color of the distribution and its confidence bands to red
plt2.get_children()[0].set_color('tab:red')     # Change the color of the fit
plt2.get_children()[1].set_color('tab:red')     # Change the color of the 50% CI
plt2.get_children()[2].set_color('tab:red')     # Change the color of the 95% CI

# Change the fontsize on both subplots
for ax in [plt1, plt2]:
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(13)

# Remove legends on both subplots
plt2.get_children()[12].remove()
plt1.get_children()[16].remove()

# Distribute everything in the available space to avoid overlapping
fig.set_tight_layout(True)


# %% [markdown]
# When all settings have been modified. The figure can be finally rendered by simply calling the ``display()`` function 
# with the modified figure object. This only works in IJupyter frameworks not as a .py script. 
