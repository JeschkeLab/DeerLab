#%% [markdown]
""" 
Saving and loading DeerLab fit results
-------------------------------------------------------------------------

DeerLab FitResult objects can be saved and loaded to file, with all relevant information (fitted model, uncertainties, regularization parameter, etc.) preserved.
This is useful for archiving results, sharing them with collaborators, or importing fits from DeerAnalysis 2026 into Python for plotting.

""" 


# %%
# First we will fit a basic model to generate a result and then save it to file.

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl
path = '../data/'
file = 'example_4pdeer_1.DTA'

tau1 = 0.3; tau2 = 4.0; tmin = 0.1; 

t,Vexp = dl.deerload(path + file)

# Pre-processing
Vexp = dl.correctphase(Vexp) # Phase correction
Vexp = Vexp/np.max(Vexp)     # Rescaling (aesthetic)
t = t - t[0]                  # Account for zerotime
t = t + tmin    
# Distance vector
r = np.arange(2.5,5,0.01) # nm

# Construct the model
Vmodel = dl.dipolarmodel(t,r, experiment = dl.ex_4pdeer(tau1,tau2, pathways=[1]))

# Fit the model to the data
results = dl.fit(Vmodel,Vexp)


# %%
# We will then save the results to an HDF5 file.
# The file can also be saved as a JSON or TOML file, if needed. A HDF5 file is easy to read using MATLAB, Origin or similar software. 

dl.save('fit_result.hdf5', results)

#%%
# The file can then be loaded back into Python, with all information preserved.
results2 = dl.load('fit_result.hdf5')

# We can check that the loaded results are the same as the original results by comparison
print('Model result is the same:', np.allclose(results.model, results2.model))

# Built in methods of the FitResult object can be used as normal, e.g. plotting the results
results2.plot();

# %%
