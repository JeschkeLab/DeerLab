#%% 
""" 
Validating multi-pathway models based on the data
-------------------------------------------------------------------

How to determine whether the dipolar signal model defined with the
set of specified dipolar pathways pathways is a proper descriptor of
the experimental data. 

This example shows how to use goodness-of-fit criteria to determine 
whether enough dipolar pathways have been accounted for in the model.  

A model that accurately describes the data must result in a residual vector
that is normally distributed, has zero mean, and has no significant 
autocorrelations. In this example, we will look at an experimental 4-pulse DEER
dataset acquired on a maltose-binding protein (MBP) and use the built-in 
goodness-of-fit tools to quantitatively validate whether the dataste is well described by a 
dipolar model with a single, two or three dipolar pathways.
""" 

import numpy as np 
import deerlab as dl 
import matplotlib.pyplot as plt 
#%% 

# File location
file = "../data/experimental_mbp_protein_4pdeer.DTA"

# Experiment information
t0 = 0.040
tau1 = 0.4 
tau2 = 3.0

# Laod and preprocess the data
t,Vexp = dl.deerload(file)
t = t[:-2]
Vexp = Vexp[:-2]
Vexp = dl.correctphase(Vexp) 
Vexp = Vexp/max(Vexp)
t = t- t[0] + t0

# Define the distance vector
r = np.arange(3,4.5,0.05)

# Loop over different dipolar models with varying number of pathways
for Npathways in [1,2,3]:
    print(f'Model with {Npathways} dipolar pathways:')

    # Construct the experiment model with different pathways
    experiment = dl.ex_4pdeer(tau1,tau2,pathways=np.arange(1,Npathways+1,1))

    # Construct the dipolar model with a non-parametric distance distribution 
    Vmodel = dl.dipolarmodel(t,r,experiment=experiment)

    # Define the compactness penalty for best results
    compactness = dl.dipolarpenalty(None,r,'compactness')

    # Fit the data to the current model
    results = dl.fit(Vmodel,Vexp,penalties=compactness)

    # Print the summary of the results
    print(results)

    # Plot the fit of the model to the data along its goodness-of-fit tests
    results.plot(axis=t, xlabel='t (μs)', gof=True)
    plt.suptitle(f'Model with {Npathways} dipolar pathways:')
    plt.show()

#%% [markdown]
# The first model is clearly underparametrized as it results in 
# non-normal residuals and strong correlations. This is supported by the large chi-squared value.
# Adding the second pathway seems to improve the description of the data, 
# as the residuals are now better distributed. However, there appears to be some autocorrelations 
# left and the chi-squared value still presents too large values. 
# Adding the third pathway results in the best description of the data, 
# with normally distributed residuals and no significant autocorrelations.

# %%
