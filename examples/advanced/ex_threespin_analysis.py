#%%
"""
Analyzing 4-pulse DEER data acquired on three-spin systems
============================================================================

As in the publication referenced below, this example will take two 4-pulse DEER signals acquired
on the same protein sample conatining three nitroxide spins with different attenuation levels of the pump 
pulse power.

For the original model and more information on these systems please refer to: 
L. Fábregas Ibáñez, M. H. Tessmer, G. Jeschke, and S. Stoll. 
Dipolar pathways in multi-spin and multi-dimensional dipolar EPR spectroscopy
Phys. Chem. Chem. Phys., 24 2022, 22645-22660
""" 
#%%
import numpy as  np 
import deerlab as dl 

# Load experimental data
#files = [f'../data/triradical_protein_deer_{dB}dB.DTA' for dB in [0,6,9]]
files = [f'D:\lufa\projects\DeerLab\DeerLab\\examples\data\\triradical_protein_deer_{dB}dB.DTA' for dB in [0,6,9]]
# Experiment information
t0  = 0.280 # Acquisition deadtime, μs
tau1 = 0.40 # First interpulse delay, μs
tau2 = 9.00 # Second interpulse delay, μs

# Construct 4-pulse DEER experiment model
my4pdeer = dl.ex_4pdeer(tau1,tau2,pathways=[1])

# Loop over the different datasets 
Vmodels,Vexps,ts,Vexps_sub,ts_sub = [],[],[],[],[]
for n,file in enumerate(files):
    # Load the dataset
    t,Vexp, descriptor = dl.deerload(file,full_output=True)
    t = t[:-80]
    Vexp = Vexp[:-80]
    # Adjust the start time
    t = t - t[0] + t0

    # Pre-processing
    Vexp = dl.correctphase(Vexp)
    Vexp /= np.max(Vexp) 

    # Store the pre-processed datasets in a list
    Vexps.append(Vexp)
    ts.append(t) 

    # Subsampling 
    # (required for efficient analysis in densely sampled datasets)
    sampling = np.arange(0,len(t),4) # Take every 4th point
    t_sub = t[sampling]
    Vexp_sub = Vexp[sampling]

    # Store the subsampled datasets in a list
    Vexps_sub.append(Vexp_sub)
    ts_sub.append(t_sub) 

    # Construct the three-spin dipolar model
    Vmodel = dl.dipolarmodel(t_sub,spins=3,experiment=my4pdeer, minamp=0.01)

    # Add dipolar model to list of models
    Vmodels.append(Vmodel)

# Construct a global dipolar model describing all datasets
Vglobal = dl.merge(*Vmodels)
Vglobal = dl.link(Vglobal,
        rmean1=[f'rmean1_{n+1}' for n in range(len(Vmodels))],
        rmean2=[f'rmean2_{n+1}' for n in range(len(Vmodels))],
        rmean3=[f'rmean3_{n+1}' for n in range(len(Vmodels))],
        chol11=[f'chol11_{n+1}' for n in range(len(Vmodels))],
        chol22=[f'chol22_{n+1}' for n in range(len(Vmodels))],    
        chol33=[f'chol33_{n+1}' for n in range(len(Vmodels))],    
        chol21=[f'chol21_{n+1}' for n in range(len(Vmodels))],    
        chol31=[f'chol31_{n+1}' for n in range(len(Vmodels))],    
        chol32=[f'chol32_{n+1}' for n in range(len(Vmodels))],    
        conc=[f'conc_{n+1}' for n in range(len(Vmodels))],    
        reftime1=[f'reftime1_{n+1}' for n in range(len(Vmodels))],  
)
# Freeze the Cholesky-factors accounting for the correlation coefficients
# to zero (not always applicable) 
Vglobal.chol21.freeze(0)
Vglobal.chol31.freeze(0)
Vglobal.chol32.freeze(0)

# Fit the model to the data 
results = dl.fit(Vglobal, Vexps_sub, reg=False, ftol=1e-5) 

# %%
# Plot the fitted datasets 
results.plot(axis=ts_sub,xlabel='time (μs)')

# Print the fit summary 
print(results)

# %%
