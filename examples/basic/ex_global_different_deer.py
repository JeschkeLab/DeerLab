# %%
""" 
Global fitting of multiple different DEER signals
-------------------------------------------------------------------------------

How to fit multiple signals from different DEER experiments to a model with a non-parametric
distribution and a homogeneous background, using Tikhonov regularization.
"""

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

#%%

# File location
path = "../data/"
file4p = "example_4pdeer_2.DTA"
file5p = "example_5pdeer_2.DTA"

# Experimental parameters (4pDEER)
tau1_4p = 0.5  # First inter-pulse delay, μs
tau2_4p = 3.5  # Second inter-pulse delay, μs
deadtime_4p = 0.1  # Acquisition deadtime, μs
# Experimental parameters (reversed 5pDEER)
tau1_5p = 2.9  # First inter-pulse delay, μs
tau2_5p = 3.3  # Second inter-pulse delay, μs
tau3_5p = 0.3  # Third inter-pulse delay, μs
deadtime_5p = 0.1  # Acquisition deadtime, μs

# Load the experimental data (4pDEER)
t4p, V4p = dl.deerload(path + file4p)
V4p = dl.correctphase(V4p)  # Phase correction
V4p = V4p / np.max(V4p)  # Rescaling (aesthetic)
t4p = t4p + deadtime_4p  # Account for deadtime
# Load the experimental data (reversed 5pDEER)
t5p, V5p = dl.deerload(path + file5p)
V5p = dl.correctphase(V5p)  # Phase correction
V5p = V5p / np.max(V5p)  # Rescaling (aesthetic)
t5p = t5p + deadtime_5p  # Account for deadtime

# Run fit
r = np.arange(2.5, 6, 0.05)

# Construct the individual dipolar signal models
V4pmodel = dl.dipolarmodel(
    t4p, r, experiment=dl.ex_4pdeer(tau1=tau1_4p, tau2=tau2_4p, pathways=[1, 2, 3])
)
V5pmodel = dl.dipolarmodel(
    t5p,
    r,
    experiment=dl.ex_rev5pdeer(
        tau1=tau1_5p, tau2=tau2_5p, tau3=tau3_5p, pathways=[1, 2]
    ),
)

# Make the joint model with the distribution as a global parameters
globalmodel = dl.merge(V4pmodel, V5pmodel, addweights=True)
globalmodel = dl.link(globalmodel, P=["P_1", "P_2"])

# Fit the model to the data (with fixed regularization parameter)
results = dl.fit(globalmodel, [V4p, V5p], weights=[1, 1])

# %%

plt.figure(figsize=[10, 7])
violet = "#4550e6"

# Extract fitted distance distribution
Pfit = results.P
Pci95 = results.PUncert.ci(95)
Pci50 = results.PUncert.ci(50)
for n, (t, V) in enumerate(zip([t4p, t5p], [V4p, V5p])):

    # Extract fitted dipolar signal
    Vfit = results.model[n]
    Vci = results.modelUncert[n].ci(95)

    plt.subplot(2, 2, 1 + n * 2)
    # Plot experimental data
    plt.plot(t, V, ".", color="grey", label="Data")
    # Plot the fitted signal
    plt.plot(t, Vfit, linewidth=3, color=violet, label="Fit")
    plt.fill_between(t, Vci[:, 0], Vci[:, 1], color=violet, alpha=0.3)
    plt.legend(frameon=False, loc="best")
    plt.xlabel("Time $t$ (μs)")
    plt.ylabel("$V(t)$ (arb.u.)")

# Plot the distance distribution
plt.subplot(1, 2, 2)
plt.plot(r, Pfit, linewidth=3, color=violet, label="Fit")
plt.fill_between(
    r,
    Pci95[:, 0],
    Pci95[:, 1],
    alpha=0.3,
    color=violet,
    label="95%-Conf. Inter.",
    linewidth=0,
)
plt.fill_between(
    r,
    Pci50[:, 0],
    Pci50[:, 1],
    alpha=0.5,
    color=violet,
    label="50%-Conf. Inter.",
    linewidth=0,
)
plt.legend(frameon=False, loc="best")
plt.autoscale(enable=True, axis="both", tight=True)
plt.xlabel("Distance $r$ (nm)")
plt.ylabel("$P(r)$ (nm$^{-1}$)")
plt.tight_layout()
plt.show()

# %%
