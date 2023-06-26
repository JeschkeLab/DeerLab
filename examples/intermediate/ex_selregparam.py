# %% [markdown]
""" 
Analysing the selection of regularisation parameter
-------------------------------------------------------------------------

This example demonstates how to generate plots from the regularisation parameter selection,
and how to use this infomation.

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
deadtime = 0.1  # Acquisition deadtime, μs

# Load the experimental data
t,Vexp = dl.deerload(path + file)

# Pre-processing
Vexp = dl.correctphase(Vexp) # Phase correction
Vexp = Vexp/np.max(Vexp)     # Rescaling (aesthetic)
t = t + deadtime             # Account for deadtime

# Distance vector
r = np.linspace(1.5,7,50) # nm

# Construct the model
Vmodel = dl.dipolarmodel(t,r, experiment=dl.ex_4pdeer(tau1,tau2, pathways=[1]))

# Fit the model to the data with compactness criterion
results= dl.fit(Vmodel,Vexp)
print(results)

# %%

fig, axs =plt.subplots(1,3, figsize=(9,4),width_ratios=(1,1,0.1))
fig.tight_layout()
alphas = results.regparam_stats['alphas_evaled'][1:]
funcs = results.regparam_stats['functional'][1:]


idx = np.argsort(alphas)

axs[0].semilogx(alphas[idx], funcs[idx],marker='.',ls='-')
axs[0].set_title(r"$\alpha$ selection functional");
axs[0].set_xlabel("Regularisation Parameter")
axs[0].set_ylabel("Functional Value ")





# Just the final L-Curve

cmap = plt.get_cmap('plasma')
import matplotlib as mpl

x = results.regparam_stats['residuals']
y = results.regparam_stats['penalties']
idx = np.argsort(x)


axs[1].loglog(x[idx],y[idx])

n_points = results.regparam_stats['alphas_evaled'].shape[-1]
lams = results.regparam_stats['alphas_evaled']
norm = mpl.colors.LogNorm(vmin=lams[1:].min(), vmax=lams.max())
for i in range(n_points):
    # plt.annotate(fr"$\lambda =$ {lams[i]:.2g}", (x[i], y[i] + 0.2))
    axs[1].plot(x[i], y[i],marker = '.', ms=8, color=cmap(norm(lams[i])))

i_optimal = np.argmin(np.abs(lams - results.regparam))
axs[1].annotate(fr"$\lambda =$ {results.regparam:.2g}", xy = (x[i_optimal],y[i_optimal]),arrowprops=dict(facecolor='black', shrink=0.05, width=5), xytext=(20, 20),textcoords='offset pixels')
fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),cax=axs[2])
axs[1].set_ylabel("Penalties")
axs[2].set_ylabel("Regularisation Parameter")
axs[1].set_xlabel("Residuals")
axs[1].set_title("L-Curve");

# %%
