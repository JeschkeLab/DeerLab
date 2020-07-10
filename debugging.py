
#%%
from deerlab import *
import numpy as np
import matplotlib.pyplot as plt
t = np.linspace(0,4,150) # us
r = np.linspace(1.5,6,500)

P = dd_wormgauss(r,[3.7, 10, 0.2])
plt.plot(r,P)
np.trapz(P,r)

# %%
