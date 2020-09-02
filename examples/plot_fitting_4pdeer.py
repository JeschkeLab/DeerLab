# %% [markdown]
""" 
Basic fitting of a 4-pulse DEER signal, non-parametric distribution
-------------------------------------------------------------------

How to fit a simple 4-pulse DEER signal with a parameter-free distribution.
""" 

import numpy as np
import matplotlib.pyplot as plt
from deerlab import *

# %% [markdown]
#Generate data
#--------------

t = np.linspace(-0.1,4,250)      # time axis, us
r = np.linspace(1.5,6,len(t))    # distance axis, ns
param = [3, 0.1, 0.2, 3.5, 0.1, 0.65, 3.8, 0.05, 0.15] # parameters for three-Gaussian model
P = dd_gauss3(r,param)          # model distance distribution
lam = 0.5                        # modulation depth
B = bg_hom3d(t,300,lam)          # background decay
K = dipolarkernel(t,r,lam,B)
np.random.seed(0)
Vexp = K@P + whitegaussnoise(t,0.01)


# %% [markdown]
# Run fit
#---------
fit = fitsignal(Vexp,t,r,'P',bg_hom3d,ex_4pdeer)
fit.plot()

