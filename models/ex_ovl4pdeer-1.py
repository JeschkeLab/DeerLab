import deerlab as dl
import matplotlib.pyplot as plt
import numpy as np
model = dl.ex_ovl4pdeer
t = np.linspace(-0.5,5,400)
r = np.linspace(2,5,200)
info = dl.dd_gauss()
par0 = info['Start']
P = dl.dd_gauss(r,par0)
info = model()
par0 = info['Start']
paths = model(par0)
K = dl.dipolarkernel(t,r,paths,lambda t,lam: dl.bg_hom3d(t,80,lam))
V = K@P
plt.figure(figsize=[6,3])
plt.plot(t,V)
plt.xlabel('t [μs]',fontsize=13)
plt.ylabel('V(t)',fontsize=13)
plt.grid(alpha=0.4)
plt.tick_params(labelsize=12)
plt.tight_layout()