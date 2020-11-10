import deerlab as dl
import matplotlib.pyplot as plt
import numpy as np
model = dl.dd_wormchain
r = np.linspace(0.5,8,400)
info = model()
par0 = info['Start']
P = model(r,par0)
plt.figure(figsize=[6,3])
plt.plot(r,P)
plt.xlabel('r [nm]',fontsize=13)
plt.ylabel('P(r) [nm$^{-1}$]',fontsize=13)
plt.grid(alpha=0.4)
plt.tick_params(labelsize=12)
plt.tick_params(labelsize=12)
plt.tight_layout()