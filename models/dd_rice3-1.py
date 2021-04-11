import deerlab as dl
import matplotlib.pyplot as plt 
import numpy as np 
model = dl.dd_rice3
r = np.linspace(2,5,400)
info = model() 
par0 = info['Start']
P = model(r,par0)
plt.figure(figsize=[6,3])
plt.plot(r,P)
plt.xlabel('r (nm)',fontsize=13)
plt.ylabel('P (nm⁻¹)',fontsize=13)
plt.grid(alpha=0.4)
plt.tick_params(labelsize=12)
plt.tick_params(labelsize=12)
plt.tight_layout()