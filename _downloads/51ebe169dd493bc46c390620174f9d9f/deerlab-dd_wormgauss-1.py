import deerlab as dl
import matplotlib.pyplot as plt 
import numpy as np 

model = dl.dd_wormgauss
r = np.linspace(1,6,400)
par0 = model.getmetadata()['par0']
P = model(r,*par0)
plt.figure(figsize=[6,3])
plt.plot(r,P,color='#4550e6')
plt.xlabel('r (nm)',fontsize=13)
plt.ylabel('P (nm⁻¹)',fontsize=13)
plt.grid(alpha=0.4)
plt.tick_params(labelsize=12)
plt.tick_params(labelsize=12)
plt.tight_layout()
plt.show()