# %% 
"""
Analyzing pseudo-titration (dose-response) curves with non-parametric distributions 
======================================================================================

How to fit a pseudo-titration curve to multiple DEER datsets, using
non-parametric distance distributions.

In this example we will simulate a protein system in their states 
A (natural) and B (changed upon addition of a ligand L) given by
the chemical equilibrium  A + L <-> B.
"""

import deerlab as dl 
import matplotlib.pyplot as plt 
import numpy as np 


def chemicalequilibrium(Kdis,L):
    """Prepare equilibrium of type: A + L <-> B"""
    Ctot = 1 # total protein concentration, µM
    Kb = 1/Kdis
    xB = np.roots(np.squeeze([Kb, -(Kb*L + Kb*Ctot + 1), Kb*L]))
    xB = xB[(xB<=1) & (xB>=0)]
    return xB

t1,V1 = np.load('../data/example_data_titration_#1.npy')
t2,V2 = np.load('../data/example_data_titration_#2.npy')
t3,V3 = np.load('../data/example_data_titration_#3.npy')
t4,V4 = np.load('../data/example_data_titration_#4.npy')
t5,V5 = np.load('../data/example_data_titration_#5.npy')

ts = [t1,t2,t3,t4,t5]
Vs = [V1,V2,V3,V4,V5]

# Total ligand concentrations used in experiments
L = [0.3, 3, 10, 30, 300] # µM

# Distance vector
r = np.linspace(1.5,6,90)

# Construct a non-parametric distance distribution that is a
# linear combination of two non-parametric distributions
PAmodel = dl.freedist(r)
PBmodel = dl.freedist(r)
Pmodel = dl.lincombine(PAmodel,PBmodel,addweights=True)

# Construct the dipolar models of the individual signals
Vmodels = [dl.dipolarmodel(t,r,Pmodel) for t in ts]

# Create the global model
titrmodel = dl.merge(*Vmodels)
# Make the two components of the distance distriution global
titrmodel = dl.link(titrmodel, 
                PA = ['P_1_1', 'P_1_2', 'P_1_3', 'P_1_4', 'P_1_5'],
                PB = ['P_2_1', 'P_2_2', 'P_2_3', 'P_2_4', 'P_2_5'])

# Functionalize the chemical equilibrium model
titrmodel.addnonlinear('Kdis',lb=3,ub=7,par0=5,description='Dissociation constant')

titrmodel = dl.relate(titrmodel, 
            weight_2_1 = lambda weight_1_1: 1-weight_1_1, weight_1_1 = lambda Kdis: chemicalequilibrium(Kdis,L[0]),
            weight_2_2 = lambda weight_1_2: 1-weight_1_2, weight_1_2 = lambda Kdis: chemicalequilibrium(Kdis,L[1]),
            weight_2_3 = lambda weight_1_3: 1-weight_1_3, weight_1_3 = lambda Kdis: chemicalequilibrium(Kdis,L[2]),
            weight_2_4 = lambda weight_1_4: 1-weight_1_4, weight_1_4 = lambda Kdis: chemicalequilibrium(Kdis,L[3]),
            weight_2_5 = lambda weight_1_5: 1-weight_1_5, weight_1_5 = lambda Kdis: chemicalequilibrium(Kdis,L[4]))
            
# Fit the model to the data
fit = dl.fit(titrmodel,Vs,regparam = 0.5)

# %%

# Evaluate the dose-response curve at the fit with confidence bands
xAfcn = lambda Kdis: np.squeeze(np.array([chemicalequilibrium(Kdis,Ln) for Ln in L]))
xBfcn = lambda Kdis: np.squeeze(np.array([1 - chemicalequilibrium(Kdis,Ln) for Ln in L]))
xAfit = xAfcn(fit.Kdis)
xBfit = xBfcn(fit.Kdis)
xAci = fit.propagate(xAfcn,lb=np.zeros_like(L),ub=np.ones_like(L)).ci(95)
xBci = fit.propagate(xBfcn,lb=np.zeros_like(L),ub=np.ones_like(L)).ci(95)

# Plot the dose-reponse curve
plt.plot(L,xAfit,'-o')
plt.fill_between(L,xAci[:,0],xAci[:,1],alpha=0.5)
plt.plot(L,xBfit,'-o')
plt.fill_between(L,xBci[:,0],xBci[:,1],alpha=0.5)
plt.xscale('log')
plt.xlabel('Ligand concentration (μM)')
plt.ylabel('Molar fraction')
plt.legend(['State A (natural)','State B (ligand)'],frameon=False,loc='best')
plt.title(r'$K_\mathrm{dis}$'+f' = {fit.Kdis:.2f} ({fit.KdisUncert.ci(95)[0]:.2f}-{fit.KdisUncert.ci(95)[1]:.2f})'+' µM$^{-1}$')
plt.autoscale(enable=True, axis='both', tight=True)
plt.show() 

# Plot the fitted signals and distance distributions
plt.figure(figsize=[10,10])

plt.subplot(121)
for n,(t,Vexp,Vfit) in enumerate(zip(ts,Vs,fit.model)):
    plt.plot(t,n/2 + Vexp,'.',color='grey')
    plt.plot(t,n/2 + Vfit,'k',linewidth=2)
plt.legend(['Data','Fit'],frameon=False,loc='best')
plt.xlabel('Time $t$ (μs)')
plt.ylabel('$V(t)$ (arb.u.)')

plt.subplot(122)
for n,(xA,xB) in enumerate(zip(xAfit,xBfit)): 

    Pfit = Pmodel(P_1=fit.PA,P_2=fit.PB,weight_1=xA,weight_2=xB)
    Pfit /= np.trapz(Pfit,r)
    if n>1: label=None
    plt.plot(r,2*n + Pfit,'k',label='Total contribution' if n<1 else None)
    plt.fill(r,2*n + xA*fit.PA,color='tab:blue',alpha=0.5,label='State A (natural)' if n<1 else None)
    plt.fill(r,2*n + xB*fit.PB,color='tab:orange',alpha=0.5,label='State B (ligand)' if n<1 else None)

plt.legend(frameon=False,loc='best')
plt.ylabel('$P(r)$')
plt.xlabel('Distance $r$ (nm)')
plt.show() 


# %%