# %% [markdown]
"""
Analyzing pseudo-titration (dose-respononse) curves with parameter-free distributions 
======================================================================================

How to use separable non-linear least squares (SNLLS)
to fit a pseudo-titration curve to multiple DEER datsets, using
parameter-free distance distributions.
"""

import numpy as np
import matplotlib.pyplot as plt
import deerlab as dl

# %% [markdown]
# Generating multiple datasets
#-----------------------------------------------------------------------------
# First, let's prepare the chemical desciption of the problem. In this example we will
# simulate a protein system in their states A (natural) and B (changed upon addition
# of a ligand L) given by the chemical equilibrium  A + L <-> B.

# %%
def chemicalequilibrium(Kdis,L):
    """Prepare equilibrium of type: A + L <-> B"""
    Ctot = 1 # total protein concentration, µM

    # # Get fraction of state B
    Kb = 1/Kdis
    xB = np.zeros_like(L)
    for q in range(len(L)):
        xB_ = np.roots([Kb*Ctot, -(Kb*L[q] + Kb*Ctot + 1), Kb*L[q]])
        try:
            xB[q] = xB_[(xB_<=1) & (xB_>=0)]
        except:
            xB[q] = np.minimum(1,np.maximum(0,xB_[0]))    
    # Get fraction of state A
    xA = 1 - xB

    return xA,xB

# %% [markdown]
# Next, we define the dipolar kernel model as the non-linear function of the
# SNLLS problem. This function needs to take the parameters and return a
# cell-array of kernels, each one for the corresponding datasets that we
# have. 
# Since we have a total distribution of the form 
#
#     ``P = xA*PA + xB*PB``
#
# we can define an augmented kernel as
#
#     ``K = [xA*KA xB*KB]``
#
# such that 
#
#     ``K@[PA PB] = V``
#
# and the vector ``[PA PB]`` constitutes the linear part fitted by SNLLS.

# %%
def Kmodel(par,ts,rA,rB,L):

    Nsignals = len(ts)

    # Unpack parameters
    lam,k,Kdis = par

    # Get fractions for given KD
    [xA,xB] = chemicalequilibrium(Kdis,L)

    Ks = [[]]*Nsignals
    # General the dipolar kernels
    for i in range(Nsignals):
        B = dl.bg_exp(ts[i],k)
        # Kernel for fraction A
        KstateA = dl.dipolarkernel(ts[i],rA,mod=lam,bg=B)
        # Kernel for fraction B
        KstateB = dl.dipolarkernel(ts[i],rB,mod=lam,bg=B)
        Ks[i] = np.concatenate((xA[i]*KstateA, xB[i]*KstateB),axis=1)

    return Ks

# %% [markdown]
# Now, we can simulate multiple signals corresponding to different concentrations
# of added ligand. 

# %%
# Time axes
ts = [[]]*7
ts[0] = np.linspace(-0.2,3,100)
ts[1] = np.linspace(-0.1,5,300)
ts[2] = np.linspace(-0.5,2,200)
ts[3] = np.linspace(-0.1,1,100)
ts[4] = np.linspace(-0.2,6,300)
ts[5] = np.linspace(-0.2,3,300)
ts[6] = np.linspace(-0.1,4,100)
Nsignals = len(ts)

# Distance axes for states A and B
rA = np.linspace(1,8,100)
rB = np.linspace(1,8,100)

# Distributions for states A and B
PstateA = dl.dd_gauss(rA,[5.5, 0.25])
PstateB = dl.dd_gauss2(rB,[4.5, 0.4, 0.4, 3.5, 0.35, 0.6])

L = [0.3, 1, 3, 10, 30, 100, 300] # total ligand concentration, µM
Kdis = 5.65  # dissociation constant, µM

# Populations of states A and B
[xA,xB] = chemicalequilibrium(Kdis,L)

# Global kernel model
Ks = Kmodel([0.25, 0.1, Kdis],ts,rA,rB,L)

# Simulate dipolar signals
Vs = [[]]*Nsignals
for i in range(Nsignals):
    Vs[i] = Ks[i]@np.concatenate((PstateA, PstateB)) + dl.whitegaussnoise(ts[i],0.01,seed=i)

# %% [markdown]
# Psuedotitration SNLLS Analysis
#-----------------------------------------------------------------------------
# For simplification, we will assume that all DEER traces have the same
# background function and modulation depth. Thus, we will fit the
# modulations depth (lam) and background decay constant (k) globally along
# the dissociation constant (KD).

# %%
# Non-linear parameters:
#       lam  k   KD
par0 = [0.5, 0.5,  5]  # start values 
lb   = [ 0,   0,   1]  # lower bounds
ub   = [ 1,   1,  10] # upper bounds

# Linear parameters:
#     |-------PA--------||--------PB--------|
lbl = np.concatenate((np.zeros_like(rA), np.zeros_like(rB))) # Non-negativity constraint
ubl = [] # Unconstrained

# Run SNLLS optimization
fit = dl.snlls(Vs,lambda p: Kmodel(p,ts,rA,rB,L),par0,lb,ub,lbl,ubl)
# Extract fit results
parfit = fit.nonlin
Pfit = fit.lin

# Extract the fitted disociation constant value and its 95#-confidence interval
Kdisfit = parfit[2]
parci = fit.nonlinUncert.ci(95)
KDci = parci[2,:]

# Print result
print(f'Kdis = {Kdisfit:.2f}({KDci[0]:.2f}-{KDci[1]:.2f})µM')

# %%
# Plot results
plt.figure(figsize=(12,12))

# Simulate fits
Ksfit = Kmodel(parfit,ts,rA,rB,L)
Vsfit = []
plt.subplot(3,2,(1,3))
for i in range(Nsignals):
    Vci = fit.modelUncert[i].ci(95)
    Vfit = fit.model[i]
    plt.plot(ts[i],Vs[i]+i/3,'.',color='grey',alpha=0.7,linestyle=None)
    plt.plot(ts[i],Vfit+i/3,'b',linewidth=1.5)
    plt.fill_between(ts[i],Vci[:,0]+i/3,Vci[:,1]+i/3,color='b',alpha=0.2)
plt.grid(alpha =0.3)
plt.xlabel('t (µs)')
plt.ylabel('V (arb.u.)')
plt.legend(['data','fit'])

xAfit,xBfit = chemicalequilibrium(Kdisfit,L)

xA_model = lambda param: chemicalequilibrium(param[2],L)[0]
xB_model = lambda param: chemicalequilibrium(param[2],L)[1]
xA_uq = fit.nonlinUncert.propagate(xA_model)
xB_uq = fit.nonlinUncert.propagate(xB_model)
xA_ci = xA_uq.ci(95)
xB_ci = xB_uq.ci(95)

plt.subplot(2,2,(2,4))
for i in range(Nsignals):
    PAfit = xAfit[i]*Pfit[0:len(rA)]
    PBfit = xBfit[i]*Pfit[len(rA):len(rB)+len(rA)]
    PAci = xAfit[i]*fit.linUncert.ci(95)[0:len(rA)]
    PBci = xBfit[i]*fit.linUncert.ci(95)[len(rA):len(rB)+len(rA)]

    plt.plot(rA,PAfit+1.2*i,'tab:red',rB,PBfit+1.2*i,'tab:blue',linewidth=1.5)
    plt.fill_between(rA,PAci[:,0]+1.2*i,PAci[:,1]+1.2*i,color='tab:red',alpha=0.2)
    plt.fill_between(rB,PBci[:,0]+1.2*i,PBci[:,1]+1.2*i,color='tab:blue',alpha=0.2)

plt.grid(alpha =0.3)
plt.xlabel('r (nm)')
plt.ylabel('P (nm⁻¹)')
plt.legend(['state A','state B'])
plt.xlim([2,7])

plt.subplot(325)
plt.plot(np.log10(L),xA,'tab:red',np.log10(L),xB,'tab:blue')
plt.plot(np.log10(L),xAfit,'o',color='tab:red')
plt.plot(np.log10(L),xBfit,'o',color='tab:blue')
plt.fill_between(np.log10(L),xA_ci[:,0],xA_ci[:,1],color='tab:red',alpha=0.2)
plt.fill_between(np.log10(L),xB_ci[:,0],xB_ci[:,1],color='tab:blue',alpha=0.2)
plt.grid(alpha =0.3)
plt.xlabel('log$_{10}$([L])')
plt.ylabel('Fractions')
plt.legend(['state A','state B'])
plt.ylim([0,1])

plt.show()
# %%
