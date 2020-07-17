import numpy as np
import deerlab as dl
from deerlab.snlls import snlls

def deerlab_default(t,V):
    def Kmodel(p,t,r):

        # Unpack parameters
        lam = p[0]
        k = p[1]

        # Get background
        B = dl.bg_exp(t,k)

        # Generate 4pDEER kernel
        K = dl.dipolarkernel(t,r,lam,B)
        return K


    t = np.squeeze(np.asarray(t))
    V = np.squeeze(np.asarray(V))

    t = t/1000

    # phase correction
    V = np.real(V)
    # zerotime correction
    t = t - t[np.argmax(V)]

    r = np.linspace(2,6,100)

    #--------------------------
    # Non-linear parameters:
    #--------------------------
    #       lam  c0
    #--------------------------
    par0 = [0.5, 0.5 ] # Start values
    lb   = [ 0 , 0.05] # lower bounds
    ub   = [ 1 ,  1  ] # upper bounds

    #--------------------------
    # Linear parameters: 
    #--------------------------
    #          Pfit
    #--------------------------
    lbl = np.zeros(len(r)) # Non-negativity constraint of P
    ubl = []

    Amodel = lambda p: Kmodel(p,t,r)

    # Run SNLLS optimization
    parfit, Pfit, uq = snlls(V,Amodel,par0,lb,ub,lbl,ubl)

    # Get fitted model
    Vfit = Kmodel(parfit,t,r)@Pfit

    V = V/np.max(Vfit)
    Vfit = Vfit/np.max(Vfit)
    dr = np.mean(np.diff(r))
    Pnorm = np.sum(Pfit)/dr
    Pfit = Pfit/Pnorm
    P95 = uq.ci(95,'lin')/Pnorm
    P50 = uq.ci(50,'lin')/Pnorm
    parci = uq.ci(95,'nonlin')

    return V,Vfit,Pfit,r,P95,P50,parfit,parci