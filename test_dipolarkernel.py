
import numpy as np
import pandas as pd
from bg_models import *
from dipolarkernel import dipolarkernel,calckernelmatrix


def test_matrixsize_fresnel():
    #=======================================================================
    "Check non-square kernel matrix construction"
    Ntime = 50
    Ndist = 70
    t = np.linspace(-0.5,5,50)
    r = np.linspace(2,6,70)
    K = dipolarkernel(t,r)
    assert K.shape == (Ntime, Ndist)
    #=======================================================================


def test_singletime():
    #=======================================================================
    "Check that one can generate a kernel with one single time-domain point"

    t = 0.5
    r = np.linspace(1,5,100)
    K = dipolarkernel(t,r)

    assert np.shape(K)[0]==1
    #=======================================================================
    

def test_singledist():
    #=======================================================================
    "Check that one can generate a kernel with one single distance-domain point"

    t = np.linspace(0,5,50)
    r = 3.5
    K = dipolarkernel(t,r)

    assert np.shape(K)[1]==1
    #=======================================================================


def test_negative_time_fresnel():
    #=======================================================================
    "Check that kernel is constructed properly for negative times using the fresnel method"

    tneg = np.linspace(-5,0,50)
    tpos = np.linspace(0,5,50) 
    r = np.linspace(2,6,70)
    Kneg = dipolarkernel(tneg,r,method='fresnel')
    Kpos = dipolarkernel(tpos,r,method='fresnel')

    delta = abs(Kneg - np.flipud(Kpos))

    assert np.all(delta<1e-12)
    #=======================================================================


def test_negative_time_grid():
    #=======================================================================
    "Check that kernel is constructed properly for negative times using the grid method"

    tneg = np.linspace(-5,0,50)
    tpos = np.linspace(0,5,50) 
    r = np.linspace(2,6,70)
    Kneg = dipolarkernel(tneg,r,method='grid')
    Kpos = dipolarkernel(tpos,r,method='grid')

    delta = abs(Kneg - np.flipud(Kpos))

    assert np.all(delta<1e-12)
    #=======================================================================


def test_negative_time_integral():
    #=======================================================================
    "Check that kernel is constructed properly for negative times using the integral method"

    tneg = np.linspace(-5,0,50)
    tpos = np.linspace(0,5,50) 
    r = np.linspace(2,6,70)
    Kneg = dipolarkernel(tneg,r,method='integral')
    Kpos = dipolarkernel(tpos,r,method='integral')

    delta = abs(Kneg - np.flipud(Kpos))

    assert np.all(delta<1e-12)
    #=======================================================================


def test_value_fresnel():
    #=======================================================================
    "Test whether kernel matrix element (calculated using Fresnel integrals) is correct."

    # Generate kernel numerically
    t = 1 # us
    r = 1 # nm
    K = dipolarkernel(t,r,method='fresnel')

    # Kernel value for 1us and 1nm computed using Mathematica (FresnelC and FresnelS) 
    # and CODATA 2018 values for ge, muB, mu0, and h.
    Kref = 0.024697819895260188

    assert abs(K-Kref) < 1e-14
    #=======================================================================


def test_value_grid():
    #=======================================================================
    "Test whether kernel matrix element (calculated using the grid method) is correct."

    # Generate kernel numerically
    t = 1 # us
    r = 1 # nm
    K = dipolarkernel(t,r,method='grid')

    # Kernel value for 1us and 1nm computed using Mathematica (FresnelC and FresnelS) 
    # and CODATA 2018 values for ge, muB, mu0, and h.
    Kref = 0.024697819895260188

    assert abs(K-Kref) < 1e-8
    #=======================================================================


def test_value_integral():
    #=======================================================================
    "Test whether kernel matrix element (calculated using the integal method) is correct."

    # Generate kernel numerically
    t = 1 # us
    r = 1 # nm
    K = dipolarkernel(t,r,method='integral')

    # Kernel value for 1us and 1nm computed using Mathematica (FresnelC and FresnelS) 
    # and CODATA 2018 values for ge, muB, mu0, and h.
    Kref = 0.024697819895260188

    assert abs(K-Kref) < 1e-8
    #=======================================================================


def test_lambda():
    #=======================================================================
    "Check that dipolar kernel with modulation depth works"

    t = np.linspace(0,4,150) # us
    r = np.linspace(1,5,150) # nm
    dr = r[1]-r[0]
    lam = 0.4
    Kref = ((1-lam) + lam*dipolarkernel(t,r)/dr)*dr
    K = dipolarkernel(t,r,lam)

    assert np.all(abs(K-Kref) < 1e-14)
    #=======================================================================

def test_background():
    #=======================================================================
    "Check that dipolar kernel builds with the background included"
    t = np.linspace(0,3,80)
    r = np.linspace(2,6,100)
    B = bg_exp(t,0.5)
    dr = np.mean(np.diff(r))
    
    lam = 0.25
    KBref = (1 - lam + lam*dipolarkernel(t,r)/dr)*B[:,np.newaxis]*dr
    KB = dipolarkernel(t,r,lam,B)
    
    assert np.all(abs(KB - KBref) < 1e-5)
    #=======================================================================

def test_multipath():
    #=======================================================================
    "Check that multi-pathway kernels are properly generated"

    r = np.linspace(2,6,50)
    t1 = np.linspace(0,10,300)
    t2 = 0.3
    tau1 = 4.24
    tau2 = 4.92
    t = (tau1 + tau2) - (t1 + t2)
    prob = 0.8
    lam = np.array([1-prob, prob**2, prob*(1-prob)])
    T0 = [np.NaN, 0, tau2-t2]

    K = dipolarkernel(t,r,np.array([lam, T0]).T)

    unmodulated = pd.isnull(T0)
    Kref = sum(lam[unmodulated])
    Krenorm = Kref
    dr = np.mean(np.diff(r))
    for p in range(len(lam)):
        if not unmodulated[p]:
            Kref = Kref + lam[p]*calckernelmatrix(t-T0[p],r,'fresnel',[],[])
            Krenorm = Krenorm + lam[p]*calckernelmatrix(-T0[p],r,'fresnel',[],[])
    Kref = Kref/Krenorm
    Kref = Kref*dr

    assert np.all(abs(K-Kref) < 1e-3)
    #=======================================================================

def test_multipath_background():
    #=======================================================================
    "Check that multi-pathway kernels are properly generated with background"

    r = np.linspace(2,6,50)
    t1 = np.linspace(0,10,300)
    t2 = 0.3
    tau1 = 4.24
    tau2 = 4.92
    t = (tau1 + tau2) - (t1 + t2)
    prob = 0.8
    lam = np.array([1-prob, prob**2, prob*(1-prob)])
    T0 = [np.NaN, 0, tau2-t2]
    kappa = 0.3
    Bmodel = lambda t,lam: bg_exp(t,kappa,lam)

    # Reference
    unmodulated = pd.isnull(T0)
    Kref = sum(lam[unmodulated])
    Krenorm = Kref
    dr = np.mean(np.diff(r))
    for p in range(len(lam)):
        if not unmodulated[p]:
            Kref = Kref + lam[p]*calckernelmatrix(t-T0[p],r,'fresnel',[],[])
            Krenorm = Krenorm + lam[p]*calckernelmatrix(-T0[p],r,'fresnel',[],[])
    Kref = Kref/Krenorm
    Kref = Kref*dr
    
    Bref = 1
    Bnorm = 1
    for p in range(len(lam)):
        if not unmodulated[p]:
            Bref = Bref*Bmodel((t-T0[p]),lam[p])
            Bnorm = Bnorm*Bmodel(-T0[p],lam[p])
    Bref = Bref/Bnorm
    KBref = Kref*Bref[:,np.newaxis]

    # Output
    KB = dipolarkernel(t,r,np.array([lam, T0]).T,Bmodel)

    assert np.all(abs(KB - KBref) < 1e-3)
    #=======================================================================