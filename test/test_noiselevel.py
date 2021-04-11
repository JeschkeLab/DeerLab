import numpy as np
from deerlab import noiselevel, whitegaussnoise, dipolarkernel
from deerlab.dd_models import dd_gauss
from deerlab.bg_models import bg_exp


def test_filtered_movmean():
#============================================================
    "Check estimation of noiselevel using a moving-mean filter"

    np.random.seed(1)
    t = np.linspace(0,3,200)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[3, 0.5])
    lam = 0.25
    B = bg_exp(t,1.5)
    noise = whitegaussnoise(t,0.05)
    V = dipolarkernel(t,r,mod=lam,bg=B)@P + noise


    truelevel = np.std(noise)
    approxlevel = noiselevel(V,'movmean',3)

    assert abs(approxlevel - truelevel) < 1e-2
#============================================================


def test_der():
#============================================================
    "Check estimation of noiselevel using a moving-mean filter"

    np.random.seed(1)
    t = np.linspace(0,3,200)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[3, 0.5])
    lam = 0.25
    B = bg_exp(t,1.5)
    noise = whitegaussnoise(t,0.05)
    V = dipolarkernel(t,r,mod=lam,bg=B)@P + noise


    truelevel = np.std(noise)
    approxlevel = noiselevel(V,'der')

    assert abs(approxlevel - truelevel) < 1e-2
#============================================================



def test_reference():
#============================================================
    "Check estimation of noiselevel using a reference signal"

    np.random.seed(1)
    t = np.linspace(0,3,200)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[3, 0.5])
    lam = 0.25
    B = bg_exp(t,1.5)
    Vref = dipolarkernel(t,r,mod=lam,bg=B)@P
    noise = whitegaussnoise(t,0.03)
    V = Vref + noise


    truelevel = np.std(noise)
    approxlevel = noiselevel(V,'reference',Vref)

    assert abs(approxlevel - truelevel) < 1e-2
#============================================================


def test_filtered_savgol():
#============================================================
    "Check estimation of noiselevel using a Savitzky-Golay filter"

    np.random.seed(1)
    t = np.linspace(0,3,200)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[3, 0.5])
    lam = 0.25
    B = bg_exp(t,1.5)
    noise = whitegaussnoise(t,0.03)
    V = dipolarkernel(t,r,mod=lam,bg=B)@P + noise

    truelevel = np.std(noise)
    approxlevel = noiselevel(V,'savgol',11,3)

    assert abs(approxlevel - truelevel) < 1e-2
#============================================================



def test_multiscan():
#============================================================
    "Check estimation of noiselevel using multiple scans of a signal"

    np.random.seed(1)
    t = np.linspace(0,5,300)
    r = np.linspace(2,6,200)
    P = dd_gauss(r,[4, 0.4])
    K = dipolarkernel(t,r)

    sigma_ref = 0.1
    N = 500
    V = np.zeros((len(t),N))
    for i in range(N):
        V[:,i] = K@P + whitegaussnoise(t,sigma_ref)

    sigma = noiselevel(V,'scans')

    assert abs(sigma - sigma_ref) < 1e-2
#============================================================


def test_complex():
#============================================================
    "Check estimation of noiselevel using a complex signal"

    t = np.linspace(0,3,200)
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[4, 0.4])
    lam = 0.25
    B = bg_exp(t,1.5)

    np.random.seed(1)
    noise = whitegaussnoise(t,0.03)
    np.random.seed(2)
    noisec = 1j*whitegaussnoise(t,0.03)
    V = dipolarkernel(t,r,mod=lam,bg=B)@P
    Vco = V*np.exp(-1j*np.pi/5)
    Vco = Vco + noise + noisec
    truelevel = np.std(noise)
    approxlevel = noiselevel(Vco,'complex')

    assert abs(truelevel - approxlevel) < 1e-2
#============================================================