
import numpy as np
from deerlab import dipolarkernel, whitegaussnoise, fitparamodel, bootan
from deerlab.dd_models import dd_gauss


def test_basics():
# ======================================================================
    "Check the basic functionality of the bootstrapping"

    t = np.linspace(0,5,200)
    r = np.linspace(2,6,300)
    P = dd_gauss(r,[4, 0.8])
    K = dipolarkernel(t,r)
    Vexp = K@P + whitegaussnoise(t,0.01,seed=1)

    par0 = [3, 0.5]
    Vmodel = lambda par: K@dd_gauss(r,par)
    fit = fitparamodel(Vexp,Vmodel,par0)
    Vfit = Vmodel(fit.param)


    def bootfcn(V):
        fit = fitparamodel(V,Vmodel,par0)
        return fit.param

    paruq = bootan(bootfcn,Vexp,Vfit,10)

    assert all(abs(paruq.mean - fit.param) < 1.5e-2)
# ======================================================================


def test_resampling():
# ======================================================================
    "Check that both bootstrap resampling method yield similar results"

    t = np.linspace(0,5,200)
    r = np.linspace(2,6,300)
    P = dd_gauss(r,[4, 0.8])
    K = dipolarkernel(t,r)
    Vexp = K@P + whitegaussnoise(t,0.01,seed=1)

    par0 = [3, 0.5]
    Vmodel = lambda par: K@dd_gauss(r,par)
    fit = fitparamodel(Vexp,Vmodel,par0)
    Vfit = Vmodel(fit.param)


    def bootfcn(V):
        fit = fitparamodel(V,Vmodel,par0)
        return fit.param

    paruq1 = bootan(bootfcn,Vexp,Vfit,5,resampling='residual')
    paruq2 = bootan(bootfcn,Vexp,Vfit,5,resampling='gaussian')


    assert all(abs(paruq1.mean - paruq2.mean) < 1.6e-2)
# ======================================================================


def test_multiple_ouputs():
# ======================================================================
    "Check that both bootstrap handles the correct number outputs"

    t = np.linspace(0,5,200)
    r = np.linspace(2,6,300)
    P = dd_gauss(r,[4, 0.8])
    K = dipolarkernel(t,r)
    Vexp = K@P + whitegaussnoise(t,0.01,seed=1)

    par0 = [3, 0.5]
    Vmodel = lambda par: K@dd_gauss(r,par)
    fit = fitparamodel(Vexp,Vmodel,par0)
    Vfit = Vmodel(fit.param)


    def bootfcn(V):
        fit = fitparamodel(V,Vmodel,par0)
        Pfit = dd_gauss(r,fit.param)
        return fit.param, Pfit

    paruq1 = bootan(bootfcn,Vexp,Vfit,4)


    assert len(paruq1)==2
# ======================================================================


def test_multiple_datasets():
# ======================================================================
    "Check bootstrapping when using multiple input datasets"

    t1 = np.linspace(0,5,200)
    t2 = np.linspace(-0.5,3,300)
    r = np.linspace(2,6,300)
    P = dd_gauss(r,[4, 0.8])
    K1 = dipolarkernel(t1,r)
    K2 = dipolarkernel(t2,r)

    Vexp1 = K1@P + whitegaussnoise(t1,0.01,seed=1)
    Vexp2 = K2@P + whitegaussnoise(t2,0.02,seed=2)

    def Vmodel(par):
        V1 = K1@dd_gauss(r,par)
        V2 = K2@dd_gauss(r,par)
        return [V1,V2]

    par0 = [3, 0.5]
    fit = fitparamodel([Vexp1,Vexp2],Vmodel,par0)
    Vfit1,Vfit2 = Vmodel(fit.param)

    def bootfcn(V):
        fit = fitparamodel(V,Vmodel,par0)
        return fit.param

    paruq = bootan(bootfcn,[Vexp1,Vexp2],[Vfit1,Vfit2],5)

    assert all(abs(paruq.mean - fit.param) < 1.5e-2)
# ======================================================================

def test_parallelization():
# ======================================================================
    "Check that bootan can run with multiple cores"

    t = np.linspace(0,5,200)
    r = np.linspace(2,6,300)
    P = dd_gauss(r,[4, 0.8])
    K = dipolarkernel(t,r)
    Vexp = K@P + whitegaussnoise(t,0.01)

    par0 = [3, 0.5]
    Vmodel = lambda par: K@dd_gauss(r,par)
    fit = fitparamodel(Vexp,Vmodel,par0)
    Vfit = Vmodel(fit.param)

    def bootfcn(V):
        fit = fitparamodel(V,Vmodel,par0)
        return fit.param

    paruq = bootan(bootfcn,Vexp,Vfit,10,cores=-1)

    assert all(abs(paruq.mean - fit.param) < 1.5e-2)
# ======================================================================

