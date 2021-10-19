from deerlab.dipolarkernel import dipolarkernel
from deerlab.whitegaussnoise import whitegaussnoise
import numpy as np
import matplotlib.pyplot as plt
from deerlab.model import Model,fit
from deerlab.dipolarmodel import dipolarmodel
from deerlab.utils import ovl
from deerlab import dd_gauss,dd_gauss2,bg_hom3d,bg_exp
import deerlab as dl 

# ======================================================================
def test_type(): 
    "Check that the function returns a valid model type"

    model = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1)

    assert isinstance(model,Model)
# ======================================================================

# ======================================================================
def test_Nparam_1(): 
    "Check that the function has the correct number of nonlinear parameters"

    model = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1)

    assert model.Nnonlin==5 and model.Nlin==1 and model.Nparam==6
# ======================================================================

# ======================================================================
def test_Nparam_2(): 
    "Check that the function has the correct number of nonlinear parameters"

    model = dipolarmodel(t,r,dd_gauss2,bg_hom3d,npathways=2)

    assert model.Nnonlin==9 and model.Nlin==2 and model.Nparam==11
# ======================================================================

# ======================================================================
def test_preservation(): 
    "Check that the inputs models are not modified by the output model"

    model = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=2)
    model.mean.par0 = 15

    assert model.mean.par0==15 and dd_gauss.mean.par0!=15
# ======================================================================

# ======================================================================
def test_names(): 
    "Check that the model has correct parameter names"

    model = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1)
    parameters = ['mean','width','conc','mod','reftime','scale']
    
    for param in parameters:
        assert hasattr(model,param)
# ======================================================================


t = np.linspace(-0.5,5,100)
r = np.linspace(2,5,50)
Bfcn = lambda t,lam: bg_hom3d(t,50,lam)
Bfcn_pheno = lambda t,_: bg_exp(t,0.1)
Pr = dd_gauss(r,3,0.2)
V1path = 1e5*dipolarkernel(t,r,mod=0.3,bg=Bfcn)@Pr
V1path_noB = 1e5*dipolarkernel(t,r,mod=0.3)@Pr
V1path_phenoB = 1e5*dipolarkernel(t,r,mod=0.3,bg=Bfcn_pheno)@Pr
V2path = 1e5*dipolarkernel(t,r,pathways=[[0.6],[0.3,0],[0.1,2]],bg=Bfcn)@Pr
V3path = 1e5*dipolarkernel(t,r,pathways=[[0.5],[0.3,0],[0.1,2],[0.1,5]],bg=Bfcn)@Pr


# ======================================================================
def test_call_positional(): 
    "Check that the model called via positional arguments responds correctly"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1)
    
    Vsim = Vmodel(0.3,0.0,50,3,0.2,1e5)

    assert np.allclose(Vsim,V1path)
# ======================================================================

# ======================================================================
def test_call_keywords(): 
    "Check that the model called via keyword arguments responds correctly"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1)
    
    Vsim = Vmodel(mod=0.3,reftime=0.0,conc=50,mean=3,width=0.2,scale=1e5)

    assert np.allclose(Vsim,V1path)
# ======================================================================

# ======================================================================
def test_phenomenological_Bmodel(): 
    "Check model generation of a dipolar signal with a phenomelogical background"

    Vmodel = dipolarmodel(t,r,dd_gauss,Bmodel=bg_exp,npathways=1)
    
    Vsim = Vmodel(mod=0.3,reftime=0.0,mean=3,width=0.2,decay=0.1,scale=1e5)

    assert np.allclose(Vsim,V1path_phenoB)
# ======================================================================

# ======================================================================
def test_no_Bmodel(): 
    "Check model generation of a dipolar signal without background"

    Vmodel = dipolarmodel(t,r,dd_gauss,None,npathways=1)
    
    Vsim = Vmodel(mod=0.3,reftime=0.0,mean=3,width=0.2,scale=1e5)

    assert np.allclose(Vsim,V1path_noB)
# ======================================================================

# ======================================================================
def test_model_1pathways(): 
    "Check that the model with one dipolar pathway is correct"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1)
    
    Vsim = Vmodel(mod=0.3,reftime=0.0,conc=50,mean=3,width=0.2,scale=1e5)

    assert np.allclose(Vsim,V1path)
# ======================================================================

# ======================================================================
def test_model_2pathways(): 
    "Check that the model with two dipolar pathways is correct"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=2)
    
    Vsim = Vmodel(lam1=0.3,reftime1=0.0,lam2=0.1,
                    reftime2=2,conc=50,mean=3,width=0.2,scale=1e5)

    assert np.allclose(Vsim,V2path)
# ======================================================================

# ======================================================================
def test_model_3pathways(): 
    "Check that the model with three dipolar pathways is correct"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=3)
    
    Vsim = Vmodel(lam1=0.3,reftime1=0.0,
                lam2=0.1,reftime2=2, lam3=0.1, reftime3=5,
                conc=50,mean=3,width=0.2,scale=1e5)

    assert np.allclose(Vsim,V3path)
# ======================================================================


# ======================================================================
def test_fit_1pathways(): 
    "Check that the model can be correctly fitted with one dipolar pathway"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1)
    
    result = fit(Vmodel,V1path,nonlin_tol=1e-3)

    assert np.allclose(result.model,V1path)
# ======================================================================

# ======================================================================
def test_fit_2pathways(): 
    "Check that the model can be correctly fitted with two dipolar pathways"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=2)
    Vmodel.reftime1.freeze(0)
    Vmodel.reftime2.freeze(2)
    
    result = fit(Vmodel,V2path,nonlin_tol=1e-3)

    assert np.allclose(result.model,V2path)
# ======================================================================


# ======================================================================
def test_fit_3pathways(): 
    "Check that the model can be correctly fitted with three dipolar pathways"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=3)
    Vmodel.reftime1.freeze(0)
    Vmodel.reftime2.freeze(2)
    Vmodel.reftime3.freeze(5)
    
    result = fit(Vmodel,V3path,nonlin_tol=1e-3)

    assert np.allclose(result.model,V3path)
# ======================================================================


V1harm = 1e5*dipolarkernel(t,r,pathways=[[0.7],[0.3,0,1]],bg=Bfcn)@Pr
V2harm = 1e5*dipolarkernel(t,r,pathways=[[0.6],[0.3,0,1],[0.1,2,2]],bg=Bfcn)@Pr
V3harm = 1e5*dipolarkernel(t,r,pathways=[[0.5],[0.3,0,1],[0.1,2,2],[0.1,5,3]],bg=Bfcn)@Pr


# ======================================================================
def test_model_1harmonics(): 
    "Check that the model with one harmonic is correct"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1,harmonics=1)
    
    Vsim = Vmodel(mod=0.3,reftime=0.0,conc=50,mean=3,width=0.2,scale=1e5)

    assert np.allclose(Vsim,V1harm)
# ======================================================================

# ======================================================================
def test_model_2harmonics(): 
    "Check that the model with two different harmonics is correct"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=2,harmonics=[1,2])
    
    Vsim = Vmodel(lam1=0.3,reftime1=0.0,lam2=0.1,
                    reftime2=2,conc=50,mean=3,width=0.2,scale=1e5)

    assert np.allclose(Vsim,V2harm)
# ======================================================================

# ======================================================================
def test_model_3harmonics(): 
    "Check that the model with three different harmonics is correct"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=3,harmonics=[1,2,3])
    
    Vsim = Vmodel(lam1=0.3,reftime1=0.0,
                lam2=0.1,reftime2=2, lam3=0.1, reftime3=5,
                conc=50,mean=3,width=0.2,scale=1e5)

    assert np.allclose(Vsim,V3harm)
# ======================================================================

# ======================================================================
def test_call_Pnonparametric(): 
    "Check that the model with one dipolar pathway is correct"

    Vmodel = dipolarmodel(t,r,Bmodel=bg_hom3d,npathways=1)
    Vsim = Vmodel(mod=0.3,reftime=0.0,conc=50,P=1e5*dd_gauss(r,3,0.2))

    assert np.allclose(Vsim,V1path)
# ======================================================================

# ======================================================================
def test_fit_Pnonparametric(): 
    "Check that the model with one dipolar pathway is correct"

    Vmodel = dipolarmodel(t,r,Bmodel=bg_hom3d,npathways=1)
    
    result = fit(Vmodel,V1path,nonlin_tol=1e-3)

    assert np.allclose(result.model,V1path,atol=1e-2) and ovl(result.P/1e5,Pr)>0.975
# ======================================================================

# ======================================================================
def test_compactness_penalty_Pnonparametric(): 
    "Check the fitting with a nonparametric distribution and the compactness penalty"

    Vmodel = dipolarmodel(t,r,compactness=True)
    Vmodel.compactness.weight.freeze(0.05)

    result = fit(Vmodel,V1path+whitegaussnoise(t,0.01,seed=1),nonlin_tol=1e-3)

    assert ovl(result.P/1e5,Pr)>0.975
# ======================================================================

# ======================================================================
def test_compactness_penalty_Pparametric(): 
    "Check the fitting with a parametric distribution and the compactness penalty"

    Vmodel = dipolarmodel(t,r,Pmodel=dd_gauss,compactness=True)
    Vmodel.compactness.weight.freeze(0.05)
    
    result = fit(Vmodel,V1path+whitegaussnoise(t,0.01,seed=1),nonlin_tol=1e-3)

    assert ovl(result.evaluate(dd_gauss,r),Pr)>0.975
# ======================================================================

# ======================================================================
def test_smoothness_penalty_Pnonparametric(): 
    "Check the fitting with a nonparametric distribution and the smoothness penalty"

    Vmodel = dipolarmodel(t,r,smoothness=True)
    Vmodel.smoothness.weight.freeze(0.00005)

    result = fit(Vmodel,V1path+whitegaussnoise(t,0.01,seed=1),nonlin_tol=1e-3)

    assert ovl(result.P/1e5,Pr)>0.975
# ======================================================================

# ======================================================================
def test_compactness_penalty_Pparametric(): 
    "Check the fitting with a parametric distribution and the smoothness penalty"

    Vmodel = dipolarmodel(t,r,Pmodel=dd_gauss,smoothness=True)
    Vmodel.smoothness.weight.freeze(0.005)
    
    result = fit(Vmodel,V1path+whitegaussnoise(t,0.01,seed=1),nonlin_tol=1e-3)

    assert ovl(result.evaluate(dd_gauss,r),Pr)>0.975
# ======================================================================
