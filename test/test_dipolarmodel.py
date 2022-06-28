from deerlab.dipolarkernel import dipolarkernel
from deerlab.utils.utils import ovl
from deerlab.whitegaussnoise import whitegaussnoise
import numpy as np
import matplotlib.pyplot as plt
from deerlab.model import Model,fit
from deerlab.dipolarmodel import ExperimentInfo,dipolarpenalty, dipolarmodel, ex_4pdeer, ex_3pdeer,ex_fwd5pdeer, ex_rev5pdeer, ex_sifter, ex_ridme
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
    parameters = ['mean','std','conc','mod','reftime','scale']
    
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
    
    Vsim = Vmodel(mod=0.3,reftime=0.0,conc=50,mean=3,std=0.2,scale=1e5)

    assert np.allclose(Vsim,V1path)
# ======================================================================

# ======================================================================
def test_phenomenological_Bmodel(): 
    "Check model generation of a dipolar signal with a phenomelogical background"

    Vmodel = dipolarmodel(t,r,dd_gauss,Bmodel=bg_exp,npathways=1)
    
    Vsim = Vmodel(mod=0.3,reftime=0.0,mean=3,std=0.2,decay=0.1,scale=1e5)

    assert np.allclose(Vsim,V1path_phenoB)
# ======================================================================

# ======================================================================
def test_no_Bmodel(): 
    "Check model generation of a dipolar signal without background"

    Vmodel = dipolarmodel(t,r,dd_gauss,None,npathways=1)
    
    Vsim = Vmodel(mod=0.3,reftime=0.0,mean=3,std=0.2,scale=1e5)

    assert np.allclose(Vsim,V1path_noB)
# ======================================================================

# ======================================================================
def test_model_1pathways(): 
    "Check that the model with one dipolar pathway is correct"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1)
    
    Vsim = Vmodel(mod=0.3,reftime=0.0,conc=50,mean=3,std=0.2,scale=1e5)

    assert np.allclose(Vsim,V1path)
# ======================================================================

# ======================================================================
def test_model_2pathways(): 
    "Check that the model with two dipolar pathways is correct"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=2)
    
    Vsim = Vmodel(lam1=0.3,reftime1=0.0,lam2=0.1,
                    reftime2=2,conc=50,mean=3,std=0.2,scale=1e5)

    assert np.allclose(Vsim,V2path)
# ======================================================================

# ======================================================================
def test_model_3pathways(): 
    "Check that the model with three dipolar pathways is correct"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=3)
    
    Vsim = Vmodel(lam1=0.3,reftime1=0.0,
                lam2=0.1,reftime2=2, lam3=0.1, reftime3=5,
                conc=50,mean=3,std=0.2,scale=1e5)

    assert np.allclose(Vsim,V3path)
# ======================================================================


# ======================================================================
def test_fit_1pathways(): 
    "Check that the model can be correctly fitted with one dipolar pathway"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1)
    
    result = fit(Vmodel,V1path,ftol=1e-4)

    assert np.allclose(result.model,V1path)
# ======================================================================

# ======================================================================
def test_fit_2pathways(): 
    "Check that the model can be correctly fitted with two dipolar pathways"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=2)
    Vmodel.reftime1.freeze(0)
    Vmodel.reftime2.freeze(2)
    
    result = fit(Vmodel,V2path,ftol=1e-4)

    assert np.allclose(result.model,V2path)
# ======================================================================

# ======================================================================
def test_fit_3pathways(): 
    "Check that the model can be correctly fitted with three dipolar pathways"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=3)
    Vmodel.reftime1.freeze(0)
    Vmodel.reftime2.freeze(2)
    Vmodel.reftime3.freeze(5)
    
    result = fit(Vmodel,V3path,ftol=1e-4)

    assert np.allclose(result.model,V3path)
# ======================================================================


V1harm = 1e5*dipolarkernel(t,r,pathways=[[0.7],[0.3,0,1]],bg=Bfcn)@Pr
V2harm = 1e5*dipolarkernel(t,r,pathways=[[0.6],[0.3,0,1],[0.1,2,2]],bg=Bfcn)@Pr
V3harm = 1e5*dipolarkernel(t,r,pathways=[[0.5],[0.3,0,1],[0.1,2,2],[0.1,5,3]],bg=Bfcn)@Pr


# ======================================================================
def test_model_1harmonics(): 
    "Check that the model with one harmonic is correct"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1,harmonics=1)
    
    Vsim = Vmodel(mod=0.3,reftime=0.0,conc=50,mean=3,std=0.2,scale=1e5)

    assert np.allclose(Vsim,V1harm)
# ======================================================================

# ======================================================================
def test_model_2harmonics(): 
    "Check that the model with two different harmonics is correct"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=2,harmonics=[1,2])
    
    Vsim = Vmodel(lam1=0.3,reftime1=0.0,lam2=0.1,
                    reftime2=2,conc=50,mean=3,std=0.2,scale=1e5)

    assert np.allclose(Vsim,V2harm)
# ======================================================================

# ======================================================================
def test_model_3harmonics(): 
    "Check that the model with three different harmonics is correct"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=3,harmonics=[1,2,3])
    
    Vsim = Vmodel(lam1=0.3,reftime1=0.0,
                lam2=0.1,reftime2=2, lam3=0.1, reftime3=5,
                conc=50,mean=3,std=0.2,scale=1e5)

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
    
    result = fit(Vmodel,V1path,ftol=1e-5)

    assert ovl(result.P,Pr)>0.975
# ======================================================================

# ======================================================================
def test_fit_Pnonparametric_normalization(): 
    "Check that the model with one dipolar pathway is correct"

    Vmodel = dipolarmodel(t,r,Bmodel=bg_hom3d,npathways=1)
    
    result = fit(Vmodel,V1path,ftol=1e-5)

    assert np.isclose(np.trapz(result.P,r),1)
# ======================================================================

tdeer = np.linspace(-0.5,5,300)
tsifter = np.linspace(-2,4,300)
tridme = np.linspace(-1,3,300)
tau1,tau2,tau3 = 1,2,3
V3pdeer = 1e5*dipolarkernel(tdeer,r,pathways=[[0.6],[0.3,0],[0.1,tau1]],bg=Bfcn)@Pr
V4pdeer = 1e5*dipolarkernel(tdeer,r,pathways=[[0.6],[0.3,tau1],[0.1,tau1+tau2]],bg=Bfcn)@Pr
Vfwd5pdeer = 1e5*dipolarkernel(tdeer,r,pathways=[[0.6],[0.3,tau3],[0.1,tau1]],bg=Bfcn)@Pr
Vrev5pdeer = 1e5*dipolarkernel(tdeer,r,pathways=[[0.6],[0.3,tau3],[0.1,tau2]],bg=Bfcn)@Pr
Vsifter = 1e5*dipolarkernel(tsifter,r,pathways=[[0.3],[0.5,tau2-tau1,1],[0.1,2*tau2,1/2],[0.1,-2*tau1,1/2]],bg=Bfcn)@Pr
Vridme  = 1e5*dipolarkernel(tridme,r,pathways=[[0.3],[0.5,0],[0.1,tau2],[0.1,-tau1]],bg=Bfcn)@Pr

# ======================================================================
def test_ex_3pdeer_type(): 
    "Check the 3-pulse DEER experimental model."

    experiment = ex_3pdeer(tau1)

    assert isinstance(experiment,ExperimentInfo) 
# ======================================================================

# ======================================================================
def test_ex_3pdeer_fit(): 
    "Check the 3-pulse DEER experimental model."

    experiment = ex_3pdeer(tau1, pathways=[1,2])
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment)
    result = fit(Vmodel,V3pdeer,ftol=1e-5)

    assert np.allclose(V3pdeer,result.model,atol=1e-1) and ovl(result.P/1e5,Pr)>0.975
# ======================================================================

# ======================================================================
def test_ex_4pdeer_type(): 
    "Check the 4-pulse DEER experimental model."

    experiment = ex_4pdeer(tau1,tau2)

    assert isinstance(experiment,ExperimentInfo) 
# ======================================================================

# ======================================================================
def test_ex_4pdeer_fit(): 
    "Check the 4-pulse DEER experimental model."

    experiment = ex_4pdeer(tau1,tau2, pathways=[1,2])
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment)
    result = fit(Vmodel,V4pdeer,ftol=1e-4)

    assert np.allclose(V4pdeer,result.model,atol=1e-2) and ovl(result.P/1e5,Pr)>0.975
# ======================================================================

# ======================================================================
def test_ex_rev5pdeer_type(): 
    "Check the reverse 5-pulse DEER experimental model."

    experiment = ex_rev5pdeer(tau1,tau2,tau3)

    assert isinstance(experiment,ExperimentInfo) 
# ======================================================================

# ======================================================================
def test_ex_rev5pdeer_fit(): 
    "Check the reverse 5-pulse DEER experimental model in fitting."

    experiment = ex_rev5pdeer(tau1,tau2,tau3, pathways=[1,2,3])
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment)
    result = fit(Vmodel,Vrev5pdeer,ftol=1e-4)

    assert np.allclose(Vrev5pdeer,result.model,atol=1e-2) and ovl(result.P/1e5,Pr)>0.975
# ======================================================================


# ======================================================================
def test_ex_fwd5pdeer_type(): 
    "Check the forward 5-pulse DEER experimental model."

    experiment = ex_fwd5pdeer(tau1,tau2,tau3)

    assert isinstance(experiment,ExperimentInfo) 
# ======================================================================

# ======================================================================
def test_ex_fwd5pdeer_fit(): 
    "Check the forward 5-pulse DEER experimental model in fitting."

    experiment = ex_fwd5pdeer(tau1,tau2,tau3, pathways=[1,2,3])
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment)
    result = fit(Vmodel,Vfwd5pdeer,ftol=1e-4)

    assert np.allclose(Vfwd5pdeer,result.model,atol=1e-2) and ovl(result.P/1e5,Pr)>0.975
# ======================================================================

# ======================================================================
def test_ex_sifter_type(): 
    "Check the SIFTER experimental model."

    experiment = ex_sifter(tau1,tau2)

    assert isinstance(experiment,ExperimentInfo) 
# ======================================================================

# ======================================================================
def test_ex_sifter_fit(): 
    "Check the SIFTER experimental model in fitting."

    experiment = ex_sifter(tau1,tau2)
    Vmodel = dipolarmodel(tsifter,r,Bmodel=bg_hom3d,experiment=experiment)
    result = fit(Vmodel,Vsifter,ftol=1e-4)

    assert np.allclose(Vsifter,result.model,atol=1e-2) and ovl(result.P/1e5,Pr)>0.975
# ======================================================================

# ======================================================================
def test_ex_ridme_type(): 
    "Check the RIDME experimental model."

    experiment = ex_ridme(tau1,tau2)

    assert isinstance(experiment,ExperimentInfo) 
# ======================================================================

# ======================================================================
def test_ex_ridme_fit(): 
    "Check the RIDME experimental model in fitting."

    experiment = ex_ridme(tau1,tau2, pathways=[1,2,3])
    Vmodel = dipolarmodel(tridme,r,Bmodel=bg_hom3d,experiment=experiment)
    result = fit(Vmodel,Vridme,ftol=1e-4)

    assert np.allclose(Vridme,result.model,atol=1e-2) and ovl(result.P/1e5,Pr)>0.975
# ======================================================================

# ======================================================================
def test_orisel(): 
    "Check that dipolar models with orientation selection work"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1)
    Vmodelorisel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1,orisel=lambda theta: np.ones_like(theta))

    Vref = Vmodel(mean=3,std=0.2,mod=0.3,reftime=0,conc=200,scale=1e2)
    Vorisel = Vmodelorisel(mean=3,std=0.2,mod=0.3,reftime=0,conc=200,scale=1e2)

    assert np.allclose(Vref,Vorisel,rtol=1e-4)
# ======================================================================


# ======================================================================
def test_excitationbandwidth(): 
    "Check that dipolar models with limited excitation bandwidth work"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1)
    Vmodelorisel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1,excbandwidth=1e8)

    Vref = Vmodel(mean=3,std=0.2,mod=0.3,reftime=0,conc=200,scale=1e2)
    Vorisel = Vmodelorisel(mean=3,std=0.2,mod=0.3,reftime=0,conc=200,scale=1e2)

    assert np.allclose(Vref,Vorisel,rtol=1e-4)
# ======================================================================

