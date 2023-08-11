from deerlab.dipolarkernel import dipolarkernel
from deerlab.utils import ovl
import numpy as np
from deerlab.model import Model
from deerlab.fit import fit
from deerlab.dipolarmodel import ExperimentInfo,dipolarpenalty, dipolarmodel, ex_4pdeer, ex_3pdeer,ex_fwd5pdeer, ex_rev5pdeer, ex_sifter, ex_ridme, ex_dqc
from deerlab import dd_gauss,dd_gauss2,bg_hom3d,bg_exp
from deerlab.utils import assert_docstring
from pytest import fixture 

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

# ======================================================================
def test_pathway_labelling(): 
    "Check that the model has correct parameter names"

    exp = ex_4pdeer(1,2,pathways=[1])
    model = dipolarmodel(t,r,dd_gauss,bg_hom3d,experiment=exp)
    parameters = ['mean','std','conc','mod','reftime','scale']
    
    for param in parameters:
        assert hasattr(model,'reftime') and not hasattr(model,'reftime1') and hasattr(model,'mod') and not hasattr(model,'lam1')
# ======================================================================

# ======================================================================
def test_pathway_labelling2(): 
    "Check that the model has correct parameter names"

    exp = ex_4pdeer(1,2,pathways=[1,2])
    model = dipolarmodel(t,r,dd_gauss,bg_hom3d,experiment=exp)
    parameters = ['mean','std','conc','mod','reftime','scale']
    
    for param in parameters:
        assert hasattr(model,'reftime1') and hasattr(model,'reftime2') and hasattr(model,'lam2') and hasattr(model,'lam2')
# ======================================================================

# ======================================================================
def test_pathway_labelling3(): 
    "Check that the model has correct parameter names"

    exp = ex_4pdeer(1,2,pathways=[2,4])
    model = dipolarmodel(t,r,dd_gauss,bg_hom3d,experiment=exp)
    parameters = ['mean','std','conc','mod','reftime','scale']
    
    for param in parameters:
        assert hasattr(model,'reftime2') and hasattr(model,'reftime4') and hasattr(model,'lam2') and hasattr(model,'lam4')
# ======================================================================

# Fixtures
# ----------------------------------------------------------------------
t = np.linspace(-0.5,5,100)
r = np.linspace(2,5,50)
@fixture(scope='module')
def Bfcn():
    return lambda t,lam: bg_hom3d(t,50,lam)
@fixture(scope='module')
def Bfcn_pheno():
    return lambda t,_: bg_exp(t,0.1)
@fixture(scope='module')
def Pr():
    return dd_gauss(r,3,0.2)
@fixture(scope='module')
def V1path(Pr,Bfcn):
    return 1e5*dipolarkernel(t,r,mod=0.3,bg=Bfcn)@Pr
@fixture(scope='module')
def V1path_noB(Pr): 
    return 1e5*dipolarkernel(t,r,mod=0.3)@Pr
@fixture(scope='module')
def V1path_phenoB(Pr,Bfcn_pheno): 
    return 1e5*dipolarkernel(t,r,mod=0.3,bg=Bfcn_pheno)@Pr 
@fixture(scope='module')
def V2path(Pr,Bfcn): 
    return 1e5*dipolarkernel(t,r,pathways=[{'amp':0.6},{'amp':0.3,'reftime':0},{'amp':0.1,'reftime':2}],bg=Bfcn)@Pr
@fixture(scope='module')
def V3path(Pr,Bfcn): 
    return 1e5*dipolarkernel(t,r,pathways=[{'amp':0.5},{'amp':0.3,'reftime':0},{'amp':0.1,'reftime':2},{'amp':0.1,'reftime':5}],bg=Bfcn)@Pr
# ----------------------------------------------------------------------

# ======================================================================
def test_call_positional(V1path): 
    "Check that the model called via positional arguments responds correctly"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1)
    
    Vsim = Vmodel(0.3,0.0,50,3,0.2,1e5)

    assert np.allclose(Vsim,V1path)
# ======================================================================

# ======================================================================
def test_call_keywords(V1path): 
    "Check that the model called via keyword arguments responds correctly"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1)
    
    Vsim = Vmodel(mod=0.3,reftime=0.0,conc=50,mean=3,std=0.2,scale=1e5)

    assert np.allclose(Vsim,V1path)
# ======================================================================

# ======================================================================
def test_phenomenological_Bmodel(V1path_phenoB): 
    "Check model generation of a dipolar signal with a phenomelogical background"

    Vmodel = dipolarmodel(t,r,dd_gauss,Bmodel=bg_exp,npathways=1)
    
    Vsim = Vmodel(mod=0.3,reftime=0.0,mean=3,std=0.2,decay=0.1,scale=1e5)

    assert np.allclose(Vsim,V1path_phenoB)
# ======================================================================

# ======================================================================
def test_no_Bmodel(V1path_noB): 
    "Check model generation of a dipolar signal without background"

    Vmodel = dipolarmodel(t,r,dd_gauss,None,npathways=1)
    
    Vsim = Vmodel(mod=0.3,reftime=0.0,mean=3,std=0.2,scale=1e5)

    assert np.allclose(Vsim,V1path_noB)
# ======================================================================

# ======================================================================
def test_model_1pathways(V1path): 
    "Check that the model with one dipolar pathway is correct"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1)
    
    Vsim = Vmodel(mod=0.3,reftime=0.0,conc=50,mean=3,std=0.2,scale=1e5)

    assert np.allclose(Vsim,V1path)
# ======================================================================

# ======================================================================
def test_model_2pathways(V2path): 
    "Check that the model with two dipolar pathways is correct"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=2)
    
    Vsim = Vmodel(lam1=0.3,reftime1=0.0,lam2=0.1,
                    reftime2=2,conc=50,mean=3,std=0.2,scale=1e5)

    assert np.allclose(Vsim,V2path)
# ======================================================================

# ======================================================================
def test_model_3pathways(V3path): 
    "Check that the model with three dipolar pathways is correct"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=3)
    
    Vsim = Vmodel(lam1=0.3,reftime1=0.0,
                lam2=0.1,reftime2=2, lam3=0.1, reftime3=5,
                conc=50,mean=3,std=0.2,scale=1e5)

    assert np.allclose(Vsim,V3path)
# ======================================================================


# ======================================================================
def test_fit_1pathways(V1path): 
    "Check that the model can be correctly fitted with one dipolar pathway"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1)
    
    result = fit(Vmodel,V1path,ftol=1e-4)

    assert np.allclose(result.model,V1path)
# ======================================================================

# ======================================================================
def test_fit_2pathways(V2path): 
    "Check that the model can be correctly fitted with two dipolar pathways"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=2)
    Vmodel.reftime1.freeze(0)
    Vmodel.reftime2.freeze(2)
    
    result = fit(Vmodel,V2path,ftol=1e-4)

    assert np.allclose(result.model,V2path)
# ======================================================================

# ======================================================================
def test_fit_3pathways(V3path): 
    "Check that the model can be correctly fitted with three dipolar pathways"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=3)
    Vmodel.reftime1.freeze(0)
    Vmodel.reftime2.freeze(2)
    Vmodel.reftime3.freeze(5)
    
    result = fit(Vmodel,V3path,ftol=1e-4)

    assert np.allclose(result.model,V3path)
# ======================================================================

# Fixtures
# ----------------------------------------------------------------------
@fixture(scope='module')
def V1harm(Pr,Bfcn):
    return 1e5*dipolarkernel(t,r,pathways=[{'amp':0.7},{'amp':0.3,'reftime':0,'harmonic':1}],bg=Bfcn)@Pr
@fixture(scope='module')
def V2harm(Pr,Bfcn): 
    return 1e5*dipolarkernel(t,r,pathways=[{'amp':0.6},{'amp':0.3,'reftime':0,'harmonic':1},{'amp':0.1,'reftime':2,'harmonic':2}],bg=Bfcn)@Pr
@fixture(scope='module')
def V3harm(Pr,Bfcn): 
    return 1e5*dipolarkernel(t,r,pathways=[{'amp':0.5},{'amp':0.3,'reftime':0,'harmonic':1},{'amp':0.1,'reftime':2,'harmonic':2},{'amp':0.1,'reftime':5,'harmonic':3}],bg=Bfcn)@Pr
# ----------------------------------------------------------------------


# ======================================================================
def test_model_1harmonics(V1harm): 
    "Check that the model with one harmonic is correct"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1,harmonics=1)
    
    Vsim = Vmodel(mod=0.3,reftime=0.0,conc=50,mean=3,std=0.2,scale=1e5)

    assert np.allclose(Vsim,V1harm)
# ======================================================================

# ======================================================================
def test_model_2harmonics(V2harm): 
    "Check that the model with two different harmonics is correct"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=2,harmonics=[1,2])
    
    Vsim = Vmodel(lam1=0.3,reftime1=0.0,lam2=0.1,
                    reftime2=2,conc=50,mean=3,std=0.2,scale=1e5)

    assert np.allclose(Vsim,V2harm)
# ======================================================================

# ======================================================================
def test_model_3harmonics(V3harm): 
    "Check that the model with three different harmonics is correct"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=3,harmonics=[1,2,3])
    
    Vsim = Vmodel(lam1=0.3,reftime1=0.0,
                lam2=0.1,reftime2=2, lam3=0.1, reftime3=5,
                conc=50,mean=3,std=0.2,scale=1e5)

    assert np.allclose(Vsim,V3harm)
# ======================================================================

# ======================================================================
def test_call_Pnonparametric(V1path): 
    "Check that the model with one dipolar pathway is correct"

    Vmodel = dipolarmodel(t,r,Bmodel=bg_hom3d,npathways=1)
    Vsim = Vmodel(mod=0.3,reftime=0.0,conc=50,P=1e5*dd_gauss(r,3,0.2))

    assert np.allclose(Vsim,V1path)
# ======================================================================

# ======================================================================
def test_fit_Pnonparametric(V1path,Pr): 
    "Check that the model with one dipolar pathway is correct"

    Vmodel = dipolarmodel(t,r,Bmodel=bg_hom3d,npathways=1)
    
    result = fit(Vmodel,V1path,ftol=1e-5)

    assert ovl(result.P,Pr)>0.975
# ======================================================================

# ======================================================================
def test_fit_Pnonparametric_normalization(V1path): 
    "Check that the model with one dipolar pathway is correct"

    Vmodel = dipolarmodel(t,r,Bmodel=bg_hom3d,npathways=1)
    
    result = fit(Vmodel,V1path,ftol=1e-5)

    assert np.isclose(np.trapz(result.P,r),1)
# ======================================================================

# Fixtures
# ----------------------------------------------------------------------
tdeer = np.linspace(-0.5,5,300)
tsifter = np.linspace(-2,4,300)
tridme = np.linspace(0,3,300)
tau1,tau2,tau3 = 1,2,3
@fixture(scope='module')
def V3pdeer(Pr,Bfcn): 
    return 1e5*dipolarkernel(tdeer,r,pathways=[{'amp':0.6},{'amp':0.3,'reftime':0},{'amp':0.1,'reftime':tau1}],bg=Bfcn)@Pr
@fixture(scope='module')
def V4pdeer(Pr,Bfcn): 
    return 1e5*dipolarkernel(tdeer,r,pathways=[{'amp':0.6},{'amp':0.3,'reftime':tau1},{'amp':0.1,'reftime':tau1+tau2}],bg=Bfcn)@Pr
@fixture(scope='module')
def Vfwd5pdeer(Pr,Bfcn): 
    return 1e5*dipolarkernel(tdeer,r,pathways=[{'amp':0.6},{'amp':0.3,'reftime':tau3},{'amp':0.1,'reftime':tau1}],bg=Bfcn)@Pr
@fixture(scope='module')
def Vrev5pdeer(Pr,Bfcn): 
    return 1e5*dipolarkernel(tdeer,r,pathways=[{'amp':0.6},{'amp':0.3,'reftime':tau3},{'amp':0.1,'reftime':tau2}],bg=Bfcn)@Pr
@fixture(scope='module')
def Vsifter(Pr,Bfcn): 
    return 1e5*dipolarkernel(tsifter,r,pathways=[{'amp':0.3},{'amp':0.5,'reftime':tau2-tau1,'harmonic':1},{'amp':0.1,'reftime':2*tau2,'harmonic':1/2},{'amp':0.1,'reftime':-2*tau1,'harmonic':1/2}],bg=Bfcn)@Pr
@fixture(scope='module')
def Vridme(Pr,Bfcn): 
    return 1e5*dipolarkernel(tridme,r,pathways=[{'amp':0.3},{'amp':0.5,'reftime':tau1},{'amp':0.1,'reftime':tau1+tau2},{'amp':0.1,'reftime':0}],bg=Bfcn)@Pr
# ----------------------------------------------------------------------

# ======================================================================
def test_ex_3pdeer_type(): 
    "Check the 3-pulse DEER experimental model."

    experiment = ex_3pdeer(tau1)

    assert isinstance(experiment,ExperimentInfo) 
# ======================================================================

# ======================================================================
def test_ex_3pdeer_fit(V3pdeer,Pr): 
    "Check the 3-pulse DEER experimental model."

    experiment = ex_3pdeer(tau1, pathways=[1,2])
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment)
    result = fit(Vmodel,V3pdeer,ftol=1e-5)

    assert np.allclose(V3pdeer,result.model,rtol=1e-3) and ovl(result.P/1e5,Pr)>0.975
# ======================================================================

# ======================================================================
def test_ex_4pdeer_type(): 
    "Check the 4-pulse DEER experimental model."

    experiment = ex_4pdeer(tau1,tau2)

    assert isinstance(experiment,ExperimentInfo) 
# ======================================================================

# ======================================================================
def test_ex_4pdeer_fit(V4pdeer,Pr): 
    "Check the 4-pulse DEER experimental model."

    experiment = ex_4pdeer(tau1,tau2, pathways=[1,2])
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment)
    result = fit(Vmodel,V4pdeer,ftol=1e-4)

    assert np.allclose(V4pdeer,result.model,rtol=1e-3) and ovl(result.P/1e5,Pr)>0.975
# ======================================================================

# ======================================================================
def test_ex_4pdeer_pulselength(): 
    "Check the 4-pulse DEER experimental model."

    pulselength = 0.1
    experiment = ex_4pdeer(tau1,tau2, pathways=[1], pulselength=pulselength)
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment)

    assert np.isclose(Vmodel.reftime.ub - Vmodel.reftime.lb,6*pulselength)
# ======================================================================

# ======================================================================
def test_ex_rev5pdeer_type(): 
    "Check the reverse 5-pulse DEER experimental model."

    experiment = ex_rev5pdeer(tau1,tau2,tau3)

    assert isinstance(experiment,ExperimentInfo) 
# ======================================================================

# ======================================================================
def test_ex_rev5pdeer_fit(Vrev5pdeer,Pr): 
    "Check the reverse 5-pulse DEER experimental model in fitting."

    experiment = ex_rev5pdeer(tau1,tau2,tau3, pathways=[1,5,3])
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment)
    result = fit(Vmodel,Vrev5pdeer,ftol=1e-4)

    assert np.allclose(Vrev5pdeer,result.model,rtol=1e-3) and ovl(result.P/1e5,Pr)>0.975
# ======================================================================

# ======================================================================
def test_ex_rev5pdeer_pulselength(): 
    "Check the 5-pulse DEER experimental model."

    pulselength = 0.1
    experiment = ex_rev5pdeer(tau1,tau2,tau3, pathways=[1], pulselength=pulselength)
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment)

    assert np.isclose(Vmodel.reftime.ub - Vmodel.reftime.lb,6*pulselength)
# ======================================================================

# ======================================================================
def test_ex_fwd5pdeer_type(): 
    "Check the forward 5-pulse DEER experimental model."

    experiment = ex_fwd5pdeer(tau1,tau2,tau3)

    assert isinstance(experiment,ExperimentInfo) 
# ======================================================================

# ======================================================================
def test_ex_fwd5pdeer_fit(Vfwd5pdeer,Pr): 
    "Check the forward 5-pulse DEER experimental model in fitting."

    experiment = ex_fwd5pdeer(tau1,tau2,tau3, pathways=[1,5,2])
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment)
    result = fit(Vmodel,Vfwd5pdeer,ftol=1e-4)

    assert np.allclose(Vfwd5pdeer,result.model,rtol=1e-3) and ovl(result.P/1e5,Pr)>0.975
# ======================================================================

# ======================================================================
def test_ex_fwd5pdeer_pulselength(): 
    "Check the forward 5-pulse DEER experimental model."

    pulselength = 0.1
    experiment = ex_fwd5pdeer(tau1,tau2,tau3, pathways=[1], pulselength=pulselength)
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment)

    assert np.isclose(Vmodel.reftime.ub - Vmodel.reftime.lb,6*pulselength)
# ======================================================================

# ======================================================================
def test_ex_sifter_type(): 
    "Check the SIFTER experimental model."

    experiment = ex_sifter(tau1,tau2)

    assert isinstance(experiment,ExperimentInfo) 
# ======================================================================

# ======================================================================
def test_ex_sifter_fit(Vsifter,Pr): 
    "Check the SIFTER experimental model in fitting."

    experiment = ex_sifter(tau1,tau2)
    Vmodel = dipolarmodel(tsifter,r,Bmodel=bg_hom3d,experiment=experiment)
    result = fit(Vmodel,Vsifter,ftol=1e-4)

    assert np.allclose(Vsifter,result.model,rtol=1e-3) and ovl(result.P/1e5,Pr)>0.975
# ======================================================================

# ======================================================================
def test_ex_sifter_pulselength(): 
    "Check the SIFTER experimental model."

    pulselength = 0.1
    experiment = ex_sifter(tau1,tau2,pathways=[1], pulselength=pulselength)
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment)

    assert np.isclose(Vmodel.reftime.ub - Vmodel.reftime.lb,6*pulselength)
# ======================================================================

# ======================================================================
def test_ex_ridme_type(): 
    "Check the RIDME experimental model."

    experiment = ex_ridme(tau1,tau2)

    assert isinstance(experiment,ExperimentInfo) 
# ======================================================================

# ======================================================================
def test_ex_ridme_fit(Vridme, Pr): 
    "Check the RIDME experimental model in fitting."

    experiment = ex_ridme(tau1, tau2, pathways=[1,2,3])
    Vmodel = dipolarmodel(tridme, r, Bmodel=bg_hom3d, experiment=experiment)
    result = fit(Vmodel, Vridme, ftol=1e-4)

    assert np.allclose(Vridme, result.model, rtol=1e-3) and ovl(result.P/1e5, Pr)>0.975
# ======================================================================

# ======================================================================
def test_ex_ridme_pulselength(): 
    "Check the RIDME experimental model."

    pulselength = 0.1
    experiment = ex_ridme(tau1,tau2,pathways=[1], pulselength=pulselength)
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment)

    assert np.isclose(Vmodel.reftime.ub - Vmodel.reftime.lb,6*pulselength)
# ======================================================================

# ======================================================================
def test_orisel(): 
    "Check that dipolar models with orientation selection work"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1)
    Vmodelorisel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1,orisel=lambda theta: np.ones_like(theta),gridsize=5000)

    Vref = Vmodel(mean=3,std=0.2,mod=0.3,reftime=0,conc=200,scale=1e2)
    Vorisel = Vmodelorisel(mean=3,std=0.2,mod=0.3,reftime=0,conc=200,scale=1e2)

    assert np.allclose(Vref,Vorisel,rtol=1e-4)
# ======================================================================


# ======================================================================
def test_excitationbandwidth(): 
    "Check that dipolar models with limited excitation bandwidth work"

    Vmodel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1)
    Vmodelorisel = dipolarmodel(t,r,dd_gauss,bg_hom3d,npathways=1,excbandwidth=1e8,gridsize=5000)

    Vref = Vmodel(mean=3,std=0.2,mod=0.3,reftime=0,conc=200,scale=1e2)
    Vorisel = Vmodelorisel(mean=3,std=0.2,mod=0.3,reftime=0,conc=200,scale=1e2)

    assert np.allclose(Vref,Vorisel,rtol=1e-4)
# ======================================================================

# ======================================================================
def test_twospin_modelstructure(): 
    "Check that the multi-spin model has the correct number of nonlinear parameters"

    model = dipolarmodel(t,Pmodel=dd_gauss,spins=2)

    assert model.Nnonlin==5 and model.Nlin==1 and model.Nparam==6
# ======================================================================

# ======================================================================
def test_threespin_modelstructure(): 
    "Check that the multi-spin model has the correct number of nonlinear parameters"

    model = dipolarmodel(t,spins=3,npathways=1)

    assert model.Nnonlin==13 and model.Nlin==1 and model.Nparam==14
# ======================================================================

# ======================================================================
def test_fourspin_modelstructure(): 
    "Check that the multi-spin model has the correct number of nonlinear parameters"

    triangles = [[0,1,5],[0,3,4],[1,2,4],[2,3,5]]

    model = dipolarmodel(t,spins=4,npathways=1,triangles=triangles)

    assert model.Nnonlin==31 and model.Nlin==1 and model.Nparam==32
# ======================================================================

# ======================================================================
def test_threespin_simulation(): 
    "Check that the multi-spin model runs fully without crashing"

    model = dipolarmodel(t,spins=3,npathways=1)
    
    Vsim = model(lam1=0.5, reftime1=0, lamu=0.7, conc=50, rmean1=3, rmean2=3.2, rmean3=2.8,
               chol11=0.3, chol22=0.2, chol33=0.3, chol21=0, chol31=0, chol32=0, scale=1)
    
    # Check that there are modulations
    assert not np.mean(np.diff(Vsim))>0.025
# ======================================================================

# ======================================================================
def test_fourspin_simulation(): 
    "Check that the multi-spin model runs fully without crashing"

    triangles = [[0,1,5],[0,3,4],[1,2,4],[2,3,5]]
    model = dipolarmodel(t,spins=4,npathways=1,triangles=triangles)
    Vsim = model(lam1=0.1, reftime1=0, lamu=0.3, conc=50, rmean1=3, rmean2=3, rmean3=4, rmean4=3.5, 
            rmean5=4.5, rmean6=4.5, chol11=0.2, chol22=0.2, chol33=0.2, chol44=0.2, chol55=0.2, chol66=0.2,
            chol21=0, chol31=0, chol41=0, chol51=0, chol61=0, chol32=0, chol42=0, chol52=0, chol62=0, 
            chol43=0, chol53=0, chol63=0, chol54=0, chol64=0, chol65=0, scale=1)
    
    # Check that there are modulations
    assert not np.mean(np.diff(Vsim))>0.025
# ======================================================================

# ======================================================================
def test_parametrization_reftimes_model(): 
    "Check that the parametrization leads to the correct model"

    experiment = ex_4pdeer(tau1,tau2, pathways=[1,2])
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment, parametrization='reftimes')

    assert hasattr(Vmodel,'reftime1') and hasattr(Vmodel,'reftime2')
# ======================================================================

# ======================================================================
def test_parametrization_reftimes_fit(V4pdeer,Pr): 
    "Check that the parametrization leads to the correct analysis results"

    experiment = ex_4pdeer(tau1,tau2, pathways=[1,2])
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment, parametrization='reftimes')
    result = fit(Vmodel,V4pdeer,ftol=1e-4)

    assert np.allclose(V4pdeer,result.model,rtol=1e-3) and ovl(result.P/1e5,Pr)>0.975
# ======================================================================

# ======================================================================
def test_parametrization_delays_model(): 
    "Check that the parametrization leads to the correct model"

    experiment = ex_4pdeer(tau1,tau2, pathways=[1,2])
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment, parametrization='delays')

    assert hasattr(Vmodel,'tau1') and hasattr(Vmodel,'tau2')
# ======================================================================

# ======================================================================
def test_parametrization_delays_fit(V4pdeer,Pr): 
    "Check that the parametrization leads to the correct analysis results"

    experiment = ex_4pdeer(tau1,tau2, pathways=[1,2])
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment, parametrization='delays')
    result = fit(Vmodel,V4pdeer,ftol=1e-4)

    assert np.allclose(V4pdeer,result.model,rtol=1e-3) and ovl(result.P/1e5,Pr)>0.975
# ======================================================================

# ======================================================================
def test_parametrization_shift_model(): 
    "Check that the parametrization leads to the correct model"

    experiment = ex_4pdeer(tau1,tau2, pathways=[1,2])
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment, parametrization='shift')

    assert hasattr(Vmodel,'tshift')
# ======================================================================

# ======================================================================
def test_parametrization_shift_fit(V4pdeer,Pr): 
    "Check that the parametrization leads to the correct analysis results"

    experiment = ex_4pdeer(tau1,tau2, pathways=[1,2])
    Vmodel = dipolarmodel(tdeer,r,Bmodel=bg_hom3d,experiment=experiment, parametrization='shift')
    result = fit(Vmodel,V4pdeer,ftol=1e-4)

    assert np.allclose(V4pdeer,result.model,rtol=1e-3) and ovl(result.P/1e5,Pr)>0.975
# ======================================================================

def test_docstring_dipolarmodel():
# ======================================================================
    "Check that the docstring includes all variables and keywords."
    assert_docstring(dipolarmodel)
# ======================================================================

def test_docstring_dipolarpenalty():
# ======================================================================
    "Check that the docstring includes all variables and keywords."
    assert_docstring(dipolarpenalty)
# ======================================================================

def test_docstring_ex_models():
# ======================================================================
    "Check that the docstring includes all variables and keywords."
    for model in [ex_3pdeer,ex_4pdeer,ex_fwd5pdeer,ex_rev5pdeer,ex_sifter,ex_ridme,ex_dqc]:
        assert_docstring(model)
# ======================================================================