import numpy as np
from deerlab import dipolarkernel,dd_gauss,snlls,whitegaussnoise
from deerlab.utils import skip_on, assert_docstring
import pytest

# Fixtures
# ----------------------------------------------------------------------

r = np.linspace(1,8,80)
dr = np.mean(np.diff(r))
t = np.linspace(0,4,200)
mean = 3
std = 0.2
noiselvl = 0.001

nonlin_param = np.array([0.2]) 
nlpar0 = [0.3]
lb = [0]
ub = [1]
lin_param = dd_gauss(r,mean,std)
lbl = np.zeros_like(r)
ubl = np.full_like(r, np.inf)

@pytest.fixture(scope='module')
def mock_Amodel():
    return lambda lam: dipolarkernel(t,r,mod=lam)

@pytest.fixture(scope='module')
def mock_data(mock_Amodel):
    return mock_Amodel(nonlin_param)@lin_param + whitegaussnoise(t,noiselvl,seed=1,rescale=True)


x = np.linspace(0,7,100)
comp_par0 = [2*np.pi/5,4,0.2]
comp_lb = [-np.pi,1,0.05]
comp_ub = [np.pi,6,5]

@pytest.fixture(scope='module')
def complex_mock_Amodel():
    def _model(param):
        y = dd_gauss(x,param[0], param[1])
        y = y*np.exp(-1j*param[2])
        return y
    return _model

@pytest.fixture(scope='module')
def complex_mock_data(complex_mock_Amodel):
    return complex_mock_Amodel([np.pi/5, 3, 0.5])


mask_par0 = [4,0.2]
mask_lb = [1,0.05]
mask_ub = [6,5]
 
@pytest.fixture(scope='module')
def mock_simple_model():
    return lambda p: dd_gauss(x,*p)

@pytest.fixture(scope='module')
def mock_mask():
    mask = np.ones_like(x).astype(bool)   
    mask[(x>2.5) & (x<3.5)] = False 
    return mask 

@pytest.fixture(scope='module')
def mock_uncorrupted_data(mock_simple_model):
    return mock_simple_model([3, 0.5])  

@pytest.fixture(scope='module')
def mock_corrupted_data(mock_uncorrupted_data, mock_mask):
    data = mock_uncorrupted_data.copy()
    data[~mock_mask] = 0 
    return data + whitegaussnoise(x,noiselvl,seed=1)
# ----------------------------------------------------------------------


# =======================================================================
def test_SNLLS_fit(mock_data,mock_Amodel):
    "Check the solution of SNLLS fit problems"
    fit = snlls(mock_data, mock_Amodel,nlpar0,lb,ub,lbl,ubl,uq=False)
    assert np.all(abs(lin_param - fit.lin) < 1e-1) and np.all(abs(nonlin_param - fit.nonlin) < 1e-1)
# =======================================================================

# =======================================================================
def test_SNLLS_global_fit(mock_data,mock_Amodel):
    "Check that the global SNLLS problem is solved correctly"
    fit = snlls([mock_data,mock_data],lambda p:[mock_Amodel(p[0]),mock_Amodel(p[1])],nlpar0+nlpar0,lb+lb,ub+ub,lbl,ubl, uq=False)
    assert np.all(abs(lin_param - fit.lin) < 1e-1) and np.all(abs(nonlin_param - fit.nonlin[0]) < 1e-1) and np.all(abs(nonlin_param - fit.nonlin[1]) < 1e-1)
# =======================================================================

# ======================================================================
def test_SNLLS_global_fit_with_weights(mock_data,mock_Amodel):
    "Check that the global weights in the global SNLLS fit properly work when specified"
    fit = snlls([mock_data,mock_data],lambda p: [mock_Amodel(p[0]),mock_Amodel(p[1])],nlpar0+nlpar0,lb+lb,ub+ub,lbl,ubl,uq=False,weights=[1,1e-10])
    assert np.all(abs(lin_param - fit.lin) < 1e-1) and np.all(abs(nonlin_param - fit.nonlin[0]) < 1e-1) and np.all(abs(nonlin_param - fit.nonlin[1]) < 1e-1)
    fit = snlls([mock_data,mock_data],lambda p: [mock_Amodel(p[0]),mock_Amodel(p[1])],nlpar0+nlpar0,lb+lb,ub+ub,lbl,ubl,uq=False,weights=[1e-10,1])
    assert np.all(abs(lin_param - fit.lin) < 1e-1) and np.all(abs(nonlin_param - fit.nonlin[0]) < 1e-1) and np.all(abs(nonlin_param - fit.nonlin[1]) < 1e-1)
# ======================================================================

# =======================================================================
@pytest.mark.parametrize('parameter',['lin','nonlin','model'])
def test_SNLLS_fit_confidence_intervals(mock_data,mock_Amodel,parameter):
    "Check that the SNLLS fit confidence intervals of the different quantities are correct"
    fit = snlls(mock_data, mock_Amodel, nlpar0,lb,ub,lbl,modeluq=True)
    pfit =  np.round(getattr(fit,parameter),6)
    p50lb,p50ub = np.round(getattr(fit,parameter+'Uncert').ci(50),6).T
    p95lb,p95ub = np.round(getattr(fit,parameter+'Uncert').ci(95),6).T
    if parameter == 'lin':
        bounds = [lbl,ubl]
    elif parameter == 'nonlin':
        bounds = [lb,ub]
    else: 
        bounds = [np.full_like(t, -np.inf), np.full_like(t, np.inf)]
    errors = []
    if not np.all(p95lb <= pfit) and not np.all(p50lb <= pfit):
        errors.append("Some fitted values are below the lower bound of the confidence intervals.")
    if not np.all(p95ub >= pfit) and not np.all(p50ub >= pfit):
        errors.append("Some fitted values are over the upper bound of the confidence intervals.")
    if not np.all(p95lb <= p50lb):
        errors.append("The 50%-CI has lower values than the 95%-CI")
    if not np.all(p95ub >= p50ub):
        errors.append("The 50%-CI has larger values than the 95%-CI")
    if not np.all(np.minimum(bounds[0],p95lb)==bounds[0]):
        errors.append("The lower bounds are not satisfied by the confidence intervals.")
    if not np.all(np.maximum(bounds[1],p95ub)==bounds[1]):
        errors.append("The upper bounds are not satisfied by the confidence intervals.")
    assert not errors, f"Errors occured:\n{chr(10).join(errors)}"
# =======================================================================

# ======================================================================
def test_SNLLS_fit_confidence_intervals_values():
    "Check that the values of the SNLLS fit confidence intervals are correct"
    np.random.seed(0)
    A0  = np.random.rand(300,2)
    A = lambda p: 1 - p[0] + A0**p[1]
    pnonlin = np.array([0.5,2])
    plin = np.array([0.75,5.5])
    y = A(pnonlin)@plin + 0.3*np.random.randn(300)
    # Reference confidence intervals (calculated using lmfit package, see #116)
    cov = [99.73,95.45,68.27]
    a_ci_ref = [[0.4636991758611843, 0.5411943256617782], 
                [0.47730341575508245, 0.5286896152296477], 
                [0.4905189733270308, 0.5161245529171482]]
    b_ci_ref = [[1.7876921926935445, 2.1089296764536907], 
                [1.8384482433779732, 2.0515346381926847], 
                [1.8899306688797555, 1.9961521871803736]]
    fit = snlls(y,A,[0.5,0.2],reg=False)
    a_ci = [fit.nonlinUncert.ci(cov[i])[0,:] for i in range(3)]
    b_ci = [fit.nonlinUncert.ci(cov[i])[1,:] for i in range(3)]
    ci_match = lambda ci,ci_ref,truth:np.max(abs(np.array(ci) - np.array(ci_ref)))/truth < 0.01
    assert ci_match(a_ci,a_ci_ref,pnonlin[0]) & ci_match(b_ci,b_ci_ref,pnonlin[1])
# ======================================================================

# ============================================================
def test_SNLLS_fit_confidence_intervals_upon_scaling(mock_data,mock_Amodel):
    "Check that the SNLLS fit confidence intervals are agnostic w.r.t. scaling"
    scale1, scale2 = 1, 1e5
    cis= []
    for scale in [scale1,scale2]:
        fit = snlls(mock_data*scale,mock_Amodel,nlpar0,lb,ub,lbl,ftol=1e-3)
        ci = fit.linUncert.ci(95)
        ci[ci==0] = 1e-16
        cis.append(ci)
    ci1,ci2 = cis
    assert np.all(abs(ci2/scale2 - ci1/scale1) < 1e-3)
# ============================================================

# ============================================================
def test_SNLLS_goodness_of_fit(mock_data,mock_Amodel):
    "Check the SNLLS goodness-of-fit statistics are correct"
    fit = snlls(mock_data,mock_Amodel,nlpar0,lb,ub,lbl,ubl,noiselvl=noiselvl, uq=False)
    assert abs(fit.stats['chi2red'] - 1) < 0.1
# ============================================================

# ============================================================
def test_SNLLS_goodness_of_fit_when_scaled(mock_data,mock_Amodel):
    "Check the SNLLS goodness-of-fit assesment is correct even with arbitrary scaling"
    scale = 1e5
    fit = snlls(scale*mock_data,mock_Amodel,nlpar0,lb,ub,lbl,ubl,noiselvl=scale*noiselvl, uq=False)
    assert abs(fit.stats['chi2red'] - 1) < 0.1
# ============================================================

# ============================================================
def test_SNLLS_input_shape_control(mock_data):
    "Check that the SNLLS function checks the size of the input and output"
    with pytest.raises(RuntimeError):
        snlls(mock_data,lambda _: np.zeros((len(t)+100, len(r))),nlpar0,lb,ub,lbl)
# ============================================================

# ======================================================================
@skip_on('_tkinter.TclError', reason="A problem with the Tk backend occured")
def test_plot(mock_data,mock_Amodel):
    "Check that the plot method works"
    fit = snlls(mock_data,mock_Amodel,par0=0.2,lb=0,ub=1,lbl=lbl,uq=False)
    fig = fit.plot()    
    assert str(fig.__class__)=="<class 'matplotlib.figure.Figure'>"
# ======================================================================

# ============================================================
def test_SNLLS_cost_value(mock_data,mock_Amodel):
    "Check that the SNLLS cost value is properly returned"
    fit = snlls(mock_data,mock_Amodel,nlpar0,lb,ub,lbl,ubl,uq=False)
    assert isinstance(fit.cost,float) and np.round(fit.cost/(np.sum(fit.residuals**2)/len(fit.residuals)),5)==1
#============================================================

# ======================================================================
def test_SNLLS_fit_with_extra_penalty(mock_data,mock_Amodel):
    "Check that an additional penalty can be passed correctly to the SNLLS functional"
    beta = 0.05
    compactness_penalty = lambda _,plin: beta*np.sqrt(plin*(r - np.trapz(plin*r,r))**2*dr)
    fit = snlls(mock_data,mock_Amodel,nlpar0,lb,ub,lbl,extrapenalty=compactness_penalty)
    assert np.all(abs(lin_param - fit.lin) < 1e-1) and np.all(abs(nonlin_param - fit.nonlin[0]) < 1e-1)
# ======================================================================

def test_SNLLS_fit_with_multiple_extra_penalties(mock_data,mock_Amodel):
# ======================================================================
    "Check that multiple additional penaltyies can be passed correctly to the SNLLS functional"
    beta = 0.05
    compactness_penalty = lambda _,plin: beta*np.sqrt(plin*(r - np.trapz(plin*r,r))**2*dr)
    fit = snlls(mock_data,mock_Amodel,nlpar0,lb,ub,lbl,extrapenalty=[compactness_penalty])
    R = 0.5
    radial_penalty = lambda pnonlin,_: 1/R**2*(np.linalg.norm((pnonlin-nonlin_param)/nonlin_param-R))**2
    fit_moved = snlls(mock_data,mock_Amodel,nlpar0,lb,ub,lbl,extrapenalty=[compactness_penalty,radial_penalty])
    assert np.linalg.norm(abs(lin_param - fit.lin)) < np.linalg.norm(abs(lin_param - fit_moved.lin))
# ======================================================================

def test_SNLLS_fit_with_complex_valued_model(complex_mock_data,complex_mock_Amodel):
# ======================================================================
    "Check that a complex-valued model can be fitted to complex-valued data with SNLLS"
    fitResult = snlls(complex_mock_data,complex_mock_Amodel,comp_par0,comp_lb,comp_ub)
    assert np.allclose(fitResult.model.real,complex_mock_data.real) and np.allclose(fitResult.model.imag, complex_mock_data.imag)
# ======================================================================

def test_SNLLS_uncertainty_with_complex_valued_model(complex_mock_data,complex_mock_Amodel):
# ======================================================================
    "Check the fit of a real-valued model to complex-valued data"
    complex_noisy_data = complex_mock_data + whitegaussnoise(x,0.01,seed=1) + 1j*whitegaussnoise(x,0.05,seed=2)
    fit = snlls(complex_noisy_data,complex_mock_Amodel,comp_par0,comp_lb,comp_ub,modeluq=True)
    ci_width = np.sum(fit.modelUncert.ci(95)[:,1] - fit.modelUncert.ci(95)[:,0])
    assert (ci_width.real < ci_width.imag).all()
# ======================================================================

def test_SNLSS_fit_with_masking(mock_uncorrupted_data,mock_corrupted_data,mock_mask, mock_simple_model):
# ======================================================================
    "Check that datapoints can be masked out"
    fit = snlls(mock_corrupted_data,mock_simple_model,mask_par0,mask_lb,mask_ub)
    fit_masked = snlls(mock_corrupted_data,mock_simple_model,mask_par0,mask_lb,mask_ub,mask=mock_mask)
    assert np.allclose(fit_masked.model,mock_uncorrupted_data,atol=1e-1) and not np.allclose(fit.model,mock_uncorrupted_data,atol=1e-1) 
# ======================================================================

def test_SNLLS_estimated_chi2red_when_masking(mock_corrupted_data,mock_mask, mock_simple_model):
# ======================================================================
    "Check that masking is accounted for by the goodness-of-fit"
    fit_masked = snlls(mock_corrupted_data, mock_simple_model,mask_par0,mask_lb,mask_ub,mask=mock_mask,noiselvl=noiselvl)
    fit = snlls(mock_corrupted_data, mock_simple_model,mask_par0,mask_lb,mask_ub,noiselvl=noiselvl)
    assert fit_masked.stats['chi2red']<fit.stats['chi2red'] and abs(fit_masked.stats['chi2red']-1)<0.3
# ======================================================================

# ======================================================================
def test_SNLLS_fit_with_frozen_parameters():
    "Check that linear and nonlinear parameters can be frozen during the optimization"
    def Amodel(p):
        mean1,mean2,std1,std2 = p
        return np.atleast_2d([dd_gauss.nonlinmodel(r,mean1,std1), dd_gauss.nonlinmodel(r,mean2,std2)]).T
    x = np.array([0.5,0.6])
    y = Amodel([3,5,0.2,0.3])@x
    nonlin_frozen = [None,5,None,None]
    lin_frozen = [0.5,None]
    fit = snlls(y,Amodel,par0=[3.2,5.2,0.2,0.3],lb=[0,0,0.01,0.01],ub=[10,10,5,5],lbl=[0,0])
    fit_frozen = snlls(y,Amodel,par0=[3.2,5.2,0.2,0.3],lb=[0,0,0.01,0.01],ub=[10,10,5,5],lbl=[0,0],
            nonlin_frozen=nonlin_frozen, lin_frozen=lin_frozen)
    assert np.allclose(fit_frozen.model,y,atol=1e-2) and np.allclose(fit_frozen.model,fit.model,atol=1e-2)
# ======================================================================

# ======================================================================
def test_SNLLS_with_frozen_parameters_Nparam():
    "Check that the correct number of linear and nonlinear parameters are return even when freezing"
    r = np.linspace(0,6,90)
    def Amodel(p):
        mean1,mean2,std1,std2 = p
        return np.atleast_2d([dd_gauss.nonlinmodel(r,mean1,std1), dd_gauss.nonlinmodel(r,mean2,std2)]).T
    x = np.array([0.5,0.6])
    y = Amodel([3,5,0.2,0.3])@x
    nonlin_frozen = [None,5,None,None]
    lin_frozen = [0.5,None]
    fit = snlls(y,Amodel,par0=[2,4,0.2,0.2],lb=[0,0,0.01,0.01],ub=[10,10,5,5],lbl=[0,0],
            nonlin_frozen=nonlin_frozen, lin_frozen=lin_frozen)
    fit_frozen = snlls(y,Amodel,par0=[2,4,0.2,0.2],lb=[0,0,0.01,0.01],ub=[10,10,5,5],lbl=[0,0])
    assert len(fit.nonlin)==4 and len(fit.lin)==2 and len(fit_frozen.nonlin)==4 and len(fit_frozen.lin)==2
# ======================================================================

# ======================================================================
def test_docstring():
    "Check that the docstring includes all variables and keywords."
    assert_docstring(snlls)
# ======================================================================
