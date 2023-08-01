
import numpy as np
from deerlab import dipolarkernel, regoperator, whitegaussnoise, snlls
from deerlab.dd_models import dd_gauss,dd_gauss2
from deerlab.utils import ovl, skip_on, assert_docstring
import pytest 

# Fixtures 
# ------------------------------------------------------------
t = np.linspace(-1,6,300)
r = np.linspace(2,4,100)
lbl = np.zeros_like(r)
noiselvl = 0.005

@pytest.fixture(scope='module')
def reference():
    return dd_gauss(r,3,0.2)

@pytest.fixture(scope='module')
def reg_matrix():
    return regoperator(np.arange(len(r)),2)

@pytest.fixture(scope='module')
def design_matrix():
    return dipolarkernel(t,r)

@pytest.fixture(scope='module')
def mock_data(design_matrix, reference):
    return design_matrix@reference + whitegaussnoise(t,noiselvl,seed=1,rescale=True)
# ------------------------------------------------------------


# ============================================================

try :
    import quadprog
    solvers = ['qp','fnnls','cvx']
except ImportError:
    solvers = ['fnnls','cvx']

@pytest.mark.parametrize('solver',solvers)
def test_RLLS_fit_solvers(mock_data, design_matrix, reference, solver):
    "Check that the RLLS problem is correctly solved by all numerical solvers"
    fit = snlls(mock_data,design_matrix,lbl=lbl,nnlsSolver=solver,uq=False)
    assert ovl(reference,fit.param) > 0.95 # more than 95% overlap
# ============================================================

# ============================================================
def test_RLLS_fit_unconstrained(mock_data, design_matrix, reg_matrix):
    "Check the solution of the RLLS unconstrained problem"
    fit = snlls(mock_data, design_matrix, regparam=0.2)
    Pfit_ref = np.linalg.solve(design_matrix.T@design_matrix + 0.2**2*reg_matrix.T@reg_matrix, design_matrix.T@mock_data)
    assert ovl(Pfit_ref,fit.param) > 0.95  
# ============================================================

# ============================================================
def test_RLLS_global_fit(mock_data, design_matrix, reference):
    "Check the solution of a global RLLS problem"
    fit = snlls([mock_data]*2,[design_matrix]*2,lbl=lbl)
    assert ovl(reference,fit.param) > 0.95  
# ============================================================

# ============================================================
def test_global_weights(mock_data, design_matrix, reference):
    "Check that the global weights properly work when specified"
    fit1 = snlls([mock_data]*2,[design_matrix]*2,lbl=np.zeros_like(r),weights=[1,0])
    fit2 = snlls([mock_data]*2,[design_matrix]*2,lbl=np.zeros_like(r),weights=[0,1])
    assert ovl(reference,fit1.param) > 0.95 and ovl(reference,fit2.param) > 0.95
# ============================================================

def test_RLLS_fit_confidence_interval_values():
# ============================================================
    "Check that the values of the confidence intervals are correct"
    np.random.seed(0)
    A = np.random.rand(300,2)
    p = np.array([0.5,3])
    y = A@p + 0.2*np.random.randn(300)
    # Reference confidence intervals (calculated using lmfit package, see #116)
    cov = [99.73,95.45,68.27]
    a_ci_ref = [[0.41918705096124576, 0.6052333453428683],
                [0.4504549830560608, 0.5739654132472912], 
                [0.48140740318342196, 0.5430129931207299]]
    b_ci_ref = [[2.8472129716943995, 3.0317740823087562], 
                [2.878254646104185, 3.00073387702068], 
                [2.9089331341860465, 2.9700584546043713]]
    fit = snlls(y,A,reg=False)
    a_ci = [fit.paramUncert.ci(cov[i])[0,:] for i in range(3)]
    b_ci = [fit.paramUncert.ci(cov[i])[1,:] for i in range(3)]
    ci_match = lambda ci,ci_ref,truth:np.max(abs(np.array(ci) - np.array(ci_ref)))/truth < 0.05
    assert ci_match(a_ci,a_ci_ref,p[0]) & ci_match(b_ci,b_ci_ref,p[1])
# ============================================================

def test_RLLS_fit_whether_scale_agnostic(mock_data, design_matrix):
#============================================================
    "Check agnosticity of the w.r.t. an arbitrary scale of the data"
    scale = 1e8
    fit = snlls(scale*mock_data,design_matrix,lbl=lbl)
    assert max(abs(scale*mock_data - fit.model)/scale) < 0.1
#============================================================

# ============================================================
def test_RLLS_goodness_of_fit(mock_data, design_matrix):
    "Check the goodness-of-fit statistics are correct"
    fit = snlls(mock_data,design_matrix,lbl=lbl,noiselvl=noiselvl)
    assert abs(fit.stats['chi2red'] - 1) < 0.2
# ============================================================

# ============================================================
def test_RLLS_goodness_of_fit_when_scaled(mock_data, design_matrix):
    "Check the goodness-of-fit statistics are correct even with arbitrary scales"
    scale = 1e5
    fit = snlls(scale*mock_data,scale*design_matrix,lbl=np.zeros_like(r),noiselvl=scale*noiselvl)
    assert abs(fit.stats['chi2red'] - 1) < 0.2
# ============================================================

# ============================================================
def test_RLLS_global_goodness_of_fit(mock_data, design_matrix):
    "Check the goodness-of-fit statistics are correct with multiple datasets"
    fit = snlls([mock_data]*2,[design_matrix]*2,lbl=lbl,noiselvl=[noiselvl]*2)
    assert (abs(fit.stats[0]['chi2red'] - 1) < 0.1) and (abs(fit.stats[1]['chi2red'] - 1) < 0.1)
# ============================================================

# ============================================================
@skip_on('_tkinter.TclError', reason="A problem with the Tk backend occured")
def test_plot(mock_data, design_matrix):
    "Check the plotting method"
    fit = snlls(mock_data,design_matrix,lbl=lbl)
    fig = fit.plot()
    assert str(fig.__class__)=="<class 'matplotlib.figure.Figure'>"
# ============================================================

#============================================================
def test_RLLS_fit_cost_value(mock_data, design_matrix):
    "Check that the cost value is properly returned"
    fit = snlls(mock_data,design_matrix,lbl=np.zeros_like(r))
    assert isinstance(fit.cost,float) and np.round(fit.cost/(np.sum(fit.residuals**2)/len(fit.residuals)),5)==1
#============================================================

#============================================================
def test_RLLS_fit_convergence_criteria(mock_data, design_matrix, reference):
    "Check that convergence criteria can be specified without crashing"
    fit = snlls(mock_data,design_matrix,lin_tol=1e-9,lin_maxiter=2e3,lbl=np.zeros_like(r))
    assert ovl(reference,fit.param) > 0.90 # more than 80% overlap
#============================================================

# ======================================================================
def test_frozen_values():
    "Check that linear parameters can be frozen during the optimization"
    x = np.linspace(0,6,100)
    def gauss(mean,std): 
        return np.exp(-(x-mean)**2/std**2/2)
    A = np.squeeze(np.atleast_2d([gauss(3,0.4), gauss(4,0.2)]).T)
    y = A@np.array([0.5,0.6])
    # Freeze the first linear parameters
    xfrozen = [0.5,None]
    fit = snlls(y,A,lbl=np.zeros(2),lin_frozen=xfrozen)
    assert np.allclose(fit.model,y)
# ======================================================================

# ======================================================================
def test_frozen_Nparam():
    "Check that the correct number of parameters are returned even with frozen parameters"
    x = np.linspace(0,6,100)
    def gauss(mean,std): 
        return np.exp(-(x-mean)**2/std**2/2)
    A = np.squeeze(np.atleast_2d([gauss(3,0.4), gauss(4,0.2)]).T)
    y = A@np.array([0.5,0.6])
    xfrozen = [0.5,None]
    fit = snlls(y,A,lbl=np.zeros(2))
    fit_frozen = snlls(y,A,lbl=np.zeros(2),lin_frozen=xfrozen)
    assert len(fit.param)==2 and len(fit_frozen.param)==2
# ======================================================================