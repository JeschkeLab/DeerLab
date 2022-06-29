
from deerlab.utils.utils import assert_docstring
import numpy as np
from deerlab import dipolarkernel, regoperator, selregparam, whitegaussnoise
from deerlab.dd_models import dd_gauss,dd_gauss2
from deerlab.utils import assert_docstring
from deerlab.solvers import cvxnnls

def test_compensate_condition():
#=======================================================================
    "Check that alpha compensates for larger condition numbers"
    
    r = np.linspace(2,6,100)
    P = dd_gauss(r,3,0.2)

    # Lower condition number    
    t1 = np.linspace(0,3,200)
    K1 = dipolarkernel(t1,r)
    V1 = K1@P
    alpha1 = selregparam(V1,K1,cvxnnls,method='aic')

    # Larger condition number
    t2 = np.linspace(0,3,400)
    K2 = dipolarkernel(t2,r)
    V2 = K2@P
    alpha2 = selregparam(V2,K2,cvxnnls,method='aic')

    assert alpha2 > alpha1
#=======================================================================

def get_alpha_from_method(method):

    t = np.linspace(0,5,500)
    r = np.linspace(2,5,80)
    P = dd_gauss(r,3.0,0.16986436005760383)
    K = dipolarkernel(t,r)
    L = regoperator(r,2,includeedges=True)
    V = K@P

    alpha = selregparam(V,K,cvxnnls,method=method,noiselvl=0,regop=L)
    return np.log10(alpha)

def test_aic_value():
#=======================================================================
    "Check that the value returned by the AIC selection method is correct"
    
    loga = get_alpha_from_method('aic')
    logaref = -6.8108 # Computed with DeerLab (0.11.0)

    assert abs(1-loga/logaref) < 0.1
#=======================================================================

def test_bic_value():
#=======================================================================
    "Check that the value returned by the BIC selection method is correct"
    
    loga = get_alpha_from_method('bic')
    logaref = -6.8089 # Computed with DeerLab (0.11.0)

    assert abs(1-loga/logaref) < 0.1
#=======================================================================

def test_aicc_value():
#=======================================================================
    "Check that the value returned by the AICc selection method is correct"
    
    loga = get_alpha_from_method('aicc')
    logaref = -6.8075 # Computed with DeerLab (0.11.0)

    assert abs(1-loga/logaref) < 0.1
#=======================================================================

def test_cv_value():
#=======================================================================
    "Check that the value returned by the CV selection method is correct"
    
    loga = get_alpha_from_method('cv')
    logaref = -5.8777 # Computed with DeerLab-Matlab (0.9.2)

    assert abs(1-loga/logaref) < 0.16
#=======================================================================

def test_gcv_value():
#=======================================================================
    "Check that the value returned by the GCV selection method is correct"
    
    loga = get_alpha_from_method('gcv')
    logaref = -5.8778 # Computed with DeerLab-Matlab (0.9.2)

    assert abs(1-loga/logaref) < 0.16
#=======================================================================

def test_rgcv_value():
#=======================================================================
    "Check that the value returned by the rGCV selection method is correct"
    
    loga = get_alpha_from_method('rgcv')
    logaref = -5.8778 # Computed with DeerLab-Matlab (0.9.2)

    assert abs(1-loga/logaref) < 0.16
#=======================================================================

def test_srgcv_value():
#=======================================================================
    "Check that the value returned by the srGCV selection method is correct"
    
    loga = get_alpha_from_method('srgcv')
    logaref = -5.8771 # Computed with DeerLab-Matlab (0.9.2)

    assert abs(1-loga/logaref) < 0.17
#=======================================================================

def test_rm_value():
#=======================================================================
    "Check that the value returned by the RM selection method is correct"
    
    loga = get_alpha_from_method('rm')
    logaref = -5.8785 # Computed with DeerLab-Matlab (0.9.2)

    assert abs(1-loga/logaref) < 0.2
#=======================================================================

def test_ee_value():
#=======================================================================
    "Check that the value returned by the EE selection method is correct"
    
    loga = get_alpha_from_method('ee')
    logaref = -5.8798 # Computed with DeerLab-Matlab (0.9.2)

    assert abs(1-loga/logaref) < 0.2
#=======================================================================

def test_ncp_value():
#=======================================================================
    "Check that the value returned by the NCP selection method is correct"
    
    loga = get_alpha_from_method('ncp')
    logaref = 1.7574 # Computed with DeerLab-Matlab (0.9.2)

    assert abs(1-loga/logaref) < 0.1
#=======================================================================

def test_mcl_value():
#=======================================================================
    "Check that the value returned by the MCL selection method is correct"
    
    loga = get_alpha_from_method('mcl')
    logaref = -5.878 # Computed with DeerLab-Matlab (0.9.2)

    assert abs(1-loga/logaref) < 0.16
#=======================================================================

def test_gml_value():
#=======================================================================
    "Check that the value returned by the GML selection method is correct"
    
    loga = get_alpha_from_method('gml')
    logaref = -7.89  # Computed with DeerLab-Matlab (0.9.2)

    assert abs(1-loga/logaref) < 0.16
#=======================================================================

def test_lr_value():
#=======================================================================
    "Check that the value returned by the LR selection method is correct"
    
    loga = get_alpha_from_method('lr')
    logaref = -7.66  # Computed with DeerLab (0.11.0)

    assert abs(1-loga/logaref) < 0.15
#=======================================================================

def test_lc_value():
#=======================================================================
    "Check that the value returned by the LC selection method is correct"
    
    loga = get_alpha_from_method('lc')
    logaref = -1.39

    assert abs(1-loga/logaref) < 0.20
#=======================================================================

def test_algorithms():
#=======================================================================
    "Check that the value returned by the the grid and Brent algorithms coincide"
    
    t = np.linspace(0,5,80)
    r = np.linspace(2,5,80)
    P = dd_gauss(r,3,0.2)
    K = dipolarkernel(t,r)
    L = regoperator(r,2)
    V = K@P + whitegaussnoise(t,0.02,seed=1)

    alpha_grid = selregparam(V,K,cvxnnls,method='aic',algorithm='grid',regop=L)
    alpha_brent = selregparam(V,K,cvxnnls,method='aic',algorithm='brent',regop=L)

    assert abs(1-alpha_grid/alpha_brent) < 0.15
#=======================================================================


def test_nonuniform_r():
#=======================================================================
    "Check the value returned when using a non-uniform distance axis"
    
    t = np.linspace(0,3,200)
    r = np.sqrt(np.linspace(1,7**2,200))
    P = dd_gauss(r,3,0.2)
    K = dipolarkernel(t,r)
    L = regoperator(r,2)
    V = K@P

    logalpha = np.log10(selregparam(V,K,cvxnnls,method='aic',regop=L))
    logalpharef = -6.8517

    assert abs(1 - logalpha/logalpharef) < 0.2 
#=======================================================================

def assert_full_output(method):

    t = np.linspace(0,5,80)
    r = np.linspace(2,5,80)
    P = dd_gauss(r,3,0.4)
    K = dipolarkernel(t,r)
    L = regoperator(r,2)
    V = K@P

    alpha,alphas_evaled,functional,residuals,penalties = selregparam(V,K,cvxnnls,method='aic',algorithm=method,full_output=True,regop=L)
    errors = []
    if np.size(alpha)!=1:
        errors.append("alphaopt is not a scalar")
    if len(functional)!=len(alphas_evaled):
        errors.append("The number of elements of functional values and evaluated alphas are different.")
    if len(residuals)!=len(penalties):
        errors.append("The number of elements of evluated residuals and penalties are different")
    if not alpha in alphas_evaled:
        errors.append("The optimal alpha is not part of the evaluated alphas")
    assert not errors, f"Errors occured:\n{chr(10).join(errors)}"

def test_full_output_brent():
#=======================================================================
    "Check that the full output argument works using the grid algorithm"

    assert_full_output('brent')
#=======================================================================

def test_full_output_grid():
#=======================================================================
    "Check that the full output argument works using the grid algorithm"

    assert_full_output('grid')
#=======================================================================

def test_unconstrained():
#=======================================================================
    "Check the algorithm works with unconstrained disributions"
    
    t = np.linspace(0,5,80)
    r = np.linspace(2,5,80)
    P = dd_gauss(r,3,0.15)
    K = dipolarkernel(t,r)
    L = regoperator(r,2,includeedges=True)
    V = K@P

    logalpha = np.log10(selregparam(V,K,np.linalg.solve,method='aic',regop=L))
    logalpharef = -8.87

    assert abs(1 - logalpha/logalpharef) < 0.1
#=======================================================================

def test_manual_candidates():
#=======================================================================
    "Check that the alpha-search range can be manually passed"
    
    t = np.linspace(0,5,80)
    r = np.linspace(2,5,80)
    P = dd_gauss(r,3,0.15)
    K = dipolarkernel(t,r)
    L = regoperator(r,2,includeedges=True)
    alphas = np.linspace(-8,2,60)
    V = K@P

    alpha_manual = np.log10(selregparam(V,K,cvxnnls,method='aic',candidates=alphas,regop=L))
    alpha_auto = np.log10(selregparam(V,K,cvxnnls,method='aic',regop=L))

    assert abs(alpha_manual-alpha_auto)<1e-4
#=======================================================================

def test_tikh_value():
#=======================================================================
    "Check that the value returned by Tikhonov regularization"
    
    np.random.seed(1)
    t = np.linspace(0,5,500)
    r = np.linspace(2,5,80)
    P = dd_gauss(r,3,0.15)
    K = dipolarkernel(t,r)
    V = K@P + whitegaussnoise(t,0.01)
    L = regoperator(r,2,includeedges=True)

    alpha = selregparam(V,K,cvxnnls,method='aic',regop=L)
    loga = np.log10(alpha)
    logaref = -3.51  # Computed with DeerLab-Matlab (0.9.2)

    assert abs(1-loga/logaref) < 0.02 # less than 2% error
#=======================================================================


def test_docstring():
# ======================================================================
    "Check that the docstring includes all variables and keywords."
    assert_docstring(selregparam)
# ======================================================================
