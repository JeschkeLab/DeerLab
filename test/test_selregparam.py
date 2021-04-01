
import numpy as np
from deerlab import dipolarkernel, regoperator, regparamrange, selregparam, whitegaussnoise
from deerlab.dd_models import dd_gauss,dd_gauss2

def test_compensate_condition():
#=======================================================================
    "Check that alpha compensates for larger condition numbers"
    
    r = np.linspace(2,6,100)
    P = dd_gauss(r,[3,0.2])

    # Lower condition number    
    t1 = np.linspace(0,3,200)
    K1 = dipolarkernel(t1,r)
    V1 = K1@P
    alpha1 = selregparam(V1,K1,r,'tikhonov','aic')

    # Larger condition number
    t2 = np.linspace(0,3,400)
    K2 = dipolarkernel(t2,r)
    V2 = K2@P
    alpha2 = selregparam(V2,K2,r,'tikhonov','aic')

    assert alpha2 > alpha1
#=======================================================================

def get_alpha_from_method(method):

    t = np.linspace(0,5,500)
    r = np.linspace(2,5,80)
    P = dd_gauss(r,[3,0.16986436005760383])
    K = dipolarkernel(t,r)
    V = K@P

    alpha = selregparam(V,K,r,'tikhonov',method,noiselvl=0)
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
    logaref = -1.7574 # Computed with DeerLab-Matlab (0.9.2)

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
    logaref = -0.50  # Computed with DeerLab (0.11.0)

    assert abs(1-loga/logaref) < 0.15
#=======================================================================

def test_lc_value():
#=======================================================================
    "Check that the value returned by the LC selection method is correct"
    
    loga = get_alpha_from_method('lc')
    logaref = -7.900  # Computed with DeerLab-Matlab (0.9.2)

    assert abs(1-loga/logaref) < 0.20
#=======================================================================

def test_algorithms():
#=======================================================================
    "Check that the value returned by the the grid and Brent algorithms coincide"
    
    t = np.linspace(0,5,80)
    r = np.linspace(2,5,80)
    P = dd_gauss(r,[3,0.16986436005760383])
    K = dipolarkernel(t,r)
    V = K@P

    alpha_grid = selregparam(V,K,r,'tikhonov','aic',algorithm='grid')
    alpha_brent = selregparam(V,K,r,'tikhonov','aic',algorithm='brent')

    assert abs(1-alpha_grid/alpha_brent) < 0.1
#=======================================================================


def test_nonuniform_r():
#=======================================================================
    "Check the value returned when using a non-uniform distance axis"
    
    t = np.linspace(0,3,200)
    r = np.sqrt(np.linspace(1,7**2,200))
    P = dd_gauss(r,[3,0.2])
    K = dipolarkernel(t,r)
    V = K@P

    logalpha = np.log10(selregparam(V,K,r,'tikhonov','aic'))
    logalpharef = -6.8517

    assert abs(1 - logalpha/logalpharef) < 0.2 
#=======================================================================


def test_tikh_global():
#=======================================================================
    "Check the value returned when using global Tikhonov regularization"

    t1 = np.linspace(0,4,50)
    t2 = np.linspace(0,3,80)
    t3 = np.linspace(0,4,70)

    r = np.linspace(2,5,80)
    P = dd_gauss2(r,[3,0.15,0.3,3.5,0.15,0.7])

    K1 = dipolarkernel(t1,r)
    S1 = K1@P + whitegaussnoise(t1,0.03)
    K2 = dipolarkernel(t2,r)
    S2 = K2@P + whitegaussnoise(t2,0.02)
    K3 = dipolarkernel(t3,r)
    S3 = K3@P + whitegaussnoise(t3,0.02)

    logalpha = np.log10(selregparam([S1,S2,S3],[K1,K2,K3],r,'tikhonov','aic',weights=[1,2,2]))
    logalpharef = -3.273

    assert abs(1 - logalpha/logalpharef) < 0.1

def test_tv_global():
#=======================================================================
    "Check the value returned when using global TV regularization"

    t1 = np.linspace(0,4,50)
    t2 = np.linspace(0,3,80)
    t3 = np.linspace(0,4,70)

    r = np.linspace(2,5,80)
    P = dd_gauss2(r,[3,0.15,0.3,3.5,0.15,0.7])

    K1 = dipolarkernel(t1,r)
    S1 = K1@P + whitegaussnoise(t1,0.03)
    K2 = dipolarkernel(t2,r)
    S2 = K2@P + whitegaussnoise(t2,0.02)
    K3 = dipolarkernel(t3,r)
    S3 = K3@P + whitegaussnoise(t3,0.02)

    logalpha = np.log10(selregparam([S1,S2,S3],[K1,K2,K3],r,'tv','aic',weights=[1,2,2]))
    logalpharef = -4.574

    assert abs(1 - logalpha/logalpharef) < 0.1

def test_huber_global():
#=======================================================================
    "Check the value returned when using global Huber regularization"

    t1 = np.linspace(0,4,50)
    t2 = np.linspace(0,3,80)
    t3 = np.linspace(0,4,70)

    r = np.linspace(2,5,80)
    P = dd_gauss2(r,[3,0.15,0.3,3.5,0.15,0.7])

    K1 = dipolarkernel(t1,r)
    S1 = K1@P + whitegaussnoise(t1,0.03)
    K2 = dipolarkernel(t2,r)
    S2 = K2@P + whitegaussnoise(t2,0.02)
    K3 = dipolarkernel(t3,r)
    S3 = K3@P + whitegaussnoise(t3,0.02)

    logalpha = np.log10(selregparam([S1,S2,S3],[K1,K2,K3],r,'huber','aic',weights=[1,2,2]))
    logalpharef = -3.114

    assert abs(1 - logalpha/logalpharef) < 0.1


def assert_full_output(method):

    t = np.linspace(0,5,80)
    r = np.linspace(2,5,80)
    P = dd_gauss(r,[3,0.4])
    K = dipolarkernel(t,r)
    V = K@P

    alpha,alphas_evaled,functional,residuals,penalties = selregparam(V,K,r,'tikhonov','aic',algorithm=method,full_output=True)
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

def test_full_output_grid():
#=======================================================================
    "Check that the full output argument works using the grid algorithm"

    assert_full_output('grid')


def test_unconstrained():
#=======================================================================
    "Check the algorithm works with unconstrained disributions"
    
    t = np.linspace(0,5,80)
    r = np.linspace(2,5,80)
    P = dd_gauss(r,[3,0.15])
    K = dipolarkernel(t,r)
    V = K@P

    logalpha = np.log10(selregparam(V,K,r,'tikhonov','aic',nonnegativity=False))
    logalpharef = -8.87

    assert abs(1 - logalpha/logalpharef) < 0.1
#=======================================================================

def test_manual_candidates():
#=======================================================================
    "Check that the alpha-search range can be manually passed"
    
    t = np.linspace(0,5,80)
    r = np.linspace(2,5,80)
    P = dd_gauss(r,[3,0.15])
    K = dipolarkernel(t,r)
    L = regoperator(r,2)
    alphas = regparamrange(K,L)
    V = K@P

    alpha_manual = np.log10(selregparam(V,K,r,'tikhonov','aic',candidates=alphas))
    alpha_auto = np.log10(selregparam(V,K,r,'tikhonov','aic'))

    assert abs(alpha_manual-alpha_auto)<1e-4
#=======================================================================

def get_alpha_from_regtype(regtype):

    np.random.seed(1)
    t = np.linspace(0,5,500)
    r = np.linspace(2,5,80)
    P = dd_gauss(r,[3,0.15])
    K = dipolarkernel(t,r)
    V = K@P + whitegaussnoise(t,0.01)

    alpha = selregparam(V,K,r,regtype,'aic')
    return np.log10(alpha)

def test_tikh_value():
#=======================================================================
    "Check that the value returned by Tikhonov regularization"
    
    loga = get_alpha_from_regtype('tikhonov')
    logaref = -3.51  # Computed with DeerLab-Matlab (0.9.2)

    assert abs(1-loga/logaref) < 0.02 # less than 2% error
#=======================================================================

def test_tv_value():
#=======================================================================
    "Check that the value returned by TV regularization"
    
    loga = get_alpha_from_regtype('tv')
    logaref = -4.3632  # Computed with DeerLab (0.11.0)

    assert abs(1-loga/logaref) < 0.02 # less than 2% error
#=======================================================================

def test_huber_value():
#=======================================================================
    "Check that the value returned by Huber regularization"
    
    loga = get_alpha_from_regtype('huber')
    logaref = -3.27  # Computed with DeerLab-Matlab (0.9.2)

    assert abs(1-loga/logaref) < 0.02 # less than 2% error
#=======================================================================