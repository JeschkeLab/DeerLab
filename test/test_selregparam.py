
from deerlab.utils import assert_docstring
import numpy as np
from deerlab import dipolarkernel, regoperator, selregparam, whitegaussnoise
from deerlab.dd_models import dd_gauss,dd_gauss2
from deerlab.utils import assert_docstring
from deerlab.solvers import cvxnnls
import pytest

# Fixtures 
# ------------------------------------------------------------------------
t = np.linspace(0,5,100)
r = np.linspace(2,5,80)

@pytest.fixture(scope='module')
def design_matrix():
    return dipolarkernel(t,r)
    
@pytest.fixture(scope='module')
def regularization_matrix():
    return regoperator(r,2,includeedges=True)

@pytest.fixture(scope='module')
def dataset(design_matrix):
    P = dd_gauss(r,3,0.16986436005760383)
    return design_matrix@P
# ------------------------------------------------------------------------

#=======================================================================
# References computed with DeerLab (1.0) grid method using quadprog
@pytest.mark.parametrize('criterion, log_reference', [
    ('aic', -8.37288),
    ('bic', -8.37288),
    ('aicc', -8.37288),
    ('cv', -8.37288),
    ('gcv', -8.37288),
    ('rgcv', -8.37288),
    ('srgcv', -8.16949),
    ('rm', -8.28318),
    ('ee', -8.37288),
    ('ncp', -9.18644),
    ('mcl', -8.37288),
    ('gml', -9.2670),
    ('lr', -8.37288),
    ('lc', -8.37288),
])
def test_selection_criteria_values(dataset, design_matrix, regularization_matrix, criterion,log_reference):
    "Check that the value returned by the selection criteria are correct"
    alpha = selregparam(dataset, design_matrix, cvxnnls, method=criterion, searchrange=[1e-10,1e-6], noiselvl=0, regop=regularization_matrix)
    assert abs(1-np.log10(alpha)/log_reference) < 0.25
#=======================================================================

#=======================================================================
def test_algorithms(dataset, design_matrix, regularization_matrix):
    "Check that the value returned by the the grid and Brent algorithms coincide"
    alpha_grid = selregparam(dataset,design_matrix,cvxnnls,method='aic',algorithm='grid',regop=regularization_matrix, searchrange=[1e-9,1e-5])
    alpha_brent = selregparam(dataset,design_matrix,cvxnnls,method='aic',algorithm='brent',regop=regularization_matrix, searchrange=[1e-9,1e-5])

    assert abs(1-alpha_grid/alpha_brent) < 0.2
#=======================================================================

#=======================================================================
def test_nonuniform_r(dataset):
    "Check the value returned when using a non-uniform distance axis"
    r = np.sqrt(np.linspace(1,7**2,200))
    K = dipolarkernel(t,r)
    L = regoperator(r,2)

    logalpha = np.log10(selregparam(dataset,K,cvxnnls,method='aic',regop=L))
    logalpharef = -4.3325

    assert abs(1 - logalpha/logalpharef) < 0.2 
#=======================================================================

#=======================================================================
@pytest.mark.parametrize('algorithm',['brent','grid'])
def test_full_output_behavior_with_algorithms(dataset,design_matrix,regularization_matrix,algorithm):
    "Check that the full output argument works using all algorithms"
    alpha,alphas_evaled,functional,residuals,penalties = selregparam(dataset,design_matrix,cvxnnls,method='aic',algorithm=algorithm,full_output=True,regop=regularization_matrix)

    errors = []
    if np.size(alpha)!=1:
        errors.append("alphaopt is not a scalar")
    if len(functional)!=len(alphas_evaled):
        errors.append("The number of elements of functional values and evaluated alphas are different.")
    if len(residuals)!=len(penalties):
        errors.append("The number of elements of evaluated residuals and penalties are different")
    if not alpha in alphas_evaled:
        errors.append("The optimal alpha is not part of the evaluated alphas")
    assert not errors, f"Errors occured:\n{chr(10).join(errors)}"
#=======================================================================

#=======================================================================
def test_unconstrained(dataset,design_matrix,regularization_matrix):
    "Check the algorithm works with unconstrained disributions"
    logalpha = np.log10(selregparam(dataset,design_matrix,np.linalg.solve,method='aic',regop=regularization_matrix))
    logalpharef = -8.87

    assert abs(1 - logalpha/logalpharef) < 0.1
#=======================================================================

#=======================================================================
def test_manual_candidates(dataset,design_matrix,regularization_matrix):
    "Check that the alpha-search range can be manually passed"
    alphas = np.linspace(-8,2,60)
    alpha_manual = np.log10(selregparam(dataset,design_matrix,cvxnnls,method='aic',candidates=alphas,regop=regularization_matrix))
    alpha_auto = np.log10(selregparam(dataset,design_matrix,cvxnnls,method='aic',regop=regularization_matrix))
    assert abs(alpha_manual-alpha_auto)<1e-4
#=======================================================================

#=======================================================================
def test_tikh_value(dataset,design_matrix,regularization_matrix):
    "Check that the value returned by Tikhonov regularization"
    alpha = selregparam(dataset + whitegaussnoise(t,0.01,seed=1),design_matrix,cvxnnls,method='aic',regop=regularization_matrix)
    loga = np.log10(alpha)
    logaref = -3.667  # Computed with DeerLab-Matlab (1.0.0)
    assert abs(1-loga/logaref) < 0.02 # less than 2% error
#=======================================================================

# ======================================================================
def test_docstring():
    "Check that the docstring includes all variables and keywords."
    assert_docstring(selregparam)
# ======================================================================
