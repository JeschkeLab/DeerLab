import numpy as np
from deerlab.solvers import fnnls, cvxnnls, _lsqcomponents
from deerlab import dipolarkernel, regoperator
from deerlab.dd_models import dd_gauss,dd_gauss2


def assert_multigauss_problem(solver):
    
    # Prepare test data
    r = np.linspace(1,8,80)
    t = np.linspace(0,4,200)
    K0 = dipolarkernel(t,r)
    par = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    r1, w1, amp1, r2, w2, amp2  = par
    # Generate basic kernel
    K0 = dipolarkernel(t,r)
    # Get Gauss basis functions
    P1 = dd_gauss(r,r1,w1)
    P2 = dd_gauss(r,r2,w2)

    amps = np.array([0.3,0.6])

    # Combine all non-linear functions into one
    A = np.zeros((len(t),2))
    A[:,0] = K0@P1
    A[:,1] = K0@P2
    V = A@amps
    KtK,KtV = _lsqcomponents(V,A)
    ampsfit = solver(KtK,KtV)
    assert max(abs(amps - ampsfit)) < 1e-8

def test_multigauss_problem_fnnls():
#=======================================================================
    "Check fnnls can solve the linear part of a multi-Gauss problem"

    assert_multigauss_problem(fnnls)
#=======================================================================

def test_multigauss_problem_cvxnnls():
#=======================================================================
    "Check cvxnnls can solve the linear part of a multi-Gauss problem"

    assert_multigauss_problem(cvxnnls)
#=======================================================================