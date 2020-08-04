import numpy as np
from deerlab.nnls import fnnls, cvxnnls, nnlsbpp
from deerlab import lsqcomponents, dipolarkernel, regoperator
from deerlab.dd_models import dd_gauss,dd_gauss2


def assert_multigauss_problem(solver):
    
    # Prepare test data
    r = np.linspace(1,8,80)
    t = np.linspace(0,4,200)
    K0 = dipolarkernel(t,r)
    par = [3.5, 0.4, 0.6, 4.5, 0.5, 0.4]
    P = dd_gauss2(r,par)
    np.random.seed(0)
    V = K0@P

    r1, w1, amp1, r2, w2, amp2  = par
    # Generate basic kernel
    K0 = dipolarkernel(t,r)
    # Get Gauss basis functions
    P1 = dd_gauss(r,[r1,w1])
    P2 = dd_gauss(r,[r2,w2])
    # Combine all non-linear functions into one
    K = np.zeros((len(t),2))
    K[:,0] = K0@P1
    K[:,1] = K0@P2
    L = regoperator(range(K.shape[1]),1)
    KtK,KtV = lsqcomponents(V,K,L,0,1)
    amps = solver(KtK,KtV)
    assert max(abs(amps - [amp1, amp2])) < 1e-8


def test_multigauss_problem_fnnls():
#=======================================================================
    "Check fnnls can solve the linear part of a multi-Gauss problem"

    assert_multigauss_problem(fnnls)
#=======================================================================

def test_multigauss_problem_nnlsbpp():
#=======================================================================
    "Check nnlsbpp can solve the linear part of a multi-Gauss problem"

    assert_multigauss_problem(nnlsbpp)
#=======================================================================

def test_multigauss_problem_cvxnnls():
#=======================================================================
    "Check nnlsbpp can solve the linear part of a multi-Gauss problem"

    assert_multigauss_problem(cvxnnls)
#=======================================================================