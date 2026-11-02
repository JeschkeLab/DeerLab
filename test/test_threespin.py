import os
import numpy as np
import pytest
from scipy.io import savemat
from scipy.stats import multivariate_normal
from deerlab.threespin import prepare_MMMx
from deerlab.utils import choleskycovmat

# Synthetic trivariate distribution in the MMMx layout: axes (r12,r13,r23) in Å
# Different axes, means and widths for each distance so that any mix-up of the axes is detected
means_A = np.array([35, 42, 30.])
cov_A = np.array([[ 9,  3, -2],
                  [ 3, 16,  5],
                  [-2,  5, 25.]])
r_axes_A = [np.arange(15,55.5,0.5), np.arange(18,66.5,0.5), np.arange(5,60.5,0.5)]

@pytest.fixture
def MMMx_file(tmp_path):
    grid = np.stack(np.meshgrid(*r_axes_A,indexing='ij'),axis=-1)
    # MMMx stores a probability density in 1/Å³
    P = multivariate_normal(means_A,cov_A).pdf(grid)
    filepath = tmp_path/'trivariate.mat'
    savemat(filepath, dict(trivariate=P, r_axis_1=r_axes_A[0], r_axis_2=r_axes_A[1], r_axis_3=r_axes_A[2]))
    return filepath

# ======================================================================
def test_prepare_MMMx_moments(MMMx_file):
    "Check that the means and covariance matrix are recovered, in nm"
    MMMx = prepare_MMMx(MMMx_file)

    assert np.allclose(MMMx['means'], means_A/10, atol=1e-4)
    assert np.allclose(MMMx['cov'], cov_A/100, atol=1e-4)
# ======================================================================

# ======================================================================
def test_prepare_MMMx_cholesky(MMMx_file):
    "Check that the Cholesky factor is lower triangular and matches the dipolarmodel parametrization"
    MMMx = prepare_MMMx(MMMx_file)
    L = MMMx['L']

    assert np.allclose(L, np.tril(L)) and np.all(np.diag(L)>0)
    assert np.allclose(L@L.T, MMMx['cov'])
    # Order of the dipolarmodel parameters: chol11, chol22, chol33, chol21, chol31, chol32
    cholfactors = [L[0,0], L[1,1], L[2,2], L[1,0], L[2,0], L[2,1]]
    assert np.allclose(choleskycovmat(3,cholfactors), MMMx['cov'])
# ======================================================================

# ======================================================================
def test_prepare_MMMx_marginals(MMMx_file):
    "Check that the distance axes and marginal distributions are consistent"
    MMMx = prepare_MMMx(MMMx_file)
    r, P_1D, P_2D = MMMx['r'], MMMx['P_1D'], MMMx['P_2D']

    for k,label in enumerate(['r12','r13','r23']):
        assert np.allclose(r[label], r_axes_A[k]/10)
        assert P_1D[label].shape == r[label].shape
        assert np.isclose(np.sum(P_1D[label]), 1)
    for pair in ['r12_r13','r12_r23','r13_r23']:
        label1, label2 = pair.split('_')
        assert P_2D[pair].shape == (len(r[label1]), len(r[label2]))
        assert np.allclose(np.sum(P_2D[pair],axis=1), P_1D[label1])
        assert np.allclose(np.sum(P_2D[pair],axis=0), P_1D[label2])
# ======================================================================

# ======================================================================
def test_prepare_MMMx_Rpo47():
    "Check the moments of an MMMx simulation of the triple-labelled Rpo47 complex (F36, F63, E123)"
    filepath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data/Rpo47_MMMx_trivariate_F36_F63_E123.mat')
    MMMx = prepare_MMMx(filepath)

    # Reference values computed independently following the MMMx conventions, in nm
    means_ref = [4.38959, 2.32370, 5.95224]
    sigmas_ref = [0.216628, 0.348898, 0.313336]
    rho_ref = [-0.043155, 0.528744, 0.345720] # (r12,r13), (r12,r23), (r13,r23)

    sigmas = np.sqrt(np.diag(MMMx['cov']))
    rho = MMMx['cov']/np.outer(sigmas,sigmas)

    assert np.allclose(MMMx['means'], means_ref, atol=1e-5)
    assert np.allclose(sigmas, sigmas_ref, atol=1e-5)
    assert np.allclose([rho[0,1], rho[0,2], rho[1,2]], rho_ref, atol=1e-5)
    assert np.allclose(MMMx['L']@MMMx['L'].T, MMMx['cov'])
# ======================================================================

# ======================================================================
def test_prepare_MMMx_invalid_file(tmp_path):
    "Check that a file without a trivariate distribution raises an error"
    filepath = tmp_path/'bivariate.mat'
    savemat(filepath, dict(bivariate=np.ones((3,3))))

    with pytest.raises(ValueError):
        prepare_MMMx(filepath)
# ======================================================================
