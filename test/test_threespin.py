import os
import numpy as np
import pytest
from scipy.io import savemat
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt
from matplotlib.table import Table
from deerlab import FitResult
from deerlab.threespin import prepare_MMMx, plot_multivariate, _extract_means_and_cholesky, _gaussian_marginals
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

# ----------------------------------------------------------------------
# plot_multivariate
# ----------------------------------------------------------------------

# Simulated trivariate Gaussian fit in nm, with the dipolarmodel parametrization (rmean{q}, chol{i}{j})
means_fit = np.array([3.5, 4.2, 3.0])
L_fit = np.linalg.cholesky(cov_A/100)

@pytest.fixture
def gauss_fit():
    "FitResult of a 3-spin dd_gauss dipolarmodel, built from simulated parameters without fitting"
    Q = len(means_fit)
    paramlist = [f'rmean{q+1}' for q in range(Q)] + [f'chol{q+1}{q+1}' for q in range(Q)]
    paramlist += [f'chol{j+1}{q+1}' for q in range(Q) for j in range(q+1,Q)]
    fit_dict = {'paramlist': paramlist}
    fit_dict.update({f'rmean{q+1}': means_fit[q] for q in range(Q)})
    fit_dict.update({f'chol{i+1}{j+1}': L_fit[i,j] for i in range(Q) for j in range(i+1)})
    return FitResult.from_dict(fit_dict)

def _table_text(fig):
    "Returns the text of all cells of the tables in a figure"
    tables = [child for ax in fig.axes for child in ax.get_children() if isinstance(child,Table)]
    return [cell.get_text().get_text() for table in tables for cell in table.get_celld().values()]

# ======================================================================
def test_extract_means_and_cholesky(gauss_fit):
    "Check that the means and Cholesky factor are recovered from a dd_gauss FitResult"
    values = _extract_means_and_cholesky(gauss_fit)

    assert values['model'] == 'dd_Gauss'
    assert np.allclose(values['means'], means_fit)
    assert np.allclose(values['L'], L_fit)
# ======================================================================

# ======================================================================
def test_extract_means_and_cholesky_invalid():
    "Check that a FitResult without multivariate distance parameters raises an error"
    fit = FitResult.from_dict({'paramlist': ['mean','std'], 'mean': 3, 'std': 0.3})

    with pytest.raises(ValueError):
        _extract_means_and_cholesky(fit)
# ======================================================================

# ======================================================================
def test_gaussian_marginals_moments():
    "Check that the marginals of a fit far from the triangle inequality boundary recover its means and covariance"
    marginals = _gaussian_marginals(means_fit, L_fit)
    r, P_1D, P_2D = marginals['r'], marginals['P_1D'], marginals['P_2D']
    Σ = L_fit@L_fit.T

    for q,label in enumerate(['r1','r2','r3']):
        assert np.isclose(np.sum(P_1D[label]), 1)
        assert np.isclose(np.sum(r[label]*P_1D[label]), means_fit[q], atol=1e-4)
        assert np.isclose(np.sum((r[label]-means_fit[q])**2*P_1D[label]), Σ[q,q], atol=1e-4)
    for i,j in [(0,1),(0,2),(1,2)]:
        Pij = P_2D[f'r{i+1}_r{j+1}']
        cov_ij = r[f'r{i+1}'] @ Pij @ r[f'r{j+1}'] - means_fit[i]*means_fit[j]
        assert np.isclose(cov_ij, Σ[i,j], atol=1e-4)
# ======================================================================

# ======================================================================
def test_gaussian_marginals_triangle():
    "Check that the marginals are truncated by the triangle inequality"
    means = np.array([2, 2, 4.])
    L = np.diag([0.3, 0.3, 0.3])
    marginals = _gaussian_marginals(means, L)
    r, P_2D = marginals['r'], marginals['P_2D']

    # Without truncation r3 would be centered at 4 nm, truncation to r3 <= r1 + r2 shifts it down
    assert np.sum(r['r3']*marginals['P_1D']['r3']) < 3.9
    # No probability where r3 >= r1 + r2 for all r1 on the grid
    violated = r['r3'][np.newaxis,:] >= r['r2'][:,np.newaxis] + r['r1'].max()
    assert np.any(violated) and np.all(P_2D['r2_r3'][violated] == 0)
# ======================================================================

# ======================================================================
def test_plot_multivariate_fit(gauss_fit):
    "Check that a 3-spin fit can be plotted with its parameter table"
    fig = plt.figure()
    plot_multivariate(Nspins=3, fit=gauss_fit, fig=fig)

    axs = {ax.get_label(): ax for ax in fig.axes}
    for label in ['r1','r2','r3']:
        assert len(axs[label].get_lines()) == 1
    for label in ['r12','r13','r23']:
        assert len(axs[label].collections) > 0
    text = _table_text(fig)
    assert 'Fit' in text
    assert all(f'{mean:.2f}' in text for mean in means_fit)
    assert all(f'{L_fit[i,j]:.3f}' in text for i in range(3) for j in range(i+1))
    plt.close(fig)
# ======================================================================

# ======================================================================
def test_plot_multivariate_fit_and_sim(gauss_fit, MMMx_file):
    "Check that a 3-spin fit can be plotted together with an MMMx simulation"
    MMMx = prepare_MMMx(MMMx_file)
    fig = plt.figure()
    plot_multivariate(Nspins=3, fit=gauss_fit, sim=MMMx, fig=fig)

    text = _table_text(fig)
    assert 'Fit' in text and 'MMMx' in text
    assert all(f'{mean:.2f}' in text for mean in MMMx['means'])
    plt.close(fig)
# ======================================================================

# ======================================================================
def test_plot_multivariate_sim(MMMx_file):
    "Check that an MMMx simulation can be plotted on its own"
    MMMx = prepare_MMMx(MMMx_file)
    fig = plt.figure()
    plot_multivariate(Nspins=3, sim=MMMx, fig=fig)

    assert 'MMMx' in _table_text(fig)
    plt.close(fig)
# ======================================================================

# ======================================================================
def test_plot_multivariate_no_table(gauss_fit):
    "Check that the parameter table is omitted when table=False"
    fig = plt.figure()
    plot_multivariate(Nspins=3, fit=gauss_fit, table=False, fig=fig)

    assert _table_text(fig) == []
    plt.close(fig)
# ======================================================================

# ======================================================================
@pytest.mark.parametrize('Nspins', [2,5])
def test_plot_multivariate_invalid_Nspins(gauss_fit, Nspins):
    "Check that unsupported numbers of spins raise an error"
    with pytest.raises(ValueError):
        plot_multivariate(Nspins=Nspins, fit=gauss_fit)
# ======================================================================
