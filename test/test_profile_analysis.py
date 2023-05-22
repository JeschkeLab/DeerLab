from numpy.random.mtrand import rand
from deerlab.classes import UQResult
import numpy as np 
from deerlab import profile_analysis, dd_gauss, whitegaussnoise, merge, link, Model
from deerlab.utils import assert_docstring
from copy import deepcopy
import pytest 

# Fixtures 
# ----------------------------------------------------------------------
x = np.linspace(2,6,300)
noiselvl = 0.1

@pytest.fixture(scope='module')
def model():
    def _model(mean,std):
        return dd_gauss(x,mean,std)
    mock_model = Model(_model)
    mock_model.mean.set(lb=dd_gauss.mean.lb,ub=dd_gauss.mean.ub,par0=dd_gauss.mean.par0)
    mock_model.std.set(lb=dd_gauss.std.lb,ub=dd_gauss.std.ub,par0=dd_gauss.std.par0)
    return mock_model 

@pytest.fixture(scope='module')
def mock_data(model):
    return model(mean=3,std=0.2) + whitegaussnoise(x,noiselvl,seed=1)
# ----------------------------------------------------------------------


def test_types(model, mock_data):
# ======================================================================
    "Check that the correct data types are returned"
    model.addlinear('scale',par0 = 1)
    profuq = profile_analysis(model,mock_data,samples=3,noiselvl=noiselvl)
    assert isinstance(profuq['mean'],UQResult) and profuq['scale'] is None and isinstance(profuq,dict)
# ======================================================================

# ======================================================================
def test_basics(model, mock_data):
    "Check the basic functionality of the profiling"
    profuq = profile_analysis(model,mock_data,samples=3,noiselvl=noiselvl)
    x,pdf = profuq['mean'].pardist()
    mean_mean = x[np.argmax(pdf)]
    x,pdf = profuq['std'].pardist()
    std_mean = x[np.argmax(pdf)]
    assert np.allclose([mean_mean,std_mean],[3,0.2],rtol=1e-2)
# ======================================================================

# ======================================================================
def test_globalmodel(model, mock_data):
    "Check the profiling works with multiple datasets"
    model = merge(model,model)
    model = link(model,
            mean=['mean_1','mean_2'],
            std=['std_1','std_2'])
    profuq = profile_analysis(model,[mock_data,mock_data],samples=3,noiselvl=noiselvl)
    x,pdf = profuq['mean'].pardist()
    mean_mean = x[np.argmax(pdf)]
    x,pdf = profuq['std'].pardist()
    width_mean = x[np.argmax(pdf)]
    assert np.allclose([mean_mean,width_mean],[3,0.2],rtol=1e-2)
# ======================================================================

# ======================================================================
def test_specific_parameters(model, mock_data):
    "Check that parameters can be specified via keywords"
    profuq = profile_analysis(model,mock_data,samples=3,noiselvl=noiselvl,parameters='mean')
    assert 'mean'in profuq.keys() and not 'std' in profuq.keys()
# ======================================================================

# ======================================================================
def test_grids(model, mock_data):
    "Check that grids can be specified"
    grid = {'mean':np.linspace(3,5,3)}
    profuq = profile_analysis(model,mock_data,parameters='mean',grids=grid,noiselvl=noiselvl)
    x = profuq['mean'].profile['x']
    assert np.allclose(x,grid['mean'])
# ======================================================================

# ======================================================================
def test_docstring():
    "Check that the docstring includes all variables and keywords."
    assert_docstring(profile_analysis)
# ======================================================================
