from numpy.random.mtrand import rand
from deerlab.classes import UQResult
import numpy as np 
from deerlab import profile_analysis, dd_gauss, whitegaussnoise, merge, link
from deerlab.utils import assert_docstring
from copy import deepcopy

def test_types():
# ======================================================================
    "Check that the correct data types are returned"

    r = np.linspace(2,6,300)
    sigma = 0.1
    model = deepcopy(dd_gauss)
    model.addlinear('scale',par0 = 1)

    y = model(r,mean=3,std=0.2,scale=1) + whitegaussnoise(r,sigma,seed=1)

    profuq = profile_analysis(model,y,r,samples=3,noiselvl=sigma)

    assert isinstance(profuq['mean'],UQResult) and profuq['scale'] is None and isinstance(profuq,dict)
# ======================================================================

# ======================================================================
def test_basics():
    "Check the basic functionality of the profiling"

    r = np.linspace(2,6,300)
    sigma = 0.1
    model = dd_gauss
    y = model(r,mean=3,std=0.2) + whitegaussnoise(r,sigma,seed=1)
    
    profuq = profile_analysis(model,y,r,samples=3,noiselvl=sigma)

    x,pdf = profuq['mean'].pardist()
    mean_mean = x[np.argmax(pdf)]
    x,pdf = profuq['std'].pardist()
    std_mean = x[np.argmax(pdf)]

    assert np.allclose([mean_mean,std_mean],[3,0.2],rtol=1e-2)
# ======================================================================


# ======================================================================
def test_globalmodel():
    "Check the profiling works with multiple datasets"

    r = np.linspace(2,6,300)
    sigma = 0.1
    modelA = deepcopy(dd_gauss)
    modelB = deepcopy(dd_gauss)
    model = merge(modelA,modelB)
    model = link(model,
            mean=['mean_1','mean_2'],
            std=['std_1','std_2'])
    y = model(r,r,mean=3,std=0.2,scale_1=1,scale_2=1)
    y[0] += whitegaussnoise(r,sigma,seed=1)
    y[1] += whitegaussnoise(r,sigma,seed=1)

    profuq = profile_analysis(model,y,r,r,samples=3,noiselvl=sigma)

    x,pdf = profuq['mean'].pardist()
    mean_mean = x[np.argmax(pdf)]
    x,pdf = profuq['std'].pardist()
    width_mean = x[np.argmax(pdf)]

    assert np.allclose([mean_mean,width_mean],[3,0.2],rtol=1e-2)
# ======================================================================

# ======================================================================
def test_specific_parameters():
    "Check that parameters can be specified via keywords"

    r = np.linspace(2,6,300)
    sigma = 0.1
    model = deepcopy(dd_gauss)

    y = model(r,mean=3,std=0.2) + whitegaussnoise(r,sigma,seed=1)

    profuq = profile_analysis(model,y,r,samples=3,noiselvl=sigma,parameters='mean')

    assert 'mean'in profuq.keys() and not 'std' in profuq.keys()
# ======================================================================


# ======================================================================
def test_docstring():
    "Check that the docstring includes all variables and keywords."
    assert_docstring(profile_analysis)
# ======================================================================


# ======================================================================
def test_grids():
    "Check that grids can be specified"

    r = np.linspace(2,6,300)
    sigma = 0.1
    model = dd_gauss
    y = model(r,mean=3,std=0.2) + whitegaussnoise(r,sigma,seed=1)
    grid = {'mean':np.linspace(3,5,3)}
    profuq = profile_analysis(model,y,r,parameters='mean',grids=grid,noiselvl=sigma)

    x = profuq['mean'].profile['x']

    assert np.allclose(x,grid['mean'])
# ======================================================================

