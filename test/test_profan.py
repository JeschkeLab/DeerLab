from numpy.random.mtrand import rand
from deerlab.classes import UQResult
import numpy as np 
from deerlab import profan, dd_gauss, whitegaussnoise, expand, link
from deerlab.utils import assert_docstring

def test_types():
# ======================================================================
    "Check that the correct data types are returned"

    r = np.linspace(2,6,300)
    sigma = 0.1
    model = dd_gauss.copy()
    model.addlinear('scale',par0 = 1)

    y = model(r,mean=3,width=0.2,scale=1) + whitegaussnoise(r,sigma,seed=1)

    profuq = profan(model,y,r,samples=3,noiselvl=sigma)

    assert isinstance(profuq['mean'],UQResult) and profuq['scale'] is None
# ======================================================================

# ======================================================================
def test_basics():
    "Check the basic functionality of the profiling"

    r = np.linspace(2,6,300)
    sigma = 0.1
    model = dd_gauss
    y = model(r,mean=3,width=0.2) + whitegaussnoise(r,sigma,seed=1)
    
    profuq = profan(model,y,r,samples=3,noiselvl=sigma)

    x,pdf = profuq['mean'].pardist()
    mean_mean = x[np.argmax(pdf)]
    x,pdf = profuq['width'].pardist()
    width_mean = x[np.argmax(pdf)]

    assert np.allclose([mean_mean,width_mean],[3,0.2],rtol=1e-2)
# ======================================================================


# ======================================================================
def test_globalmodel():
    "Check the profiling works with multiple datasets"

    r = np.linspace(2,6,300)
    sigma = 0.1
    modelA = dd_gauss.copy()
    modelB = dd_gauss.copy()
    model = expand(modelA,modelB)
    model = link(model,
            mean=[model.mean_1,model.mean_2],
            width=[model.width_1,model.width_2])
    y = model(r,r,mean=3,width=0.2,scale_1=1,scale_2=1)
    y[0] += whitegaussnoise(r,sigma,seed=1)      
    y[1] += whitegaussnoise(r,sigma,seed=1)

    profuq = profan(model,y,r,r,samples=3,noiselvl=sigma)

    x,pdf = profuq['mean'].pardist()
    mean_mean = x[np.argmax(pdf)]
    x,pdf = profuq['width'].pardist()
    width_mean = x[np.argmax(pdf)]

    assert np.allclose([mean_mean,width_mean],[3,0.2],rtol=1e-2)
# ======================================================================



# ======================================================================
def test_types():
    "Check that the correct data types are returned"

    r = np.linspace(2,6,300)
    sigma = 0.1
    model = dd_gauss.copy()
    model.addlinear('scale',par0 = 1)

    y = model(r,mean=3,width=0.2,scale=1) + whitegaussnoise(r,sigma,seed=1)

    profuq = profan(model,y,r,samples=3,noiselvl=sigma)

    assert isinstance(profuq['mean'],UQResult) and profuq['scale'] is None
# ======================================================================


# ======================================================================
def test_specific_parameters():
    "Check that parameters can be specified via keywords"

    r = np.linspace(2,6,300)
    sigma = 0.1
    model = dd_gauss.copy()

    y = model(r,mean=3,width=0.2) + whitegaussnoise(r,sigma,seed=1)

    profuq = profan(model,y,r,samples=3,noiselvl=sigma,parameters='mean')

    assert 'mean'in profuq.keys() and not 'width' in profuq.keys()
# ======================================================================


# ======================================================================
def test_docstring():
    "Check that the docstring includes all variables and keywords."
    assert_docstring(profan)
# ======================================================================
