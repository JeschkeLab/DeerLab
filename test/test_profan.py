import numpy as np 
from deerlab import profan, dd_gauss, whitegaussnoise
from deerlab.utils import assert_docstring

def test_nonlinearmodel():
# ======================================================================
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



def test_docstring():
# ======================================================================
    "Check that the docstring includes all variables and keywords."
    assert_docstring(profan)
# ======================================================================
