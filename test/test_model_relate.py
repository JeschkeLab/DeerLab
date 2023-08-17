import numpy as np
from deerlab.dd_models import dd_gauss,dd_gauss2,dd_gauss3
from deerlab.model import Model,relate
from deerlab.fit import fit
from deerlab.utils import store_pickle, read_pickle 
from copy import deepcopy
import os 
import pytest 

# Fixtures 
# ----------------------------------------------------------------------
x = np.linspace(0,10,400)

@pytest.fixture(scope='module')
def gauss():
    return deepcopy(dd_gauss)

@pytest.fixture(scope='module')
def gauss2():
    return deepcopy(dd_gauss2)

@pytest.fixture(scope='module')
def gauss3():
    return deepcopy(dd_gauss3)

# ----------------------------------------------------------------------

# ======================================================================
def test_type(gauss): 
    "Check that the function returns a valid model type"
    newmodel = relate(gauss, std = lambda mean: mean/2)
    assert isinstance(newmodel,Model)
# ======================================================================

# ======================================================================
def test_preserve_original(gauss): 
    "Check that the original input model is not changed by the function"
    _ = relate(gauss, std = lambda mean: mean/2)
    assert gauss._parameter_list() == ['mean','std']
# ======================================================================

# ======================================================================
def test_relate_signature(gauss2): 
    "Check the model signature after relate"
    newmodel = relate(gauss2,mean1= lambda mean2: mean2)
    assert newmodel.signature==['r', 'mean1', 'std1', 'std2', 'amp1', 'amp2']
# ======================================================================

# ======================================================================
def test_Nparam_nonlin(gauss): 
    "Check that the output model has the right number of parameters"
    newmodel = relate(gauss, std = lambda mean: mean/2)
    assert newmodel.Nnonlin == gauss.Nnonlin - 1
# ======================================================================

# ======================================================================
def test_Nparam_lin(gauss2): 
    "Check that the related model has the right number of parameters"
    newmodel = relate(gauss2, std1 = lambda mean1: mean1/2)
    assert newmodel.Nlin == gauss2.Nlin
# ======================================================================

# ======================================================================
def test_Nparam(gauss2): 
    "Check that the related model has the right number of parameters"
    newmodel = relate(gauss2, std1 = lambda mean1: mean1/2)
    assert newmodel.Nparam == gauss2.Nparam - 1
# ======================================================================

# ======================================================================
def test_Nparam_list(gauss3): 
    "Check that the related model has the right number of parameters"
    newmodel = relate(gauss3, std1 = lambda mean1: mean1/2)
    assert newmodel.Nparam == len(newmodel._parameter_list())
# ======================================================================

# ======================================================================
def test_param_names(gauss): 
    "Check that the related model has the adjusted parameter names"
    newmodel = relate(gauss, std = lambda mean: mean/2)
    assert all([ str in newmodel._parameter_list() for str in ['mean'] ])
# ======================================================================

# ======================================================================
def test_relate_one_nonlinear(gauss): 
    "Check that the method works correctly for one related parameter"
    ref = gauss(x,4,0.2)
    newmodel = relate(gauss, std = lambda mean: mean/20)
    response = newmodel(x,4)
    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_relate_two_nonlinear(gauss2): 
    "Check that the method works correctly for two related parameters"
    ref = gauss2(x,3,0.2,6,0.1,1,1)
    newmodel = relate(gauss2, 
                        mean1 = lambda mean2: mean2/2,
                        std1 = lambda std2: std2*2)
    response = newmodel(r=x,mean2=6,std2=0.1,amp1=1,amp2=1)
    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_relate_two_nonlinear_2(gauss2): 
    "Check that the method works correctly for two related parameters"
    ref = gauss2(x,3,0.2,3,0.1,1,1)
    newmodel = relate(gauss2, 
                        mean1 = lambda std2, std1: 10*(std1+std2),
                        mean2 = lambda std2, std1: 10*(std1+std2),)
    response = newmodel(r=x,std1=0.2,std2=0.1,amp1=1,amp2=1)
    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_relate_fit(gauss): 
    "Check that the method works correctly for one related parameter"
    ref = gauss(x,4,0.2)
    newmodel = relate(gauss, std = lambda mean: mean/20)
    result = fit(newmodel,ref,x)
    assert np.allclose(result.model,ref)
# ======================================================================

# ======================================================================
def test_relate_addednonlinear(gauss): 
    "Check that a parameter can be related with an added nonlinear parameter"
    newgauss = deepcopy(gauss)
    newgauss.addlinear('amp')
    ref = newgauss(x,4,0.2,5)
    newgauss.addnonlinear('factor')
    newmodel = relate(newgauss, std = lambda factor: factor/5)
    response = newmodel(r=x,mean=4,factor=1,amp=5)

    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_relate_conflict_order(gauss2): 
    "Check that even if the user specifies a conflicted order, it runs"
    ref = gauss2(x,4,0.2,5,0.4,1,1)
    # No conflict, there parameters are ordered correctly
    newmodel = relate(gauss2, 
                    std2=lambda std1: std1*2,
                    std1=lambda mean1: mean1/20,)
    # Conflict: std1 is deleted first, but std2 depends on it
    conflictmodel = relate(gauss2, 
                    std1=lambda mean1: mean1/20,
                    std2=lambda std1: std1*2)
    
    response1 = newmodel(r=x,mean1=4,mean2=5,amp1=1,amp2=1)
    response2 = conflictmodel(r=x,mean1=4,mean2=5,amp1=1,amp2=1)
    assert np.allclose(response1,ref) and np.allclose(response2,ref)
# ======================================================================


def test_pickle_linked_model(gauss):
# ======================================================================
    "Check that the model object can be pickled"
    ref = gauss(x,4,0.2)
    newmodel = relate(gauss, std = lambda mean: mean/20)
    store_pickle(newmodel,'pickled_model')
    try:
        pickled_model = read_pickle('pickled_model')
        result = fit(pickled_model,ref,x)

        os.remove("pickled_model.pkl") 
        assert np.allclose(result.model,ref)
    except Exception as exception:
        os.remove("pickled_model.pkl") 
        raise exception
# ======================================================================

