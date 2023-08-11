import numpy as np
import deerlab as dl 
from deerlab.model import link,Model
from deerlab.fit import fit
from deerlab.utils import store_pickle, read_pickle 
from copy import deepcopy
import os  
import pytest 

# Fixtures
# ----------------------------------------------------------------------
x = np.linspace(0,10,80)

@pytest.fixture(scope='module')
def gauss1():
    return deepcopy(dl.dd_gauss)

@pytest.fixture(scope='module')
def gauss2():
    return deepcopy(dl.dd_gauss2)

@pytest.fixture(scope='module')
def gauss3():
    return deepcopy(dl.dd_gauss3)

@pytest.fixture(scope='module')
def model_double_vec():
    model =  Model(lambda shift, scale: shift + scale*np.eye(80))
    model.addlinear('vec1',vec=40,lb=0,par0=0)
    model.addlinear('vec2',vec=40,lb=0,par0=0)
    return model 
# ----------------------------------------------------------------------


# ======================================================================
def test_link_name(gauss2): 
    "Check that the parameters are linked to the proper name"
    linkedmodel = link(gauss2,mean=['mean1','mean2'])
    assert hasattr(linkedmodel,'mean') and not hasattr(linkedmodel,'mean1') and not hasattr(linkedmodel,'mean2')
# ======================================================================

# ======================================================================
def test_link_Nparam(gauss2): 
    "Check that the parameters are linked resulting in proper number of parameters"
    linkedmodel = link(gauss2,mean=['mean1','mean2'])

    assert linkedmodel.Nparam == gauss2.Nparam - 1
# ======================================================================

# ======================================================================
def test_link_Nparam_list(gauss2): 
    "Check that the combined model has the right number of parameters"
    linkedmodel = link(gauss2,mean=['mean1','mean2'])
    assert linkedmodel.Nparam == len(linkedmodel._parameter_list())
# ======================================================================

# ======================================================================
def test_link_signature(gauss2): 
    "Check the model signature after linking"
    linkedmodel = link(gauss2,mean=['mean1','mean2'])
    assert linkedmodel.signature==['r', 'mean', 'std1', 'std2', 'amp1', 'amp2']
# ======================================================================

# ======================================================================
def test_link_multiple(gauss3): 
    "Check that multiple links can be established"
    linkedmodel = link(gauss3,
            mean1=['mean1','mean2'],
            std2=['std2','std3'])
    assert linkedmodel.Nparam == gauss3.Nparam - 2
# ======================================================================

# ======================================================================
def test_link_call(gauss2): 
    "Check that linked parameter models return the correct responses"
    linkedmodel = link(gauss2,mean=['mean1','mean2'])
    ref = gauss2(x,4,0.5,4,0.2,1,1)
    response = linkedmodel(x,4,0.5,0.2,1,1)
    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_link_fit(gauss2): 
    "Check that linked parameter models can be properly fitted"
    linkedmodel = link(gauss2,
                mean=['mean1','mean2'],
                std=['std1','std2'])
    linkedmodel.mean.par0 = 3
    linkedmodel.std.par0 = 0.5
    ref = gauss2(x,3,0.5,3,0.5,1,1)
    result = fit(gauss2,ref,x)
    assert np.allclose(result.model,ref,atol=1e-5)
# ======================================================================

# ======================================================================
def test_link_linear_name(gauss3): 
    "Check that the parameters are linked to the proper name"
    linkedmodel = link(gauss3,amp=['amp1','amp3'])
    assert hasattr(linkedmodel,'amp') and not hasattr(linkedmodel,'amp1') and not hasattr(linkedmodel,'amp3')
# ======================================================================

# ======================================================================
def test_link_linear_order(gauss3): 
    "Check that the parameters are linked to the proper name"
    linkedmodel1 = link(gauss3,amp=['amp1','amp3'])
    linkedmodel2 = link(gauss3,amp=['amp3','amp1'])
    for linkedmodel in [linkedmodel1,linkedmodel2]:
        assert hasattr(linkedmodel,'amp') and not hasattr(linkedmodel,'amp1') and not hasattr(linkedmodel,'amp3')
# ======================================================================

# ======================================================================
def test_link_linear_call(gauss2): 
    "Check that linked parameter models return the correct responses"
    linkedmodel = link(gauss2,amp=['amp1','amp2'])
    ref = gauss2(x,4,0.5,2,0.2,0.9,0.9)
    response = linkedmodel(x,4,0.5,2,0.2,0.9)
    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_vec_link_name(model_double_vec): 
    "Check that the parameters are linked to the proper name"
    linkedmodel = link(model_double_vec,vec=['vec1','vec2'])
    assert hasattr(linkedmodel,'vec') and not hasattr(linkedmodel,'vec1') and not hasattr(linkedmodel,'vec2')
# ======================================================================

# ======================================================================
def test_vec_link_Nparam(model_double_vec): 
    "Check that the parameters are linked resulting in proper number of parameters"
    linkedmodel = link(model_double_vec,vec=['vec1','vec2'])
    assert linkedmodel.Nparam == model_double_vec.Nparam - 40
# ======================================================================

# ======================================================================
def test_vec_link_call(model_double_vec,gauss1): 
    "Check that linked parameter models return the correct responses"
    linkedmodel = link(model_double_vec,vec=['vec1','vec2'])
    x = np.linspace(0,10,40)
    ref = model_double_vec(1,2,gauss1(x,3,0.2),gauss1(x,3,0.2))
    response = linkedmodel(1,2,gauss1(x,3,0.2))
    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_vec_link_fit(model_double_vec,gauss1): 
    "Check that linked parameter models can be properly fitted"
    linkedmodel = link(model_double_vec,vec=['vec1','vec2'])
    linkedmodel.shift.par0 = 1
    linkedmodel.scale.par0 = 2
    linkedmodel.shift.freeze(1)
    linkedmodel.scale.freeze(2)
    x = np.linspace(0,10,40)
    ref = model_double_vec(1,2,gauss1(x,3,0.6),gauss1(x,3,0.6))
    result = fit(linkedmodel,ref,ftol=1e-3)
    assert np.allclose(result.model,ref,atol=1e-2)
# ======================================================================

def test_pickle_linked_model(gauss2):
# ======================================================================
    "Check that the model object can be pickled"
    ref = gauss2(x,3,0.5,3,0.5,1,1)
    linkedmodel = link(gauss2,
                mean=['mean1','mean2'],
                std=['std1','std2'])
    linkedmodel.mean.par0 = 3
    linkedmodel.std.par0 = 0.5
    store_pickle(linkedmodel,'pickled_model')
    try:
        pickled_model = read_pickle('pickled_model')
        os.remove("pickled_model.pkl") 
        assert pickled_model.signature == ['r', 'mean', 'std', 'amp1', 'amp2']
    except Exception as exception:
        os.remove("pickled_model.pkl") 
        raise exception
# ======================================================================
