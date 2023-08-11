import numpy as np
import deerlab as dl 
from deerlab.model import Model,lincombine
from deerlab.fit import fit
from deerlab.utils import store_pickle, read_pickle  
from copy import deepcopy
import os 
import pytest

# Fixtures 
# ----------------------------------------------------------------------
@pytest.fixture(scope='module')
def modelA():
    return deepcopy(dl.dd_gauss)

@pytest.fixture(scope='module')
def modelB():
    return deepcopy(dl.dd_gauss2)

@pytest.fixture(scope='module')
def modelC():
    return deepcopy(dl.dd_gauss3)

x = np.linspace(0,10,100)

@pytest.fixture(scope='module')
def model_vec():
    model = Model(lambda r: np.eye(len(r)),constants='r')
    model.addlinear('Pvec',vec=len(x),lb=0)
    return model

@pytest.fixture(scope='module')
def model_vec_normalized():
    model = Model(lambda r: np.eye(len(r)),constants='r')
    model.addlinear('Pvec',vec=len(x),lb=0,normalization=lambda Pvec: Pvec/np.sum(Pvec))
    return model 
# ----------------------------------------------------------------------

# ======================================================================
def test_type(modelA,modelB): 
    "Check that the function returns a valid model type"
    model = lincombine(modelA,modelB)
    assert isinstance(model,Model)
# ======================================================================

# ======================================================================
def test_preserve_original(modelA): 
    "Check that the original models are not changed by the function"
    _ = lincombine(modelA,modelA)
    assert modelA._parameter_list() == ['mean','std'] and modelA._parameter_list() == ['mean','std']
# ======================================================================

# ======================================================================
def test_twomodels_Nparam_nonlin(modelA,modelB): 
    "Check that the combined model has the right number of parameters"
    model = lincombine(modelA,modelB)
    assert model.Nnonlin == modelA.Nnonlin + modelB.Nnonlin
# ======================================================================

# ======================================================================
def test_twomodels_Nparam_lin(modelB,modelC): 
    "Check that the combined model has the right number of parameters"
    model = lincombine(modelB,modelC)
    assert model.Nlin == modelB.Nlin + modelC.Nlin
# ======================================================================

# ======================================================================
def test_twomodels_Nparam(modelB,modelC): 
    "Check that the combined model has the right number of parameters"
    model = lincombine(modelB,modelC)
    assert model.Nparam == modelB.Nparam + modelC.Nparam
# ======================================================================

# ======================================================================
def test_twomodels_Nparam_list(modelB,modelC): 
    "Check that the combined model has the right number of parameters"
    model = lincombine(modelB,modelC)
    assert model.Nparam == len(model._parameter_list())
# ======================================================================

# ======================================================================
def test_twomodels_param_names(modelA): 
    "Check that the combined model has the adjusted parameter names"
    model = lincombine(modelA,modelA)
    assert all([ str in model._parameter_list() for str in ['mean_1','mean_2','std_1','std_2'] ])
# ======================================================================

# ======================================================================
def test_twomodels_default_linear(modelA): 
    """Check that the default linear scaling parameter is added if there 
    are no linear parameters on one of the models"""
    model = lincombine(modelA,modelA)
    assert 'scale_2' in model._parameter_list()
# ======================================================================

# ======================================================================
def test_twomodels_addweights(modelA): 
    "Check that that weights can be introduced properly"
    model = lincombine(modelA,modelA,addweights=True)
    x = np.linspace(0,10,400)
    ref = modelA(x,3,0.2) + modelA(x,4,0.5)
    response = model(r_1=x,r_2=x,mean_1=3,std_1=0.2,
                         mean_2=4,std_2=0.5,
                         scale_1=1,scale_2=1,
                         weight_1=1,weight_2=1)
    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_twomodels_addweights_values(modelA): 
    "Check that that weights values work properly"
    model = lincombine(modelA,modelA,addweights=True)
    x = np.linspace(0,10,400)
    ref1 = modelA(x,3,0.2) 
    ref2 = modelA(x,4,0.5)
    response1 = model(r_1=x,r_2=x,mean_1=3,std_1=0.2,
                         mean_2=4,std_2=0.5,
                         scale_1=1,scale_2=1,
                         weight_1=1,weight_2=0)
    response2 = model(r_1=x,r_2=x,mean_1=3,std_1=0.2,
                         mean_2=4,std_2=0.5,
                         scale_1=1,scale_2=1,
                         weight_1=0,weight_2=1)
    assert np.allclose(response1,ref1) and np.allclose(response2,ref2)
# ======================================================================


# ======================================================================
def test_twomodels_call(modelA): 
    "Check that that lincombine works correctly for two models"
    model = lincombine(modelA,modelA)
    x = np.linspace(0,10,400)
    ref1 = modelA(x,3,0.2)
    ref2 = modelA(x,4,0.5)
    ref = ref1 + ref2
    response = model(x,x,3,0.2,4,0.5,1,1)
    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_threemodels_Nparam_nonlin(modelB,modelC): 
    "Check that the combined model has the right number of parameters"
    model = lincombine(modelB,modelB,modelC)
    assert model.Nnonlin == modelB.Nnonlin + modelB.Nnonlin + modelC.Nnonlin
# ======================================================================

# ======================================================================
def test_threemodels_Nparam_lin(modelB,modelC): 
    "Check that the combined model has the right number of parameters"
    model = lincombine(modelB,modelB,modelC)
    assert model.Nlin == modelB.Nlin + modelB.Nlin + modelC.Nlin
# ======================================================================

# ======================================================================
def test_threemodels_Nparam(modelB,modelC): 
    "Check that the combined model has the right number of parameters"
    model = lincombine(modelB,modelB,modelC)
    assert model.Nparam == modelB.Nparam + modelB.Nparam + modelC.Nparam
# ======================================================================

# ======================================================================
def test_threemodels_Nparam_list(modelB,modelC): 
    "Check that the combined model has the right number of parameters"
    model = lincombine(modelB,modelB,modelC)
    assert model.Nparam == len(model._parameter_list())
# ======================================================================

# ======================================================================
def test_threemodels_param_names(): 
    "Check that the combined model has the adjusted parameter names"
    modelA = dl.dd_gauss
    modelB = dl.dd_gauss
    modelC = dl.dd_gauss
    model = lincombine(modelA,modelB,modelC)

    assert all([ str in model._parameter_list() for str in ['mean_1','mean_2','mean_3','std_1','std_2','std_3'] ])
# ======================================================================

# ======================================================================
def test_threemodels_default_linear(modelA,modelB,modelC): 
    """Check that the default linear scaling parameter is added if there 
    are no linear parameters on one of the models"""
    model = lincombine(modelA,modelB,modelC)
    assert [scale_par in model._parameter_list() for scale_par in ['scale_1','scale_2','scale_3']]
# ======================================================================

# ======================================================================
def test_threemodels_addweights(modelA): 
    "Check that that weights can be introduced properly"
    modelB = dl.dd_rice
    modelC = dl.dd_wormchain
    model = lincombine(modelA,modelB,modelC,addweights=True)
    x = np.linspace(0,10,400)
    ref = modelA(x,3,0.2) + modelB(x,4,0.5) + modelC(x,3.7,10)

    response = model(r_1=x,r_2=x,r_3=x,
                    mean_1=3,std_1=0.2,
                    location_2=4,spread_2=0.5,
                    contour_3=3.7,persistence_3=10,
                    scale_1=1,scale_2=1,scale_3=1,
                    weight_1=1,weight_2=1,weight_3=1)

    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_threemodels_addweights_values(modelA): 
    "Check that that weights values work properly"
    modelB = dl.dd_rice
    modelC = dl.dd_gauss
    model = lincombine(modelA,modelB,modelC,addweights=True)
    x = np.linspace(0,10,400)
    ref1 = modelA(x,3,0.2) 
    ref2 = modelB(x,4,0.5)
    ref3 = modelC(x,5,0.1)

    response1 = model(r_1=x,r_2=x,r_3=x,
                         mean_1=3,std_1=0.2,
                         location_2=4,spread_2=0.5,
                         mean_3=5,std_3=0.1, 
                         scale_1=1,scale_2=1,scale_3=1,
                         weight_1=1,weight_2=0,weight_3=0)
    response2 = model(r_1=x,r_2=x,r_3=x,
                         mean_1=3,std_1=0.2,
                         location_2=4,spread_2=0.5,
                         mean_3=5,std_3=0.1, 
                         scale_1=1,scale_2=1,scale_3=1,
                         weight_1=0,weight_2=1,weight_3=0)
    response3 = model(r_1=x,r_2=x,r_3=x,
                         mean_1=3,std_1=0.2,
                         location_2=4,spread_2=0.5,
                         mean_3=5,std_3=0.1, 
                         scale_1=1,scale_2=1,scale_3=1,
                         weight_1=0,weight_2=0,weight_3=1)

    assert all([np.allclose(response,ref) for response,ref in zip([response1,response2,response3],[ref1,ref2,ref3]) ])
# ======================================================================

# ======================================================================
def test_threemodels_call(modelA): 
    "Check that that lincombine works correctly for three models"
    model = lincombine(modelA,modelA,modelA)
    x = np.linspace(0,10,400)
    ref1 = modelA(x,3,0.2)
    ref2 = modelA(x,4,0.5)
    ref3 = modelA(x,2,0.1)
    ref = ref1+ref2+ref3
    response = model(x,x,x,3,0.2,4,0.5,2,0.1,1,1,1)
    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_fit_model(modelA): 
    "Check that that lincombine works correctly for two models"
    modelB = dl.dd_rice
    model = lincombine(modelA,modelB)
    x = np.linspace(0,10,400)
    truth = modelA(x,3,0.2)+modelB(x,4,0.5)
    model.mean_1.par0=3
    model.location_2.par0=4
    result = fit(model,truth,x,x)
    assert np.allclose(result.model,truth)
# ======================================================================

# ======================================================================
def test_vec_Nparam_nonlin(modelA,model_vec): 
    "Check that the combined model with a vector-form parameter has the right number of parameters"
    model = lincombine(modelA,model_vec)
    assert model.Nnonlin == modelA.Nnonlin + model_vec.Nnonlin
# ======================================================================

# ======================================================================
def test_vec_Nparam_lin(modelB,model_vec): 
    "Check that the combined model with a vector-form parameter has the right number of parameters"
    model = lincombine(modelB,model_vec)
    assert model.Nlin == modelB.Nlin + model_vec.Nlin
# ======================================================================

# ======================================================================
def test_vec_Nparam(modelB,model_vec): 
    "Check that the combined model a vector-form parameter has the right number of parameters"
    model = lincombine(modelB,model_vec)
    assert model.Nparam == modelB.Nparam + model_vec.Nparam
# ======================================================================

# ======================================================================
def test_vec_param_names(modelA,model_vec): 
    "Check that the combined model has the adjusted parameter names"
    model = lincombine(modelA,model_vec)
    assert all([ str in model._parameter_list() for str in ['mean_1','std_1','Pvec_2'] ])
# ======================================================================

# ======================================================================
def test_vec_addweights(model_vec): 
    "Check that that weights can be introduced properly"
    model = lincombine(model_vec,model_vec,addweights=True)
    ref1 = dl.dd_gauss(x,3,0.2)
    ref2 = dl.dd_gauss(x,4,0.2)
    ref = ref1+ref2
    response = model(r_1=x,r_2=x,Pvec_1=ref1,Pvec_2=ref2,
                         weight_1=1,weight_2=1)
    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_vec_two_models(modelA,model_vec): 
    "Check that that lincombine works correctly for two models"
    model = lincombine(modelA,model_vec)
    ref1 = modelA(x,3,0.2)
    ref2 = model_vec(r=x,Pvec=modelA(x,4,0.3))
    ref = ref1 + ref2
    response = model(x,x,3,0.2,1,modelA(x,4,0.3))
    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_vec_twomodels_normalization(model_vec_normalized): 
    """Check that the normalization of linear parameter is maintained"""
    model = lincombine(model_vec_normalized,model_vec_normalized)
    ref1 = model_vec_normalized(r=x,Pvec=dl.dd_gauss(x,4,0.3))
    ref2 = model_vec_normalized(r=x,Pvec=dl.dd_gauss(x,4,0.3))
    ref = ref1 + ref2
    results = fit(model,ref,x,x)
    assert hasattr(results,'Pvec_1_scale') and hasattr(results,'Pvec_2_scale')
# ======================================================================

# ======================================================================
def test_vec_twomodels_normalization_values(modelA,model_vec_normalized): 
    """Check that the normalization of linear parameter is maintained"""
    model = lincombine(modelA,model_vec_normalized)
    model.scale_1.freeze(1)
    model.mean_1.freeze(5)
    model.std_1.freeze(0.3)
    x = np.linspace(0,10,100)
    scale2 = 79
    dist2 = dl.dd_gauss(x,5,0.3)
    ref1 = modelA(r=x,mean=2,std=0.3)
    ref2 = scale2*model_vec_normalized(r=x,Pvec=dist2/np.sum(dist2))
    ref = ref1 + ref2
    results = fit(model,ref,x,x)
    assert np.isclose(results.Pvec_2_scale,scale2)
# ======================================================================

def test_pickle_lincombined_model(modelA):
# ======================================================================
    "Check that the model object can be pickled"
    x = np.linspace(0,10,100)
    truth = modelA(x,3,0.2) + modelA(x,4,0.5)
    model = lincombine(modelA,modelA)
    model.mean_1.par0=3
    model.mean_2.par0=4
    store_pickle(model,'pickled_model')
    try:
        pickled_model =read_pickle('pickled_model')
        result = fit(pickled_model,truth,x,x)

        os.remove("pickled_model.pkl") 
        assert np.allclose(result.model,truth)
    except Exception as exception:
        os.remove("pickled_model.pkl") 
        raise exception
# ======================================================================
