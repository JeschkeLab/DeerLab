import numpy as np
import deerlab as dl 
from deerlab.model import Model,fit,lincombine
from deerlab.utils import store_pickle, read_pickle  
import os 

# ======================================================================
def test_type(): 
    "Check that the function returns a valid model type"
    model1 = dl.dd_gauss
    model2 = dl.bg_hom3d

    model = lincombine(model1,model2)

    assert isinstance(model,Model)
# ======================================================================

# ======================================================================
def test_preserve_original(): 
    "Check that the original models are not changed by the function"
    model1 = dl.dd_gauss
    model2 = dl.bg_hom3d

    _ = lincombine(model1,model2)
    assert model1._parameter_list() == ['mean','std'] and model2._parameter_list() == ['conc','lam']
# ======================================================================

# ======================================================================
def test_twomodels_Nparam_nonlin(): 
    "Check that the combined model has the right number of parameters"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice

    model = lincombine(model1,model2)
    assert model.Nnonlin == model1.Nnonlin + model2.Nnonlin
# ======================================================================

# ======================================================================
def test_twomodels_Nparam_lin(): 
    "Check that the combined model has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = dl.dd_gauss3

    model = lincombine(model1,model2)
    assert model.Nlin == model1.Nlin + model2.Nlin
# ======================================================================

# ======================================================================
def test_twomodels_Nparam(): 
    "Check that the combined model has the right number of parameters"
    model1 = dl.dd_gauss3
    model2 = dl.dd_gauss2

    model = lincombine(model1,model2)
    assert model.Nparam == model1.Nparam + model2.Nparam
# ======================================================================

# ======================================================================
def test_twomodels_Nparam_list(): 
    "Check that the combined model has the right number of parameters"
    model1 = dl.dd_gauss3
    model2 = dl.dd_gauss2

    model = lincombine(model1,model2)
    assert model.Nparam == len(model._parameter_list())
# ======================================================================

# ======================================================================
def test_twomodels_param_names(): 
    "Check that the combined model has the adjusted parameter names"
    model1 = dl.dd_gauss
    model2 = dl.dd_gauss

    model = lincombine(model1,model2)
    assert all([ str in model._parameter_list() for str in ['mean_1','mean_2','std_1','std_2'] ])
# ======================================================================

# ======================================================================
def test_twomodels_default_linear(): 
    """Check that the default linear scaling parameter is added if there 
    are no linear parameters on one of the models"""
    model1 = dl.dd_gauss2
    model2 = dl.dd_gauss

    model = lincombine(model1,model2)
    assert 'scale_2' in model._parameter_list()
# ======================================================================

# ======================================================================
def test_twomodels_addweights(): 
    "Check that that weights can be introduced properly"
    model1 = dl.dd_gauss
    model2 = dl.dd_gauss
    model = lincombine(model1,model2,addweights=True)
    x = np.linspace(0,10,400)
    ref = model1(x,3,0.2) + model2(x,4,0.5)

    response = model(r_1=x,r_2=x,mean_1=3,std_1=0.2,
                         mean_2=4,std_2=0.5,
                         scale_1=1,scale_2=1,
                         weight_1=1,weight_2=1)

    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_twomodels_addweights_values(): 
    "Check that that weights values work properly"
    model1 = dl.dd_gauss
    model2 = dl.dd_gauss
    model = lincombine(model1,model2,addweights=True)
    x = np.linspace(0,10,400)
    ref1 = model1(x,3,0.2) 
    ref2 = model2(x,4,0.5)

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
def test_twomodels_call(): 
    "Check that that lincombine works correctly for two models"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    model = lincombine(model1,model2)
    x = np.linspace(0,10,400)
    ref1 = model1(x,3,0.2)
    ref2 = model2(x,4,0.5)

    ref = ref1 + ref2
    response = model(x,x,3,0.2,4,0.5,1,1)

    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_threemodels_Nparam_nonlin(): 
    "Check that the combined model has the right number of parameters"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    model3 = dl.dd_gauss2

    model = lincombine(model1,model2,model3)
    assert model.Nnonlin == model1.Nnonlin + model2.Nnonlin + model3.Nnonlin
# ======================================================================

# ======================================================================
def test_threemodels_Nparam_lin(): 
    "Check that the combined model has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = dl.dd_gauss3
    model3 = dl.dd_gauss2

    model = lincombine(model1,model2,model3)
    assert model.Nlin == model1.Nlin + model2.Nlin + model3.Nlin
# ======================================================================

# ======================================================================
def test_threemodels_Nparam(): 
    "Check that the combined model has the right number of parameters"
    model1 = dl.dd_gauss3
    model2 = dl.dd_gauss2
    model3 = dl.dd_gauss2

    model = lincombine(model1,model2,model3)
    assert model.Nparam == model1.Nparam + model2.Nparam + model3.Nparam
# ======================================================================

# ======================================================================
def test_threemodels_Nparam_list(): 
    "Check that the combined model has the right number of parameters"
    model1 = dl.dd_gauss3
    model2 = dl.dd_gauss2
    model3 = dl.dd_gauss2

    model = lincombine(model1,model2,model3)
    assert model.Nparam == len(model._parameter_list())
# ======================================================================

# ======================================================================
def test_threemodels_param_names(): 
    "Check that the combined model has the adjusted parameter names"
    model1 = dl.dd_gauss
    model2 = dl.dd_gauss
    model3 = dl.dd_gauss
    model = lincombine(model1,model2,model3)

    assert all([ str in model._parameter_list() for str in ['mean_1','mean_2','mean_3','std_1','std_2','std_3'] ])
# ======================================================================

# ======================================================================
def test_threemodels_default_linear(): 
    """Check that the default linear scaling parameter is added if there 
    are no linear parameters on one of the models"""
    model1 = dl.dd_gauss
    model2 = dl.dd_gauss
    model3 = dl.dd_gauss
    model = lincombine(model1,model2,model3)

    assert [scale_par in model._parameter_list() for scale_par in ['scale_1','scale_2','scale_3']]
# ======================================================================

# ======================================================================
def test_threemodels_addweights(): 
    "Check that that weights can be introduced properly"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    model3 = dl.dd_wormchain
    model = lincombine(model1,model2,model3,addweights=True)
    x = np.linspace(0,10,400)
    ref = model1(x,3,0.2) + model2(x,4,0.5) + model3(x,3.7,10)

    response = model(r_1=x,r_2=x,r_3=x,
                    mean_1=3,std_1=0.2,
                    location_2=4,spread_2=0.5,
                    contour_3=3.7,persistence_3=10,
                    scale_1=1,scale_2=1,scale_3=1,
                    weight_1=1,weight_2=1,weight_3=1)

    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_threemodels_addweights_values(): 
    "Check that that weights values work properly"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    model3 = dl.dd_gauss
    model = lincombine(model1,model2,model3,addweights=True)
    x = np.linspace(0,10,400)
    ref1 = model1(x,3,0.2) 
    ref2 = model2(x,4,0.5)
    ref3 = model3(x,5,0.1)

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
def test_threemodels_call(): 
    "Check that that lincombine works correctly for three models"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    model3 = dl.dd_shell
    model = lincombine(model1,model2,model3)
    x = np.linspace(0,10,400)
    ref1 = model1(x,3,0.2)
    ref2 = model2(x,4,0.5)
    ref3 = model3(x,2,0.1)

    ref = ref1+ref2+ref3
    response = model(x,x,x,3,0.2,4,0.5,2,0.1,1,1,1)

    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_fit_model(): 
    "Check that that lincombine works correctly for two models"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    model = lincombine(model1,model2)
    x = np.linspace(0,10,400)
    truth = model1(x,3,0.2)+model2(x,4,0.5)


    model.mean_1.par0=3
    model.location_2.par0=4
    result = fit(model,truth,x,x)

    assert np.allclose(result.model,truth)
# ======================================================================

model_vec = Model(lambda r: np.eye(len(r)),constants='r')
model_vec.addlinear('Pvec',vec=100,lb=0)

model_vec_normalized = Model(lambda r: np.eye(len(r)),constants='r')
model_vec_normalized.addlinear('Pvec',vec=100,lb=0,normalization=lambda Pvec: Pvec/np.sum(Pvec))

# ======================================================================
def test_vec_Nparam_nonlin(): 
    "Check that the combined model with a vector-form parameter has the right number of parameters"
    model1 = dl.dd_gauss
    model2 = model_vec

    model = lincombine(model1,model2)
    assert model.Nnonlin == model1.Nnonlin + model2.Nnonlin
# ======================================================================

# ======================================================================
def test_vec_Nparam_lin(): 
    "Check that the combined model with a vector-form parameter has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = model_vec

    model = lincombine(model1,model2)
    assert model.Nlin == model1.Nlin + model2.Nlin
# ======================================================================

# ======================================================================
def test_vec_Nparam(): 
    "Check that the combined model a vector-form parameter has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = model_vec

    model = lincombine(model1,model2)
    assert model.Nparam == model1.Nparam + model2.Nparam
# ======================================================================

# ======================================================================
def test_vec_param_names(): 
    "Check that the combined model has the adjusted parameter names"
    model1 = dl.dd_gauss
    model2 = model_vec

    model = lincombine(model1,model2)
    assert all([ str in model._parameter_list() for str in ['mean_1','std_1','Pvec_2'] ])
# ======================================================================

# ======================================================================
def test_vec_addweights(): 
    "Check that that weights can be introduced properly"
    model1 = model_vec
    model2 = model_vec
    model = lincombine(model1,model2,addweights=True)
    x = np.linspace(0,10,100)
    ref1 = dl.dd_gauss(x,3,0.2)
    ref2 = dl.dd_gauss(x,4,0.2)
    ref = ref1+ref2

    response = model(r_1=x,r_2=x,Pvec_1=ref1,Pvec_2=ref2,
                         weight_1=1,weight_2=1)

    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_vec_two_models(): 
    "Check that that lincombine works correctly for two models"
    model1 = dl.dd_gauss
    model2 = model_vec
    model = lincombine(model1,model2)
    x = np.linspace(0,10,100)
    ref1 = model1(x,3,0.2)
    ref2 = model2(r=x,Pvec=model1(x,4,0.3))
    ref = ref1 + ref2

    response = model(x,x,3,0.2,1,model1(x,4,0.3))

    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_vec_twomodels_normalization(): 
    """Check that the normalization of linear parameter is maintained"""
    model1 = model_vec_normalized
    model2 = model_vec_normalized
    model = lincombine(model1,model2)

    x = np.linspace(0,10,100)
    ref1 = model1(r=x,Pvec=dl.dd_gauss(x,4,0.3))
    ref2 = model2(r=x,Pvec=dl.dd_gauss(x,4,0.3))
    ref = ref1 + ref2

    results = fit(model,ref,x,x)

    assert hasattr(results,'Pvec_1_scale') and hasattr(results,'Pvec_2_scale')
# ======================================================================

# ======================================================================
def test_vec_twomodels_normalization_values(): 
    """Check that the normalization of linear parameter is maintained"""
    model1 = dl.dd_gauss
    model2 = model_vec_normalized
    model = lincombine(model1,model2)
    model.scale_1.freeze(1)
    model.mean_1.freeze(5)
    model.std_1.freeze(0.3)

    x = np.linspace(0,10,100)
    scale2 = 79
    dist2 = dl.dd_gauss(x,5,0.3)
    ref1 = model1(r=x,mean=2,std=0.3)
    ref2 = scale2*model2(r=x,Pvec=dist2/np.sum(dist2))
    ref = ref1 + ref2

    results = fit(model,ref,x,x)

    assert np.isclose(results.Pvec_2_scale,scale2)
# ======================================================================

def test_pickle_lincombined_model():
# ======================================================================
    "Check that the model object can be pickled"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    x = np.linspace(0,10,100)
    truth = model1(x,3,0.2)+model2(x,4,0.5)
    model = lincombine(model1,model2)
    model.mean_1.par0=3
    model.location_2.par0=4
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
