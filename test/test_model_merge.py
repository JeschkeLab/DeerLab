import numpy as np
import deerlab as dl 
from deerlab.model import Model,fit,merge
from deerlab.utils import store_pickle, read_pickle
import os 

# ======================================================================
def test_type(): 
    "Check that the function returns a valid model type"
    model1 = dl.dd_gauss
    model2 = dl.bg_hom3d

    model = merge(model1,model2)

    assert isinstance(model,Model)
# ======================================================================

# ======================================================================
def test_preserve_original(): 
    "Check that the original models are not changed by the function"
    model1 = dl.dd_gauss
    model2 = dl.bg_hom3d

    _ = merge(model1,model2)
    assert model1._parameter_list() == ['mean','std'] and model2._parameter_list() == ['conc','lam']
# ======================================================================

# ======================================================================
def test_twomodels_Nparam_nonlin(): 
    "Check that the merged model has the right number of parameters"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice

    model = merge(model1,model2)
    assert model.Nnonlin == model1.Nnonlin + model2.Nnonlin
# ======================================================================

# ======================================================================
def test_twomodels_Nparam_lin(): 
    "Check that the merged model has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = dl.dd_gauss3

    model = merge(model1,model2)
    assert model.Nlin == model1.Nlin + model2.Nlin
# ======================================================================

# ======================================================================
def test_twomodels_Nparam(): 
    "Check that the merged model has the right number of parameters"
    model1 = dl.dd_gauss3
    model2 = dl.dd_gauss2

    model = merge(model1,model2)
    assert model.Nparam == model1.Nparam + model2.Nparam
# ======================================================================

# ======================================================================
def test_twomodels_Nparam_list(): 
    "Check that the merged model has the right number of parameters"
    model1 = dl.dd_gauss3
    model2 = dl.dd_gauss2

    model = merge(model1,model2)
    assert model.Nparam == len(model._parameter_list())
# ======================================================================

# ======================================================================
def test_twomodels_param_names(): 
    "Check that the merged model has the adjusted parameter names"
    model1 = dl.dd_gauss
    model2 = dl.dd_gauss

    model = merge(model1,model2)
    assert all([ str in model._parameter_list() for str in ['mean_1','mean_2','std_1','std_2'] ])
# ======================================================================

# ======================================================================
def test_twomodels_default_linear(): 
    """Check that the default linear scaling parameter is added if there 
    are no linear parameters on one of the models"""
    model1 = dl.dd_gauss2
    model2 = dl.dd_gauss

    model = merge(model1,model2)
    assert 'scale_2' in model._parameter_list()
# ======================================================================

# ======================================================================
def test_twomodels_addweights(): 
    "Check that that weights can be introduced properly"
    model1 = dl.dd_gauss
    model2 = dl.dd_gauss
    model = merge(model1,model2,addweights=True)
    x = np.linspace(0,10,400)
    ref1 = model1(x,3,0.2)
    ref2 = model2(x,4,0.5)

    response = model(r_1=x,r_2=x,mean_1=3,std_1=0.2,
                         mean_2=4,std_2=0.5,
                         scale_1=1,scale_2=1,
                         weight_1=1,weight_2=1)

    assert all([np.allclose(response[n],ref) for n,ref in enumerate([ref1,ref2])])
# ======================================================================

# ======================================================================
def test_twomodels_call(): 
    "Check that that merge works correctly for two models"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    model = merge(model1,model2)
    x = np.linspace(0,10,400)
    ref1 = model1(x,3,0.2)
    ref2 = model2(x,4,0.5)

    response = model(x,x,3,0.2,4,0.5,1,1)

    assert all([np.allclose(response[n],ref) for n,ref in enumerate([ref1,ref2])])
# ======================================================================

# ======================================================================
def test_threemodels_Nparam_nonlin(): 
    "Check that the merged model has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = dl.dd_gauss2
    model3 = dl.dd_gauss3

    model = merge(model1,model2,model3)

    assert model.Nnonlin == model1.Nnonlin + model2.Nnonlin + model3.Nnonlin
# ======================================================================

# ======================================================================
def test_threemodels_Nparam_lin(): 
    "Check that the merged model has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = dl.dd_gauss2
    model3 = dl.dd_gauss3

    model = merge(model1,model2,model3)
    
    assert model.Nlin == model1.Nlin + model2.Nlin + model3.Nlin
# ======================================================================

# ======================================================================
def test_threemodels_Nparam(): 
    "Check that the merged model has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = dl.dd_gauss2
    model3 = dl.dd_gauss3

    model = merge(model1,model2,model3)

    assert model.Nparam == model1.Nparam + model2.Nparam + model3.Nparam
# ======================================================================

# ======================================================================
def test_threemodels_Nparam_list(): 
    "Check that the merged model has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = dl.dd_gauss2
    model3 = dl.dd_gauss3

    model = merge(model1,model2,model3)

    assert model.Nparam == len(model._parameter_list())
# ======================================================================

# ======================================================================
def test_threemodels_param_names(): 
    "Check that the merged model has the adjusted parameter names"
    model1 = dl.dd_gauss
    model2 = dl.dd_gauss
    model3 = dl.dd_gauss

    model = merge(model1,model2,model3)

    assert all([ str in model._parameter_list() for str in ['mean_1','mean_2','mean_3','std_1','std_2','std_3'] ])
# ======================================================================

# ======================================================================
def test_threemodels_default_linear(): 
    """Check that the default linear scaling parameter is added if there 
    are no linear parameters on one of the models"""
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    model3 = dl.dd_shell

    model = merge(model1,model2,model3)

    assert [scaleparam in model._parameter_list() for scaleparam in ['scale_1','scale_2','scale_3']]
# ======================================================================

# ======================================================================
def test_threemodels_addweights(): 
    "Check that that weights can be introduced properly"
    model1 = dl.dd_gauss
    model2 = dl.dd_gauss
    model = merge(model1,model2,addweights=True)
    x = np.linspace(0,10,400)
    ref1 = model1(x,3,0.2)
    ref2 = model2(x,4,0.5)

    response = model(r_1=x,r_2=x,mean_1=3,std_1=0.2,
                         mean_2=4,std_2=0.5,
                         scale_1=1,scale_2=1,
                         weight_1=1,weight_2=1)

    assert all([np.allclose(response[n],ref) for n,ref in enumerate([ref1,ref2])])
# ======================================================================

# ======================================================================
def test_threemodels_call(): 
    "Check that that merge works correctly for three models"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    model3 = dl.dd_shell
    model = merge(model1,model2,model3)
    x = np.linspace(0,10,400)
    ref1 = model1(x,3,0.2)
    ref2 = model2(x,4,0.5)
    ref3 = model3(x,2,0.1)

    response = model(x,x,x,3,0.2,4,0.5,2,0.1,1,1,1)

    assert all([np.allclose(response[n],ref) for n,ref in enumerate([ref1,ref2,ref3])])
# ======================================================================

# ======================================================================
def test_fit_model(): 
    "Check that that merge works correctly for two models"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    model = merge(model1,model2)
    x = np.linspace(0,10,400)
    truth = [model1(x,3,0.2),
             model2(x,4,0.5)]
    result = fit(model,truth,x,x,weights=[1,1])

    assert np.allclose(result.model[0],truth[0]) and np.allclose(result.model[1],truth[1])
# ======================================================================

model_vec = Model(lambda r: np.eye(len(r)),constants='r')
model_vec.addlinear('Pvec',vec=100,lb=0)

model_vec_normalized = Model(lambda r: np.eye(len(r)),constants='r')
model_vec_normalized.addlinear('Pvec',vec=100,lb=0,normalization = lambda Pvec: Pvec/np.sum(Pvec))

# ======================================================================
def test_vec_Nparam_nonlin(): 
    "Check that the merged model with a vector-form parameter has the right number of parameters"
    model1 = dl.dd_gauss
    model2 = model_vec

    model = merge(model1,model2)
    assert model.Nnonlin == model1.Nnonlin + model2.Nnonlin
# ======================================================================

# ======================================================================
def test_vec_Nparam_lin(): 
    "Check that the merged model with a vector-form parameter has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = model_vec

    model = merge(model1,model2)
    assert model.Nlin == model1.Nlin + model2.Nlin
# ======================================================================

# ======================================================================
def test_vec_Nparam(): 
    "Check that the merged model a vector-form parameter has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = model_vec

    model = merge(model1,model2)
    assert model.Nparam == model1.Nparam + model2.Nparam
# ======================================================================

# ======================================================================
def test_vec_param_names(): 
    "Check that the merged model has the adjusted parameter names"
    model1 = dl.dd_gauss
    model2 = model_vec

    model = merge(model1,model2)
    assert all([ str in model._parameter_list() for str in ['mean_1','std_1','Pvec_2'] ])
# ======================================================================

# ======================================================================
def test_vec_twomodels(): 
    "Check that that merge works correctly for two vector-based models"
    model1 = model_vec
    model2 = model_vec
    model = merge(model1,model2)
    x = np.linspace(0,10,100)
    ref1 = model1(r=x,Pvec=dl.dd_gauss(x,5,0.2))
    ref2 = model2(r=x,Pvec=dl.dd_gauss(x,4,0.3))

    response = model(x,x,dl.dd_gauss(x,5,0.2),dl.dd_gauss(x,4,0.3))

    assert all([np.allclose(response[n],ref) for n,ref in enumerate([ref1,ref2])])
# ======================================================================

# ======================================================================
def test_vec_twomodels_mixed(): 
    "Check that that merge works correctly for two vector-based models"
    model1 = dl.dd_gauss
    model2 = model_vec
    model = merge(model1,model2)
    x = np.linspace(0,10,100)
    ref1 = model1(x,3,0.2)
    ref2 = model2(r=x,Pvec=dl.dd_gauss(x,4,0.3))

    response = model(x,x,3,0.2,1,dl.dd_gauss(x,4,0.3))

    assert all([np.allclose(response[n],ref) for n,ref in enumerate([ref1,ref2])])
# ======================================================================

# ======================================================================
def test_vec_twomodels_mixed_normalized(): 
    "Check that that merge works correctly for mixed models"
    model1 = dl.dd_gauss
    model2 = model_vec_normalized
    model = merge(model1,model2)
    x = np.linspace(0,10,100)
    ref1 = model1(x,3,0.2)
    ref2 = model2(r=x,Pvec=dl.dd_gauss(x,4,0.3))

    response = model(x,x,3,0.2,1,dl.dd_gauss(x,4,0.3))

    assert all([np.allclose(response[n],ref) for n,ref in enumerate([ref1,ref2])])
# ======================================================================

# ======================================================================
def test_twomodels_normalization(): 
    """Check that the normalization of linear parameter is maintained"""
    model1 = model_vec_normalized
    model2 = model_vec_normalized
    model = merge(model1,model2)

    x = np.linspace(0,10,100)
    ref1 = model1(r=x,Pvec=dl.dd_gauss(x,4,0.3))
    ref2 = model2(r=x,Pvec=dl.dd_gauss(x,4,0.3))

    results = fit(model,[ref1,ref2],x,x)

    assert hasattr(results,'Pvec_1_scale') and hasattr(results,'Pvec_2_scale')
# ======================================================================

# ======================================================================
def test_twomodels_normalization_values(): 
    """Check that the normalization of linear parameter is maintained"""
    model1 = model_vec_normalized
    model2 = model_vec_normalized
    model = merge(model1,model2)

    x = np.linspace(0,10,100)
    scale1,scale2 = 55,79
    dist = dl.dd_gauss(x,4,0.3)
    ref1 = scale1*model1(r=x,Pvec=dist/np.sum(dist))
    ref2 = scale2*model2(r=x,Pvec=dist/np.sum(dist))

    results = fit(model,[ref1,ref2],x,x)

    assert np.isclose(results.Pvec_1_scale,scale1) and np.isclose(results.Pvec_2_scale,scale2)
# ======================================================================

# ======================================================================
def test_merge_linked():
    "Check that that merge works correctly for models with linked parameters"
    submodel1 = dl.dd_gauss
    submodel2 = dl.dd_gauss2
    submodel2 = dl.link(submodel2, mean=['mean1','mean2'], std=['std1','std2'])

    x = np.linspace(0,10,100)
    ref1 = submodel1(x,3,0.2)
    ref2 = submodel1(x,4,0.3)

    # Create a global model
    model = dl.merge(submodel1,submodel2, addweights=True)

    response = model(x,x,mean_1=3,std_1=0.2,mean_2=4,std_2=0.3,amp1_2=1,amp2_2=1,scale_1=1, weight_1=1, weight_2=1)
    assert all([np.allclose(response[n],ref) for n,ref in enumerate([ref1,ref2])])
# ======================================================================

# ======================================================================
def test_twomodels_vec_equal_global(): 
    """Check that fitting two equal models results in the same fit as the local model"""
    model = model_vec
    globalmodel = merge(model,model)

    x = np.linspace(0,10,100)
    ref = model(r=x,Pvec=dl.dd_gauss(x,4,0.3))

    results = fit(model,ref,x)
    globalresults = fit(globalmodel,[ref,ref],x,x)

    assert np.allclose(results.model,globalresults.model) and np.allclose(results.Pvec,globalresults.Pvec_1,globalresults.Pvec_2)
# ======================================================================

# ======================================================================
def test_twomodels_equal_global(): 
    "Check that that merge works correctly for two models"
    model = dl.dd_gauss
    globalmodel = merge(model,model)
    x = np.linspace(0,10,400)
    truth = model(x,3,0.2)

    result = fit(model,truth,x)
    globalresult = fit(globalmodel,[truth,truth],x,x)

    assert np.allclose(result.model,truth,globalresult.model) and np.allclose(result.mean,globalresult.mean_1)
# ======================================================================

def test_pickle_merged_model():
# ======================================================================
    "Check that the model object can be pickled"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    x = np.linspace(0,10,100)
    truth = [model1(x,3,0.2),
            model2(x,4,0.5)]
    model = merge(model1,model2)
    store_pickle(model,'pickled_model')
    try:
        pickled_model =read_pickle('pickled_model')
        result = fit(pickled_model,truth,x,x,weights=[1,1])

        os.remove("pickled_model.pkl") 
        assert np.allclose(result.model[0],truth[0]) and np.allclose(result.model[1],truth[1])
    except Exception as exception:
        os.remove("pickled_model.pkl") 
        raise exception
# ======================================================================
