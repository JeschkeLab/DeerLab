import numpy as np
import deerlab as dl 
from deerlab.model import Model,fit,expand

# ======================================================================
def test_type(): 
    "Check that the function returns a valid model type"
    model1 = dl.dd_gauss
    model2 = dl.bg_hom3d

    model = expand(model1,model2)

    assert isinstance(model,Model)
# ======================================================================

# ======================================================================
def test_preserve_original(): 
    "Check that the original models are not changed by the function"
    model1 = dl.dd_gauss
    model2 = dl.bg_hom3d

    _ = expand(model1,model2)
    assert model1._parameter_list() == ['mean','width'] and model2._parameter_list() == ['conc','lam']
# ======================================================================

# ======================================================================
def test_twomodels_Nparam_nonlin(): 
    "Check that the expandd model has the right number of parameters"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice

    model = expand(model1,model2)
    assert model.Nnonlin == model1.Nnonlin + model2.Nnonlin
# ======================================================================

# ======================================================================
def test_twomodels_Nparam_lin(): 
    "Check that the expandd model has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = dl.dd_gauss3

    model = expand(model1,model2)
    assert model.Nlin == model1.Nlin + model2.Nlin
# ======================================================================

# ======================================================================
def test_twomodels_Nparam(): 
    "Check that the expandd model has the right number of parameters"
    model1 = dl.dd_gauss3
    model2 = dl.dd_gauss2

    model = expand(model1,model2)
    assert model.Nparam == model1.Nparam + model2.Nparam
# ======================================================================

# ======================================================================
def test_twomodels_Nparam_list(): 
    "Check that the expandd model has the right number of parameters"
    model1 = dl.dd_gauss3
    model2 = dl.dd_gauss2

    model = expand(model1,model2)
    assert model.Nparam == len(model._parameter_list())
# ======================================================================

# ======================================================================
def test_twomodels_param_names(): 
    "Check that the expandd model has the adjusted parameter names"
    model1 = dl.dd_gauss
    model2 = dl.dd_gauss

    model = expand(model1,model2)
    assert all([ str in model._parameter_list() for str in ['mean_1','mean_2','width_1','width_2'] ])
# ======================================================================

# ======================================================================
def test_twomodels_default_linear(): 
    """Check that the default linear scaling parameter is added if there 
    are no linear parameters on one of the models"""
    model1 = dl.dd_gauss2
    model2 = dl.dd_gauss

    model = expand(model1,model2)
    assert 'scale_2' in model._parameter_list()
# ======================================================================

# ======================================================================
def test_twomodels_addweights(): 
    "Check that that weights can be introduced properly"
    model1 = dl.dd_gauss
    model2 = dl.dd_gauss
    model = expand(model1,model2,addweights=True)
    x = np.linspace(0,10,400)
    ref1 = model1(x,3,0.2)
    ref2 = model2(x,4,0.5)

    response = model(r_1=x,r_2=x,mean_1=3,width_1=0.2,
                         mean_2=4,width_2=0.5,
                         scale_1=1,scale_2=1,
                         weight_1=1,weight_2=1)

    assert all([np.allclose(response[n],ref) for n,ref in enumerate([ref1,ref2])])
# ======================================================================

# ======================================================================
def test_twomodels_call(): 
    "Check that that expand works correctly for two models"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    model = expand(model1,model2)
    x = np.linspace(0,10,400)
    ref1 = model1(x,3,0.2)
    ref2 = model2(x,4,0.5)

    response = model(x,x,3,0.2,4,0.5,1,1)

    assert all([np.allclose(response[n],ref) for n,ref in enumerate([ref1,ref2])])
# ======================================================================

# ======================================================================
def test_threemodels_Nparam_nonlin(): 
    "Check that the expandd model has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = dl.dd_gauss2
    model3 = dl.dd_gauss3

    model = expand(model1,model2,model3)

    assert model.Nnonlin == model1.Nnonlin + model2.Nnonlin + model3.Nnonlin
# ======================================================================

# ======================================================================
def test_threemodels_Nparam_lin(): 
    "Check that the expandd model has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = dl.dd_gauss2
    model3 = dl.dd_gauss3

    model = expand(model1,model2,model3)
    
    assert model.Nlin == model1.Nlin + model2.Nlin + model3.Nlin
# ======================================================================

# ======================================================================
def test_threemodels_Nparam(): 
    "Check that the expandd model has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = dl.dd_gauss2
    model3 = dl.dd_gauss3

    model = expand(model1,model2,model3)

    assert model.Nparam == model1.Nparam + model2.Nparam + model3.Nparam
# ======================================================================

# ======================================================================
def test_threemodels_Nparam_list(): 
    "Check that the expandd model has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = dl.dd_gauss2
    model3 = dl.dd_gauss3

    model = expand(model1,model2,model3)

    assert model.Nparam == len(model._parameter_list())
# ======================================================================

# ======================================================================
def test_threemodels_param_names(): 
    "Check that the expandd model has the adjusted parameter names"
    model1 = dl.dd_gauss
    model2 = dl.dd_gauss
    model3 = dl.dd_gauss

    model = expand(model1,model2,model3)

    assert all([ str in model._parameter_list() for str in ['mean_1','mean_2','mean_3','width_1','width_2','width_3'] ])
# ======================================================================

# ======================================================================
def test_threemodels_default_linear(): 
    """Check that the default linear scaling parameter is added if there 
    are no linear parameters on one of the models"""
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    model3 = dl.dd_shell

    model = expand(model1,model2,model3)

    assert [scaleparam in model._parameter_list() for scaleparam in ['scale_1','scale_2','scale_3']]
# ======================================================================

# ======================================================================
def test_threemodels_addweights(): 
    "Check that that weights can be introduced properly"
    model1 = dl.dd_gauss
    model2 = dl.dd_gauss
    model = expand(model1,model2,addweights=True)
    x = np.linspace(0,10,400)
    ref1 = model1(x,3,0.2)
    ref2 = model2(x,4,0.5)

    response = model(r_1=x,r_2=x,mean_1=3,width_1=0.2,
                         mean_2=4,width_2=0.5,
                         scale_1=1,scale_2=1,
                         weight_1=1,weight_2=1)

    assert all([np.allclose(response[n],ref) for n,ref in enumerate([ref1,ref2])])
# ======================================================================

# ======================================================================
def test_threemodels_call(): 
    "Check that that expand works correctly for three models"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    model3 = dl.dd_shell
    model = expand(model1,model2,model3)
    x = np.linspace(0,10,400)
    ref1 = model1(x,3,0.2)
    ref2 = model2(x,4,0.5)
    ref3 = model3(x,2,0.1)

    response = model(x,x,x,3,0.2,4,0.5,2,0.1,1,1,1)

    assert all([np.allclose(response[n],ref) for n,ref in enumerate([ref1,ref2,ref3])])
# ======================================================================

# ======================================================================
def test_fit_model(): 
    "Check that that expand works correctly for two models"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    model = expand(model1,model2)
    x = np.linspace(0,10,400)
    truth = [model1(x,3,0.2),model2(x,4,0.5)]
    result = fit(model,truth,x,x)

    assert np.allclose(result.model[0],truth[0]) and np.allclose(result.model[1],truth[1])
# ======================================================================

model_vec = Model(lambda r: np.eye(len(r)),constants='r')
model_vec.addlinear('Pvec',vec=100,lb=0)

# ======================================================================
def test_vec_Nparam_nonlin(): 
    "Check that the expandd model with a vector-form parameter has the right number of parameters"
    model1 = dl.dd_gauss
    model2 = model_vec

    model = expand(model1,model2)
    assert model.Nnonlin == model1.Nnonlin + model2.Nnonlin
# ======================================================================

# ======================================================================
def test_vec_Nparam_lin(): 
    "Check that the expandd model with a vector-form parameter has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = model_vec

    model = expand(model1,model2)
    assert model.Nlin == model1.Nlin + model2.Nlin
# ======================================================================

# ======================================================================
def test_vec_Nparam(): 
    "Check that the expandd model a vector-form parameter has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = model_vec

    model = expand(model1,model2)
    assert model.Nparam == model1.Nparam + model2.Nparam
# ======================================================================

# ======================================================================
def test_vec_param_names(): 
    "Check that the expandd model has the adjusted parameter names"
    model1 = dl.dd_gauss
    model2 = model_vec

    model = expand(model1,model2)
    assert all([ str in model._parameter_list() for str in ['mean_1','width_1','Pvec_2'] ])
# ======================================================================

# ======================================================================
def test_vec_twomodels(): 
    "Check that that expand works correctly for two vector-based models"
    model1 = model_vec
    model2 = model_vec
    model = expand(model1,model2)
    x = np.linspace(0,10,100)
    ref1 = model1(r=x,Pvec=dl.dd_gauss(x,5,0.2))
    ref2 = model2(r=x,Pvec=dl.dd_gauss(x,4,0.3))

    response = model(x,x,dl.dd_gauss(x,5,0.2),dl.dd_gauss(x,4,0.3))

    assert all([np.allclose(response[n],ref) for n,ref in enumerate([ref1,ref2])])
# ======================================================================

# ======================================================================
def test_vec_twomodels_mixed(): 
    "Check that that expand works correctly for two vector-based models"
    model1 = dl.dd_gauss
    model2 = model_vec
    model = expand(model1,model2)
    x = np.linspace(0,10,100)
    ref1 = model1(x,3,0.2)
    ref2 = model2(r=x,Pvec=dl.dd_gauss(x,4,0.3))

    response = model(x,x,3,0.2,1,dl.dd_gauss(x,4,0.3))

    assert all([np.allclose(response[n],ref) for n,ref in enumerate([ref1,ref2])])
# ======================================================================

