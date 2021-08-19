import numpy as np
import deerlab as dl 
from deerlab.model import Model,fit,combine

# ======================================================================
def test_type(): 
    "Check that the function returns a valid model type"
    model1 = dl.dd_gauss
    model2 = dl.bg_hom3d

    model = combine(model1,model2)

    assert isinstance(model,Model)
# ======================================================================

# ======================================================================
def test_preserve_original(): 
    "Check that the original models are not changed by the function"
    model1 = dl.dd_gauss
    model2 = dl.bg_hom3d

    _ = combine(model1,model2)
    assert model1._parameter_list() == ['mean','width'] and model2._parameter_list() == ['conc','lam']
# ======================================================================

# ======================================================================
def test_Nparam_nonlin(): 
    "Check that the combined model has the right number of parameters"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice

    model = combine(model1,model2)
    assert model.Nnonlin == model1.Nnonlin + model2.Nnonlin
# ======================================================================

# ======================================================================
def test_Nparam_lin(): 
    "Check that the combined model has the right number of parameters"
    model1 = dl.dd_gauss2
    model2 = dl.dd_gauss3

    model = combine(model1,model2)
    assert model.Nlin == model1.Nlin + model2.Nlin
# ======================================================================

# ======================================================================
def test_Nparam(): 
    "Check that the combined model has the right number of parameters"
    model1 = dl.dd_gauss3
    model2 = dl.dd_gauss2

    model = combine(model1,model2)
    assert model.Nparam == model1.Nparam + model2.Nparam
# ======================================================================

# ======================================================================
def test_Nparam_list(): 
    "Check that the combined model has the right number of parameters"
    model1 = dl.dd_gauss3
    model2 = dl.dd_gauss2

    model = combine(model1,model2)
    assert model.Nparam == len(model._parameter_list())
# ======================================================================

# ======================================================================
def test_param_names(): 
    "Check that the combined model has the adjusted parameter names"
    model1 = dl.dd_gauss
    model2 = dl.dd_gauss

    model = combine(model1,model2)
    assert all([ str in model._parameter_list() for str in ['mean_1','mean_2','width_1','width_2'] ])
# ======================================================================

# ======================================================================
def test_default_linear(): 
    """Check that the default linear scaling parameter is added if there 
    are no linear parameters on one of the models"""
    model1 = dl.dd_gauss2
    model2 = dl.dd_gauss

    model = combine(model1,model2)
    assert 'scale_2' in model._parameter_list()
# ======================================================================

# ======================================================================
def test_two_models(): 
    "Check that that combine works correctly for two models"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    model = combine(model1,model2)
    x = np.linspace(0,10,400)
    ref1 = model1(x,3,0.2)
    ref2 = model2(x,4,0.5)

    response = model(x,x,3,0.2,4,0.5,1,1)

    assert all([np.allclose(response[n],ref) for n,ref in enumerate([ref1,ref2])])
# ======================================================================

# ======================================================================
def test_three_models(): 
    "Check that that combine works correctly for three models"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    model3 = dl.dd_shell
    model = combine(model1,model2,model3)
    x = np.linspace(0,10,400)
    ref1 = model1(x,3,0.2)
    ref2 = model2(x,4,0.5)
    ref3 = model3(x,2,0.1)

    response = model(x,x,x,3,0.2,4,0.5,2,0.1,1,1,1)

    assert all([np.allclose(response[n],ref) for n,ref in enumerate([ref1,ref2,ref3])])
# ======================================================================

# ======================================================================
def test_fit_model(): 
    "Check that that combine works correctly for two models"
    model1 = dl.dd_gauss
    model2 = dl.dd_rice
    model = combine(model1,model2)
    x = np.linspace(0,10,400)
    truth = [model1(x,3,0.2),model2(x,4,0.5)]
    result = fit(model,truth,x,x)

    assert np.allclose(result.model[0],truth[0]) and np.allclose(result.model[1],truth[1])
# ======================================================================

# ======================================================================
def test_link_same_name(): 
    "Check that the parameters with the same name are linked properly"
    model1 = dl.dd_gauss
    model2 = dl.dd_gauss

    model = combine(model1,model2,mean=[model1.mean,model2.mean])

    assert hasattr(model,'mean') and not hasattr(model,'mean_1') and not hasattr(model,'mean_2')
# ======================================================================

# ======================================================================
def test_link_name(): 
    "Check that the parameters are linked to the proper name"
    model1 = dl.dd_gauss
    model2 = dl.dd_gauss2

    model = combine(model1,model2,mean=[model1.mean,model2.mean1])

    assert hasattr(model,'mean') and not hasattr(model,'mean_1') and not hasattr(model,'mean1_2')
# ======================================================================

# ======================================================================
def test_link_Nparam(): 
    "Check that the parameters are linked resulting in proper number of parameters"
    model1 = dl.dd_gauss3
    model2 = dl.dd_gauss2

    model = combine(model1,model2,mean1=[model1.mean2,model2.mean1])

    assert model.Nparam == model1.Nparam + model2.Nparam - 1
# ======================================================================

# ======================================================================
def test_link_Nparam_list(): 
    "Check that the combined model has the right number of parameters"
    model1 = dl.dd_gauss3
    model2 = dl.dd_gauss2

    model = combine(model1,model2,mean1=[model1.mean2,model2.mean1])
    assert model.Nparam == len(model._parameter_list())
# ======================================================================

# ======================================================================
def test_link_multiple(): 
    "Check that multiple links can be established"
    model1 = dl.dd_gauss3
    model2 = dl.dd_gauss2

    model = combine(model1,model2,
            mean1=[model1.mean1,model2.mean1],
            mean2=[model1.mean2,model2.mean2],
            width2=[model1.width2,model2.width2])

    assert model.Nparam == model1.Nparam + model2.Nparam - 3
# ======================================================================

# ======================================================================
def test_link_call(): 
    "Check that linked parameter models return the correct responses"
    model1 = dl.dd_gauss
    model2 = dl.dd_gauss2

    model = combine(model1,model2,mean=[model1.mean,model2.mean2])

    x = np.linspace(0,10,400)
    ref1 = model1(x,3,0.2)
    ref2 = model2(x,4,0.5,3,0.4,1,1)

    response = model(x,x,3,0.2,4,0.5,0.4,1,1,1)

    assert all([np.allclose(response[n],ref) for n,ref in enumerate([ref1,ref2])])
# ======================================================================
test_link_call()

# ======================================================================
def test_link_fit(): 
    "Check that linked parameter models can be properly fitted"
    model1 = dl.dd_gauss
    model2 = dl.dd_gauss
    model3 = dl.dd_gauss2

    model = combine(model1,model2,model3,
                    mean1=[model1.mean,model3.mean1],
                    mean2=[model2.mean,model3.mean2])
    model.mean1.par0 = 3
    model.mean2.par0 = 5

    x = np.linspace(0,10,400)
    ref1 = model1(x,3,0.2)
    ref2 = model2(x,5,0.2)
    ref3 = model3(x,5,0.5,3,0.4,1,1)

    result = fit(model,[ref1,ref2,ref3],x,x,x)
    
    assert all([np.allclose(result.model[n],ref) for n,ref in enumerate([ref1,ref2,ref3])])
# ======================================================================

# ======================================================================
def test_link_single(): 
    "Check that the parameters can be linked within a single model"
    model1 = dl.dd_gauss2

    model = combine(model1,means=[model1.mean1,model1.mean2])

    assert model.Nparam == model1.Nparam - 1
# ======================================================================
