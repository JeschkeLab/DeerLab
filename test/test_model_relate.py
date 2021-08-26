import numpy as np
import deerlab as dl 
from deerlab.model import Model,fit,relate
from copy import deepcopy
# ======================================================================
def test_type(): 
    "Check that the function returns a valid model type"
    model = dl.dd_gauss

    newmodel = relate(model, width = lambda mean: mean/2)

    assert isinstance(newmodel,Model)
# ======================================================================

# ======================================================================
def test_preserve_original(): 
    "Check that the original input model is not changed by the function"
    model = dl.dd_gauss

    _ = relate(model, width = lambda mean: mean/2)
    assert model._parameter_list() == ['mean','width']
# ======================================================================

# ======================================================================
def test_Nparam_nonlin(): 
    "Check that the output model has the right number of parameters"
    model = dl.dd_gauss

    newmodel = relate(model, width = lambda mean: mean/2)
    assert newmodel.Nnonlin == model.Nnonlin - 1
# ======================================================================

# ======================================================================
def test_Nparam_lin(): 
    "Check that the combined model has the right number of parameters"
    model = dl.dd_gauss2

    newmodel = relate(model, width1 = lambda mean1: mean1/2)

    assert newmodel.Nlin == model.Nlin
# ======================================================================

# ======================================================================
def test_Nparam(): 
    "Check that the combined model has the right number of parameters"
    model = dl.dd_gauss2

    newmodel = relate(model, width1 = lambda mean1: mean1/2)
    assert newmodel.Nparam == model.Nparam - 1
# ======================================================================

# ======================================================================
def test_Nparam_list(): 
    "Check that the combined model has the right number of parameters"
    model = dl.dd_gauss3

    newmodel = relate(model, width1 = lambda mean1: mean1/2)
    assert newmodel.Nparam == len(newmodel._parameter_list())
# ======================================================================

# ======================================================================
def test_param_names(): 
    "Check that the combined model has the adjusted parameter names"
    model = dl.dd_gauss

    newmodel = relate(model, width = lambda mean: mean/2)
    assert all([ str in newmodel._parameter_list() for str in ['mean'] ])
# ======================================================================

# ======================================================================
def test_relate_one_nonlinear(): 
    "Check that that combine works correctly for one related parameter"
    model = dl.dd_gauss
    x = np.linspace(0,10,400)
    ref = model(x,4,0.2)

    newmodel = relate(model, width = lambda mean: mean/20)
    response = newmodel(x,4)

    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_relate_two_nonlinear(): 
    "Check that that combine works correctly for two related parameters"
    model = dl.dd_gauss2
    x = np.linspace(0,10,400)
    ref = model(x,3,0.2,6,0.1,1,1)

    newmodel = relate(model, 
                        mean1 = lambda mean2: mean2/2,
                        width1 = lambda width2: width2*2)
    response = newmodel(r=x,mean2=6,width2=0.1,amp1=1,amp2=1)

    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_relate_two_nonlinear_2(): 
    "Check that that combine works correctly for two related parameters"
    model = dl.dd_gauss2
    x = np.linspace(0,10,400)
    ref = model(x,3,0.2,3,0.1,1,1)

    newmodel = relate(model, 
                        mean1 = lambda width2, width1: 10*(width1+width2),
                        mean2 = lambda width2, width1: 10*(width1+width2),)
    response = newmodel(r=x,width1=0.2,width2=0.1,amp1=1,amp2=1)

    assert np.allclose(response,ref)
# ======================================================================


# ======================================================================
def test_relate_fit(): 
    "Check that that combine works correctly for one related parameter"
    model = dl.dd_gauss
    x = np.linspace(0,10,400)
    ref = model(x,4,0.2)

    newmodel = relate(model, width = lambda mean: mean/20)

    result = fit(newmodel,ref,x)

    assert np.allclose(result.model,ref)
# ======================================================================

# ======================================================================
def test_relate_addednonlinear(): 
    "Check that a parameter can be related with an added nonlinear parameter"
    model = deepcopy(dl.dd_gauss)
    model.addlinear('amp')
    x = np.linspace(0,10,400)
    ref = model(x,4,0.2,5)
    
    model.addnonlinear('factor')
    newmodel = relate(model, width = lambda factor: factor/5)
    response = newmodel(r=x,mean=4,factor=1,amp=5)

    assert np.allclose(response,ref)
# ======================================================================

