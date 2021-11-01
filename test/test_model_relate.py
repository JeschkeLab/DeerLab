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
def test_relate_signature(): 
    "Check the model signature after relate"
    model = dl.dd_gauss2

    newmodel = relate(model,mean1= lambda mean2: mean2)

    assert newmodel.signature==['r', 'mean1', 'width1', 'width2', 'amp1', 'amp2']
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
    "Check that the related model has the right number of parameters"
    model = dl.dd_gauss2

    newmodel = relate(model, width1 = lambda mean1: mean1/2)

    assert newmodel.Nlin == model.Nlin
# ======================================================================

# ======================================================================
def test_Nparam(): 
    "Check that the related model has the right number of parameters"
    model = dl.dd_gauss2

    newmodel = relate(model, width1 = lambda mean1: mean1/2)
    assert newmodel.Nparam == model.Nparam - 1
# ======================================================================

# ======================================================================
def test_Nparam_list(): 
    "Check that the related model has the right number of parameters"
    model = dl.dd_gauss3

    newmodel = relate(model, width1 = lambda mean1: mean1/2)
    assert newmodel.Nparam == len(newmodel._parameter_list())
# ======================================================================

# ======================================================================
def test_param_names(): 
    "Check that the related model has the adjusted parameter names"
    model = dl.dd_gauss

    newmodel = relate(model, width = lambda mean: mean/2)
    assert all([ str in newmodel._parameter_list() for str in ['mean'] ])
# ======================================================================

# ======================================================================
def test_relate_one_nonlinear(): 
    "Check that the method works correctly for one related parameter"
    model = dl.dd_gauss
    x = np.linspace(0,10,400)
    ref = model(x,4,0.2)

    newmodel = relate(model, width = lambda mean: mean/20)
    response = newmodel(x,4)

    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_relate_two_nonlinear(): 
    "Check that the method works correctly for two related parameters"
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
    "Check that the method works correctly for two related parameters"
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
    "Check that the method works correctly for one related parameter"
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

# ======================================================================
def test_relate_conflict_order(): 
    "Check that even if the user specifies a conflicted order, it runs"
    model = dl.dd_gauss2
    x = np.linspace(0,10,400)
    ref = model(x,4,0.2,5,0.4,1,1)

    # No conflict, there parameters are ordered correctly
    newmodel = relate(model, 
                    width2=lambda width1: width1*2,
                    width1=lambda mean1: mean1/20,)

    # Conflict: width1 is deleted first, but width2 depends on it
    conflictmodel = relate(model, 
                    width1=lambda mean1: mean1/20,
                    width2=lambda width1: width1*2)
    
    response1 = newmodel(r=x,mean1=4,mean2=5,amp1=1,amp2=1)
    response2 = conflictmodel(r=x,mean1=4,mean2=5,amp1=1,amp2=1)

    assert np.allclose(response1,ref) and np.allclose(response2,ref)
# ======================================================================

