from collections import namedtuple
from deerlab.whitegaussnoise import whitegaussnoise
from deerlab.model import Model, fit
import numpy as np 

# Simple non-linear function for testing
x = np.linspace(0,5,100)
def gauss(mean,width): 
    return np.exp(-(x-mean)**2/width**2/2)

# Non-linear definition
def gauss2(mean1,mean2,width1,width2,amp1,amp2):
    return amp1*gauss(mean1,width1) + amp2*gauss(mean2,width2)

# Linear + Non-linear definition
def gauss2_design(mean1,mean2,width1,width2):
    return np.atleast_2d([gauss(mean1,width1), gauss(mean2,width2)]).T

# Linear + Non-linear definition
def gauss2_identity():
    return np.eye(len(x))

# Linear + Non-linear definition
def gauss2_scaled(scale):
    return scale*np.eye(len(x))

mock_data = gauss2(mean1=3,mean2=4,width1=0.5,width2=0.2,amp1=0.5,amp2=0.6)


def test_construction_length():
#================================================================
    "Check that the model is contructed correctly with the appropiate number of parameters"
    model = Model(gauss)

    assert model.Nparam==2
#================================================================

def test_construction_names():
#================================================================
    "Check that the model is contructed correctly with the correct parameter names"
    model = Model(gauss)

    assert 'mean' in model.__dict__ and 'width' in model.__dict__
#================================================================

def test_parameters_set():
#================================================================
    "Check that attributes of the parameters are editable"
    model = Model(gauss)
    model.mean.set(lb=0,ub=10)

    assert getattr(model.mean,'lb')==0 and getattr(model.mean,'ub')==10 
#================================================================

def test_call_keywords():
#================================================================
    "Check that calling the model with parameter returns the correct response"
    model = Model(gauss)
    response  = model(mean=3,width=0.5)
    reference = gauss(mean=3,width=0.5)

    assert np.allclose(response,reference)
#================================================================

def test_call_positional():
#================================================================
    "Check that calling the model with parameter returns the correct response"
    model = Model(gauss)
    response  = model(3,0.5)
    reference = gauss(3,0.5)
    
    assert np.allclose(response,reference)
#================================================================


def test_addlinear_length():
#================================================================
    "Check that the model is contructed correctly with the appropiate number of parameters"
    model = Model(gauss2_design)
    model.addlinear('amp1',lb=0)
    model.addlinear('amp2',lb=0)

    assert model.Nparam==6
#================================================================

def test_addlinear_names():
#================================================================
    "Check that linear parameters can be properly added"
    model = Model(gauss2_design)
    model.addlinear('amp1',lb=0)
    model.addlinear('amp2',lb=0)

    assert 'amp1' in model.__dict__ and 'amp2' in model.__dict__
#================================================================

def test_addlinear_set():
#================================================================
    "Check that attributes of the linear parameters are editable"
    model = Model(gauss2_design)
    model.addlinear('amp1')
    model.amp1.set(lb=0,ub=10)

    assert getattr(model.amp1,'lb')==0 and getattr(model.amp1,'ub')==10 
#================================================================

def test_addlinear_call_keywords():
#================================================================
    "Check that calling the model with parameters returns the correct response"
    model = Model(gauss2_design)
    model.addlinear('amp1')
    model.addlinear('amp2')
    response = model(mean1=3,mean2=4,width1=0.2,width2=0.3,amp1=0.5,amp2=0.4)
    reference = gauss2(mean1=3,mean2=4,width1=0.2,width2=0.3,amp1=0.5,amp2=0.4)

    assert np.allclose(response,reference)
#================================================================

def test_addlinear_call_positional():
#================================================================
    "Check that calling the model with parameters returns the correct response"
    model = Model(gauss2_design)
    model.addlinear('amp1')
    model.addlinear('amp2')
    response = model(3,4,0.2,0.3,0.5,0.4)
    reference = gauss2(3,4,0.2,0.3,0.5,0.4)

    assert np.allclose(response,reference)
#================================================================

#================================================================
def test_addlinear_vector_length():
    "Check that linear parameters can be defined as vectors"
    model = Model(gauss2_identity)
    model.addlinear('gaussian', vec=100)

    assert model.Nparam==100
#================================================================

#================================================================
def test_addlinear_vector_names():
    "Check that linear parameters can be defined as vectors"
    model = Model(gauss2_identity)
    model.addlinear('gaussian', vec=100)

    assert 'gaussian' in model.__dict__ 
#================================================================

def test_addlinear_vector_set():
#================================================================
    "Check that attributes of the vector linear parameters are editable"
    model = Model(gauss2_identity)
    model.addlinear('gaussian', vec=100)
    model.gaussian.set(lb=np.zeros(100))

    assert np.allclose(getattr(model.gaussian,'lb'),np.zeros(100)) 
#================================================================

def test_addlinear_vector_call_keywords():
#================================================================
    "Check that calling the model with scalar and vector parameters returns the correct response"
    model = Model(gauss2_identity)
    model.addlinear('gaussian', vec=len(x))
    reference = gauss2(mean1=3,mean2=4,width1=0.2,width2=0.3,amp1=0.5,amp2=0.4)
    response = model(gaussian=reference)

    assert np.allclose(response,reference)
#================================================================

def test_addlinear_vector_call_positional():
#================================================================
    "Check that calling the model with scalar and vector parameters returns the correct response"
    model = Model(gauss2_identity)
    model.addlinear('gaussian', vec=len(x))
    reference = gauss2(3,4,0.2,0.3,0.5,0.4)
    response = model(reference)

    assert np.allclose(response,reference)
#================================================================


def test_mixed_vector_length():
#================================================================
    "Check the definition of scalar nonlinear parameters and vector linear parameters"
    model = Model(gauss2_scaled)
    model.addlinear('gaussian', vec=100)

    assert model.Nparam==101
#================================================================

def test_mixed_vector_names():
#================================================================
    "Check the definition of scalar nonlinear parameters and vector linear parameters"
    model = Model(gauss2_scaled)
    model.addlinear('gaussian', vec=100)

    assert 'gaussian' in model.__dict__  and 'scale' in model.__dict__ 
#================================================================

def test_mixed_vector_call_keywords():
#================================================================
    "Check that calling the model with scalar and vector parameters returns the correct response"
    model = Model(gauss2_scaled)
    model.addlinear('gaussian', vec=len(x))
    reference = 5*gauss2(mean1=3,mean2=4,width1=0.2,width2=0.3,amp1=0.5,amp2=0.4)
    response = model(scale=5,gaussian=gauss2(mean1=3,mean2=4,width1=0.2,width2=0.3,amp1=0.5,amp2=0.4))

    assert np.allclose(response,reference)
#================================================================

def test_mixed_vector_call_positional():
#================================================================
    "Check that calling the model with scalar and vector parameters returns the correct response"
    model = Model(gauss2_scaled)
    model.addlinear('gaussian', vec=len(x))
    reference = 5*gauss2(3,4,0.2,0.3,0.5,0.4)
    response = model(5,gauss2(3,4,0.2,0.3,0.5,0.4))

    assert np.allclose(response,reference)
#================================================================


def test_addnonlinear_length():
#================================================================
    "Check that the model is contructed correctly with the appropiate number of parameters"
    model = Model(gauss)
    model.addnonlinear('trivial1',lb=0)
    model.addnonlinear('trivial2',lb=0)

    assert model.Nparam==4
#================================================================

def test_addnonlinear_names():
#================================================================
    "Check that linear parameters can be properly added"
    model = Model(gauss)
    model.addnonlinear('trivial1',lb=0)
    model.addnonlinear('trivial2',lb=0)

    assert hasattr(model,'trivial1') and hasattr(model,'trivial2')
#================================================================

def test_addnonlinear_set():
#================================================================
    "Check that attributes of the linear parameters are editable"
    model = Model(gauss)
    model.addnonlinear('trivial1')
    model.trivial1.set(lb=0,ub=10)

    assert getattr(model.trivial1,'lb')==0 and getattr(model.trivial1,'ub')==10 
#================================================================

def test_addnonlinear_call_keywords():
#================================================================
    "Check that calling the model with parameters returns the correct response"
    model = Model(gauss)
    model.addnonlinear('trivial1')
    model.addnonlinear('trivial2')
    response = model(mean=3,width=0.2,trivial1=1,trivial2=1)
    reference = gauss(mean=3,width=0.2)

    assert np.allclose(response,reference)
#================================================================

def test_addnonlinear_call_positional():
#================================================================
    "Check that calling the model with parameters returns the correct response"
    model = Model(gauss)
    model.addnonlinear('trivial1')
    model.addnonlinear('trivial2')
    response = model(3,0.2,1,1)
    reference = gauss(3,0.2)

    assert np.allclose(response,reference)
#================================================================

#----------------------------------------------------------------
def _getmodel(type):
    if type=='parametric':
        model = Model(gauss2)
        model.mean1.set(lb=0, ub=10, par0=2)
        model.mean2.set(lb=0, ub=10, par0=4)
        model.width1.set(lb=0.01, ub=5, par0=0.2)
        model.width2.set(lb=0.01, ub=5, par0=0.2)
        model.amp1.set(lb=0, ub=5, par0=1)
        model.amp2.set(lb=0, ub=5, par0=1)
    elif type=='semiparametric': 
        model = Model(gauss2_design)
        model.mean1.set(lb=0, ub=10, par0=2)
        model.mean2.set(lb=0, ub=10, par0=4)
        model.width1.set(lb=0.01, ub=5, par0=0.2)
        model.width2.set(lb=0.01, ub=5, par0=0.2)
        model.addlinear('amp1',lb=0, ub=5)
        model.addlinear('amp2',lb=0, ub=5)
    elif type=='nonparametric':
        model = Model(gauss2_design(3,4,0.5,0.2))
        model.addlinear('amp1',lb=0)
        model.addlinear('amp2',lb=0)
    return model
#----------------------------------------------------------------

def test_fit_parametric(): 
#================================================================
    "Check that a parametric model can be correctly fitted"
    model = _getmodel('parametric')

    fitResult = fit(model,mock_data)
    
    assert np.allclose(fitResult.model,mock_data)
#================================================================

def test_fit_semiparametric(): 
#================================================================
    "Check that a semiparametric model can be correctly fitted"
    model = _getmodel('semiparametric')

    fitResult = fit(model,mock_data)
    
    assert np.allclose(fitResult.model,mock_data)
#================================================================

def test_fit_nonparametric(): 
#================================================================
    "Check that a semiparametric model can be correctly fitted"
    model = _getmodel('nonparametric')

    fitResult = fit(model,mock_data)
    
    assert np.allclose(fitResult.model,mock_data,atol=1e-3)
#================================================================

def test_freeze():
#================================================================
    "Check that a model parameter can be frozen to a fixed value"
    model = Model(gauss)
    model.mean.freeze(3)

    assert model.mean.value==3 and model.mean.frozen==True
#================================================================

def test_unfreeze():
#================================================================
    "Check that a model parameter can be frozen and then reversed"
    model = Model(gauss)
    model.mean.freeze(3)
    model.mean.unfreeze()

    assert model.mean.frozen==False and model.mean.value==None
#================================================================


def test_fit_parametric_frozen(): 
#================================================================
    "Check that a parametric model can be correctly fitted"
    model = _getmodel('parametric')

    model.mean1.freeze(3)

    fitResult = fit(model,mock_data)
    
    assert np.allclose(fitResult.model,mock_data)
#================================================================

def test_fit_semiparametric_frozen(): 
#================================================================
    "Check that a semiparametric model can be correctly fitted"
    model = _getmodel('semiparametric')

    model.mean1.freeze(3)

    fitResult = fit(model,mock_data)
    
    assert np.allclose(fitResult.model,mock_data)
#================================================================

def test_fit_nonparametric_frozen(): 
#================================================================
    "Check that a semiparametric model can be correctly fitted"
    model = _getmodel('nonparametric')

    model.amp1.freeze(0.5)
    fitResult = fit(model,mock_data)
    
    assert np.allclose(fitResult.model,mock_data)
#================================================================

#----------------------------------------------------------------
def assert_attributes_cis(fitobject,attributes):
    for attr in attributes: 
        parfit = getattr(fitobject,attr)
        parci = getattr(fitobject,f'{attr}Uncert').ci(95)
        ci_lower = parci[0]
        ci_upper = parci[1]
        if getattr(fitobject,f'{attr}Uncert').type=='bootstrap':
            assert np.allclose(parfit,getattr(fitobject,f'{attr}Uncert').median)
        assert parfit<=ci_upper and parfit>=ci_lower 
#----------------------------------------------------------------

def test_CIs_parametric(): 
#================================================================
    "Check the default confidence intervals of the fitted parameters"
    model = _getmodel('parametric')

    fitResult = fit(model,mock_data)
    
    assert_attributes_cis(fitResult,['mean1','mean2','width1','width2','amp1','amp2'])
#================================================================

def test_CIs_semiparametric(): 
#================================================================
    "Check the default confidence intervals of the fitted parameters"
    model = _getmodel('semiparametric')

    fitResult = fit(model,mock_data)
    
    assert_attributes_cis(fitResult,['mean1','mean2','width1','width2','amp1','amp2'])
#================================================================

def test_CIs_nonparametric(): 
#================================================================
    "Check the default confidence intervals of the fitted parameters"
    model = _getmodel('semiparametric')

    fitResult = fit(model,mock_data)
    
    assert_attributes_cis(fitResult,['amp1','amp2'])
#================================================================


def test_bootCIs_parametric(): 
#================================================================
    "Check the bootstrapped confidence intervals of the fitted parameters"
    model = _getmodel('parametric')

    noisydata = mock_data + whitegaussnoise(0.01,seed=1)
    fitResult = fit(model,noisydata,bootstrap=3)
    
    assert_attributes_cis(fitResult,['mean1','mean2','width1','width2','amp1','amp2'])
#================================================================

def test_bootCIs_semiparametric(): 
#================================================================
    "Check the bootstrapped confidence intervals of the fitted parameters"
    model = _getmodel('semiparametric')

    noisydata = mock_data + whitegaussnoise(0.01,seed=1)
    fitResult = fit(model,noisydata,bootstrap=3)
    
    assert_attributes_cis(fitResult,['mean1','mean2','width1','width2','amp1','amp2'])
#================================================================

def test_bootCIs_nonparametric(): 
#================================================================
    "Check the bootstrapped confidence intervals of the fitted parameters"
    model = _getmodel('semiparametric')

    noisydata = mock_data + whitegaussnoise(0.01,seed=1)
    fitResult = fit(model,noisydata,bootstrap=3)
    
    assert_attributes_cis(fitResult,['amp1','amp2'])
#================================================================

# Simple non-linear function for testing
def gauss_axis(axis,mean,width): 
    return np.exp(-(axis-mean)**2/width**2/2)

# Non-linear definition
def gauss2_axis(axis,mean1,mean2,width1,width2,amp1,amp2):
    return amp1*gauss_axis(axis,mean1,width1) + amp2*gauss_axis(axis,mean2,width2)

# Linear + Non-linear definition
def gauss2_design_axis(axis,mean1,mean2,width1,width2):
    return np.atleast_2d([gauss_axis(axis,mean1,width1), gauss_axis(axis,mean2,width2)]).T

# Linear + Non-linear definition
def gauss2_identity_axis(axis):
    return np.eye(len(axis))

mock_data_fcn = lambda axis: gauss2_axis(axis,mean1=3,mean2=4,width1=0.5,width2=0.2,amp1=0.5,amp2=0.6)


def test_model_with_constant_positional(): 
#================================================================
    "Check that a model with axis can be defined and called"
    model = Model(gauss_axis,constants='axis')

    x = np.linspace(0,10,300)
    reference = gauss_axis(x,3,0.5)
    response = model(x,3,0.5)
        
    assert np.allclose(reference,response)
#================================================================

def test_model_with_constant_keywords(): 
#================================================================
    "Check that a model with axis can be defined and called via keywords"
    model = Model(gauss_axis,constants='axis')

    x = np.linspace(0,10,300)
    reference = gauss_axis(x,3,0.5)
    response = model(axis=x, mean=3, width=0.5)
        
    assert np.allclose(reference,response)
#================================================================


#----------------------------------------------------------------
def _getmodel_axis(type):
    if type=='parametric':
        model = Model(gauss2_axis,constants='axis')
        model.mean1.set(lb=0, ub=10, par0=2)
        model.mean2.set(lb=0, ub=10, par0=4)
        model.width1.set(lb=0.01, ub=5, par0=0.2)
        model.width2.set(lb=0.01, ub=5, par0=0.2)
        model.amp1.set(lb=0, ub=5, par0=1)
        model.amp2.set(lb=0, ub=5, par0=1)
    elif type=='semiparametric': 
        model = Model(gauss2_design_axis,constants='axis')
        model.mean1.set(lb=0, ub=10, par0=2)
        model.mean2.set(lb=0, ub=10, par0=4)
        model.width1.set(lb=0.01, ub=5, par0=0.2)
        model.width2.set(lb=0.01, ub=5, par0=0.2)
        model.addlinear('amp1',lb=0, ub=5)
        model.addlinear('amp2',lb=0, ub=5)
    elif type=='nonparametric':
        model = Model(lambda x: gauss2_design_axis(x,3,4,0.5,0.2),constants='x')
        model.addlinear('amp1',lb=0)
        model.addlinear('amp2',lb=0)
    return model
#----------------------------------------------------------------

def test_fit_parametric_constant(): 
#================================================================
    "Check that a parametric model can be correctly fitted while specifying an axis"
    model = _getmodel_axis('parametric')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x)
    
    assert np.allclose(fitResult.model,mock_data_fcn(x))
#================================================================

def test_fit_semiparametric_constant(): 
#================================================================
    "Check that a semiparametric model can be correctly fitted while specifying an axis"
    model = _getmodel_axis('semiparametric')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x)
    
    assert np.allclose(fitResult.model,mock_data_fcn(x))
#================================================================

def test_fit_nonparametric_constant(): 
#================================================================
    "Check that a semiparametric model can be correctly fitted while specifying an axis"
    model = _getmodel_axis('nonparametric')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x)
    
    assert np.allclose(fitResult.model,mock_data_fcn(x),atol=1e-3)
#================================================================


def gauss_multiaxis(axis1,axis2,mean,width): 
    return np.exp(-(axis1-mean)**2/width**2/2)

def gauss2_multiaxis(axis1,axis2,mean1,mean2,width1,width2,amp1,amp2):
    return amp1*gauss_axis(axis1,mean1,width1) + amp2*gauss_axis(axis2,mean2,width2)


def test_model_with_multiple_constants(): 
#================================================================
    "Check that a model with axis can be defined and called"
    model = Model(gauss_multiaxis,constants=['axis1','axis2'])

    x1 = np.linspace(0,5,300)
    x2 = np.linspace(5,10,300)
    reference = gauss_multiaxis(x1,x2,3,0.5)
    response = model(x1,x2,3,0.5)
        
    assert np.allclose(reference,response)
#================================================================


def test_model_with_multiple_constants_fit(): 
#================================================================
    "Check that a model with axis can be defined and called"
    model = Model(gauss2_multiaxis,constants=['axis1','axis2'])
    model.mean1.set(lb=0, ub=10, par0=2)
    model.mean2.set(lb=0, ub=10, par0=4)
    model.width1.set(lb=0.01, ub=5, par0=0.2)
    model.width2.set(lb=0.01, ub=5, par0=0.2)
    model.amp1.set(lb=0, ub=5, par0=1)
    model.amp2.set(lb=0, ub=5, par0=1)
    
    x = np.linspace(0,10,300)
    x1 = np.linspace(0,10,300)
    x2 = np.linspace(0,10,300)
    fitResult = fit(model,mock_data_fcn(x),x1,x2)
    assert np.allclose(fitResult.model,mock_data_fcn(x))        
#================================================================

def test_model_constant_same_values_keywords(): 
#================================================================
    "Check that a model with axis can be defined and called with same values as parameters"
    model = Model(gauss_axis,constants='axis')

    reference = gauss_axis(axis=3,mean=3,width=0.5)
    response = model(axis=3,mean=3,width=0.5)
        
    assert np.allclose(reference,response)
#================================================================

def test_model_constant_same_values_positional(): 
#================================================================
    "Check that a model with axis can be defined and called with same values as parameters"
    model = Model(gauss_axis,constants='axis')

    reference = gauss_axis(3,3,0.5)
    response = model(3,3,0.5)
        
    assert np.allclose(reference,response)
#================================================================