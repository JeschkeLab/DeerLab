from collections import namedtuple
from deerlab.whitegaussnoise import whitegaussnoise
from deerlab.model import Model, fit
import numpy as np 
import pytest
import os 
from deerlab.utils import store_pickle, read_pickle

# Simple non-linear function for testing
x = np.linspace(0,5,100)
def gauss(mean,std): 
    return np.exp(-(x-mean)**2/std**2/2)

# Non-linear definition
def gauss2(mean1,mean2,std1,std2,amp1,amp2):
    return amp1*gauss(mean1,std1) + amp2*gauss(mean2,std2)

# Linear + Non-linear definition
def gauss2_design(mean1,mean2,std1,std2):
    return np.atleast_2d([gauss(mean1,std1), gauss(mean2,std2)]).T

# Linear + Non-linear definition
def gauss2_identity():
    return np.eye(len(x))

# Linear + Non-linear definition
def gauss2_scaled(scale):
    return scale*np.eye(len(x))

mock_data = gauss2(mean1=3,mean2=4,std1=0.5,std2=0.2,amp1=0.5,amp2=0.6)


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

    assert 'mean' in model.__dict__ and 'std' in model.__dict__
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
    response  = model(mean=3,std=0.5)
    reference = gauss(mean=3,std=0.5)

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

def test_call_mixed():
#================================================================
    "Check that calling the model with parameter returns the correct response"
    model = Model(gauss)
    response  = model(3,std=0.5)
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

def test_addlinear_normalization():
#================================================================
    "Check that linear parameters can be properly added with a normalization function"
    model = Model(gauss2_design)
    model.addlinear('amp1',lb=0, normalization=lambda x: 2*x)
    model.addlinear('amp2',lb=0, normalization=lambda x: 3*x)

    assert 'amp1' in model.__dict__ and 'amp2' in model.__dict__
#================================================================

def test_addlinear_call_keywords():
#================================================================
    "Check that calling the model with parameters returns the correct response"
    model = Model(gauss2_design)
    model.addlinear('amp1')
    model.addlinear('amp2')
    response = model(mean1=3,mean2=4,std1=0.2,std2=0.3,amp1=0.5,amp2=0.4)
    reference = gauss2(mean1=3,mean2=4,std1=0.2,std2=0.3,amp1=0.5,amp2=0.4)

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


def test_addlinear_call_mixed():
#================================================================
    "Check that calling the model with parameters returns the correct response"
    model = Model(gauss2_design)
    model.addlinear('amp1')
    model.addlinear('amp2')
    response = model(3,4,0.2,std2=0.3,amp1=0.5,amp2=0.4)
    reference = gauss2(3,4,0.2,std2=0.3,amp1=0.5,amp2=0.4)

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
    reference = gauss2(mean1=3,mean2=4,std1=0.2,std2=0.3,amp1=0.5,amp2=0.4)
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
    reference = 5*gauss2(mean1=3,mean2=4,std1=0.2,std2=0.3,amp1=0.5,amp2=0.4)
    response = model(scale=5,gaussian=gauss2(mean1=3,mean2=4,std1=0.2,std2=0.3,amp1=0.5,amp2=0.4))

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
    response = model(mean=3,std=0.2,trivial1=1,trivial2=1)
    reference = gauss(mean=3,std=0.2)

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
        model.std1.set(lb=0.1, ub=5, par0=0.2)
        model.std2.set(lb=0.1, ub=5, par0=0.2)
        model.amp1.set(lb=0, ub=5, par0=1)
        model.amp2.set(lb=0, ub=5, par0=1)
    elif type=='semiparametric': 
        model = Model(gauss2_design)
        model.mean1.set(lb=0, ub=10, par0=2)
        model.mean2.set(lb=0, ub=10, par0=4)
        model.std1.set(lb=0.1, ub=5, par0=0.2)
        model.std2.set(lb=0.1, ub=5, par0=0.2)
        model.addlinear('amp1',lb=0, ub=5)
        model.addlinear('amp2',lb=0, ub=5)
    elif type=='nonparametric':
        model = Model(gauss2_design(3,4,0.5,0.2))
        model.addlinear('amp1',lb=0)
        model.addlinear('amp2',lb=0)
    return model
#----------------------------------------------------------------

# ======================================================================
def test_preserve_original(): 
    "Check that the original model is not changed by the function"
    model = Model(gauss)
    model.mean.par0 = 3
    model.std.par0 = 0.2
    
    _ = fit(model,mock_data)
    assert model._parameter_list() == ['mean','std']
# ======================================================================


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

def test_freeze_outofbounds():
#================================================================
    "Check that a parameter cannot be frozen outside of the bounds"
    model = Model(gauss)
    model.mean.set(lb=0,ub=10)
    with pytest.raises(ValueError):
        model.mean.freeze(-10)
#================================================================

def test_freeze_vec_outofbounds():
#================================================================
    "Check that a parameter cannot be frozen outside of the bounds"
    model = Model(gauss2_identity)
    model.addlinear('gaussian', vec=100, lb=0)
    with pytest.raises(ValueError):
        model.gaussian.freeze(np.full(100,-10))
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

def test_fit_fullyfrozen_linear(): 
#================================================================
    "Check that a model with all linear parameters frozen can be fitted"

    model = _getmodel('semiparametric')
    model.amp1.freeze(0.5)
    model.amp2.freeze(0.6)

    fitResult = fit(model,mock_data)
    
    assert np.allclose(fitResult.model,mock_data)
#================================================================

def test_fit_fullyfrozen_nonlinear(): 
#================================================================
    "Check that a model with all nonlinear parameters frozen can be fitted"

    model = _getmodel('semiparametric')
    model.mean1.freeze(3)
    model.mean2.freeze(4)
    model.std1.freeze(0.5)
    model.std2.freeze(0.2)

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
    
    assert_attributes_cis(fitResult,['mean1','mean2','std1','std2','amp1','amp2'])
#================================================================

def test_CIs_semiparametric(): 
#================================================================
    "Check the default confidence intervals of the fitted parameters"
    model = _getmodel('semiparametric')

    fitResult = fit(model,mock_data)
    
    assert_attributes_cis(fitResult,['mean1','mean2','std1','std2','amp1','amp2'])
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

    noisydata = mock_data + whitegaussnoise(x,0.01,seed=1)
    fitResult = fit(model,noisydata,bootstrap=3)
    
    assert_attributes_cis(fitResult,['mean1','mean2','std1','std2','amp1','amp2'])
#================================================================

def test_bootCIs_semiparametric(): 
#================================================================
    "Check the bootstrapped confidence intervals of the fitted parameters"
    model = _getmodel('semiparametric')

    noisydata = mock_data + whitegaussnoise(x,0.01,seed=1)
    fitResult = fit(model,noisydata,bootstrap=3)
    
    assert_attributes_cis(fitResult,['mean1','mean2','std1','std2','amp1','amp2'])
#================================================================

def test_bootCIs_nonparametric(): 
#================================================================
    "Check the bootstrapped confidence intervals of the fitted parameters"
    model = _getmodel('semiparametric')

    noisydata = mock_data + whitegaussnoise(x,0.01,seed=1)
    fitResult = fit(model,noisydata,bootstrap=3)
    
    assert_attributes_cis(fitResult,['amp1','amp2'])
#================================================================

# Simple non-linear function for testing
def gauss_axis(axis,mean,std): 
    return np.exp(-(axis-mean)**2/std**2/2)

# Non-linear definition
def gauss2_axis(axis,mean1,mean2,std1,std2,amp1,amp2):
    return amp1*gauss_axis(axis,mean1,std1) + amp2*gauss_axis(axis,mean2,std2)

# Linear + Non-linear definition
def gauss2_design_axis(axis,mean1,mean2,std1,std2):
    return np.atleast_2d([gauss_axis(axis,mean1,std1), gauss_axis(axis,mean2,std2)]).T

# Linear + Non-linear definition
def gauss2_identity_axis(axis):
    return np.eye(len(axis))

mock_data_fcn = lambda axis: gauss2_axis(axis,mean1=3,mean2=4,std1=0.5,std2=0.2,amp1=0.5,amp2=0.6)


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
    response = model(axis=x, mean=3, std=0.5)
        
    assert np.allclose(reference,response)
#================================================================

def test_model_with_constant_mixed(): 
#================================================================
    "Check that a model with axis can be defined and called via keywords"
    model = Model(gauss_axis,constants='axis')

    x = np.linspace(0,10,300)
    reference = gauss_axis(x,3,0.5)
    response = model(x, mean=3, std=0.5)
        
    assert np.allclose(reference,response)
#================================================================

#----------------------------------------------------------------
def _getmodel_axis(type,vec=50):
    if type=='parametric':
        model = Model(gauss2_axis,constants='axis')
        model.mean1.set(lb=0, ub=10, par0=2)
        model.mean2.set(lb=0, ub=10, par0=4)
        model.std1.set(lb=0.01, ub=5, par0=0.2)
        model.std2.set(lb=0.01, ub=5, par0=0.2)
        model.amp1.set(lb=0, ub=5, par0=1)
        model.amp2.set(lb=0, ub=5, par0=1)
    elif type=='semiparametric': 
        model = Model(gauss2_design_axis,constants='axis')
        model.mean1.set(lb=0, ub=10, par0=2)
        model.mean2.set(lb=0, ub=10, par0=4)
        model.std1.set(lb=0.01, ub=5, par0=0.2)
        model.std2.set(lb=0.01, ub=5, par0=0.2)
        model.addlinear('amp1',lb=0, ub=5)
        model.addlinear('amp2',lb=0, ub=5)
    elif type=='semiparametric_vec': 
        model = Model(gauss2_design_axis,constants='axis')
        model.mean1.set(lb=0, ub=10, par0=2)
        model.mean2.set(lb=0, ub=10, par0=4)
        model.std1.set(lb=0.01, ub=5, par0=0.2)
        model.std2.set(lb=0.01, ub=5, par0=0.2)
        model.addlinear('amps',lb=0, ub=5, vec=2)
    elif type=='nonparametric':
        model = Model(lambda x: gauss2_design_axis(x,3,4,0.5,0.2),constants='x')
        model.addlinear('amp1',lb=0)
        model.addlinear('amp2',lb=0)
    elif type=='nonparametric_vec':
        model = Model(lambda x: np.eye(len(x)),constants='x')
        model.addlinear('dist',lb=0,vec=vec)
    elif type=='nonparametric_vec_normalized':
        model = Model(lambda x: np.eye(len(x)),constants='x')
        model.addlinear('dist',lb=0,vec=vec,normalization=lambda dist:dist/np.sum(dist))
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
    "Check that a vectorized nonparametric model can be correctly fitted while specifying an axis"
    model = _getmodel_axis('nonparametric')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x)
    
    assert np.allclose(fitResult.model,mock_data_fcn(x),atol=1e-3)
#================================================================

def test_fit_nonparametric_vec_constant(): 
#================================================================
    "Check that a vectorized nonparametric model can be correctly fitted while specifying an axis"
    model = _getmodel_axis('nonparametric_vec',vec=80)

    x = np.linspace(0,10,80)
    fitResult = fit(model,mock_data_fcn(x),x,noiselvl=1e-10)
    
    assert np.allclose(fitResult.model,mock_data_fcn(x),atol=1e-3)
#================================================================

def test_fit_nonparametric_vec_normalized_constant(): 
#================================================================
    "Check that a vectorized nonparametric model can be correctly fitted while specifying an axis"
    model = _getmodel_axis('nonparametric_vec_normalized',vec=80)

    x = np.linspace(0,10,80)
    fitResult = fit(model,mock_data_fcn(x),x,noiselvl=1e-10)
    
    assert np.allclose(fitResult.model,mock_data_fcn(x),atol=1e-3)
#================================================================

def test_fit_semiparametric_vec_constant(): 
#================================================================
    "Check that a vectorized semiparametric model can be correctly fitted while specifying an axis"
    model = _getmodel_axis('semiparametric_vec')

    x = np.linspace(0,10,100)
    fitResult = fit(model,mock_data_fcn(x),x,noiselvl=1e-10)
    
    assert np.allclose(fitResult.model,mock_data_fcn(x),atol=1e-3)
#================================================================

def gauss_multiaxis(axis1,axis2,mean,std): 
    return np.exp(-(axis1-mean)**2/std**2/2)

def gauss2_multiaxis(axis1,axis2,mean1,mean2,std1,std2,amp1,amp2):
    return amp1*gauss_axis(axis1,mean1,std1) + amp2*gauss_axis(axis2,mean2,std2)


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
    model.std1.set(lb=0.01, ub=5, par0=0.2)
    model.std2.set(lb=0.01, ub=5, par0=0.2)
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

    reference = gauss_axis(axis=3,mean=3,std=0.5)
    response = model(axis=3,mean=3,std=0.5)
        
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


def model(phase, center, std):
    y = gauss(center, std)
    y = y*np.exp(-1j*phase)
    return y
mymodel = Model(model)
mymodel.phase.set(par0=2*np.pi/5, lb=-np.pi, ub=np.pi)
mymodel.center.set(par0=4, lb=1, ub=6)
mymodel.std.set(par0=0.2, lb=0.05, ub=5)

y = mymodel(phase=np.pi/5, center=3, std=0.5)     

def test_complex_model_complex_data():
# ======================================================================
    "Check the fit of a complex-valued model to complex-valued data"
    fitResult = fit(mymodel,y)

    assert np.allclose(fitResult.model.real,y.real) and np.allclose(fitResult.model.imag,y.imag)
# ======================================================================

def test_complex_model_real_data():
# ======================================================================
    "Check the fit of a complex-valued model to real-valued data"
    fitResult = fit(mymodel,y.real)

    assert np.allclose(fitResult.model.real,y.real) and np.allclose(fitResult.model.imag,np.zeros_like(y.real))
# ======================================================================

def test_real_model_complex_data():
# ======================================================================
    "Check the fit of a real-valued model to complex-valued data"
    mymodel = Model(gauss)
    mymodel.mean.set(par0=4, lb=1, ub=6)
    mymodel.std.set(par0=0.2, lb=0.05, ub=5)
    fitResult = fit(mymodel,y.real + 1j*np.zeros_like(y.imag))

    assert np.allclose(fitResult.model.real,y.real) and np.allclose(fitResult.model.imag,np.zeros_like(y.real))
# ======================================================================

def test_fit_evaluate_parametric(): 
#================================================================
    "Check the evaluate method of the fitResult object for evaluation of a parametric model"
    model = _getmodel_axis('parametric')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x)
    
    response = fitResult.scale*fitResult.evaluate(model,x)

    assert np.allclose(response,mock_data_fcn(x))
#================================================================

def test_fit_evaluate_nonparametric(): 
#================================================================
    "Check the evaluate method of the fitResult object for evaluation of a nonparametric model"
    model = _getmodel_axis('nonparametric')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x)
    
    response = fitResult.evaluate(model,x)

    assert np.allclose(response,mock_data_fcn(x))
#================================================================

def test_fit_evaluate_nonparametric_vec(): 
#================================================================
    "Check the evaluate method of the fitResult object for evaluation of a vectorized nonparametric model"
    model = _getmodel_axis('nonparametric_vec',vec=80)

    x = np.linspace(0,10,80)
    fitResult = fit(model,mock_data_fcn(x),x)
    
    response = fitResult.evaluate(model,x)

    assert np.allclose(response,mock_data_fcn(x))
#================================================================

def test_fit_evaluate_nonparametric_vec_normalized(): 
#================================================================
    "Check the evaluate method of the fitResult object for evaluation of a vectorized nonparametric model"
    model = _getmodel_axis('nonparametric_vec_normalized',vec=80)

    x = np.linspace(0,10,80)
    fitResult = fit(model,mock_data_fcn(x),x)
    
    response = fitResult.evaluate(model,x)

    assert np.allclose(response,mock_data_fcn(x))
#================================================================

def test_fit_evaluate_semiparametric(): 
#================================================================
    "Check the evaluate method of the fitResult object for evaluation of a semiparametric model"
    model = _getmodel_axis('semiparametric')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x)
    
    response = fitResult.evaluate(model,x)

    assert np.allclose(response,mock_data_fcn(x))
#================================================================

def test_fit_evaluate_semiparametric_vec(): 
#================================================================
    "Check the evaluate method of the fitResult object for evaluation of a vectorized semiparametric model"
    model = _getmodel_axis('semiparametric_vec')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x)
    
    response = fitResult.evaluate(model,x)

    assert np.allclose(response,mock_data_fcn(x))
#================================================================

def test_fit_evaluate_callable(): 
#================================================================
    "Check the evaluate method of the fitResult object for evaluation of a callable function"
    model = _getmodel_axis('semiparametric_vec')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x)
    
    fcn = lambda mean1,mean2,std1,std2,amps: gauss2_design_axis(x,mean1,mean2,std1,std2)@amps
    response = fitResult.evaluate(fcn)

    assert np.allclose(response,mock_data_fcn(x))
#================================================================

#----------------------------------------------------------------
def assert_cis(argUncert):
    parci = argUncert.ci(95)
    ci_lower = parci[:,0]
    ci_upper = parci[:,1]

    assert np.less_equal(ci_lower,ci_upper).all()
#----------------------------------------------------------------


def test_fit_propagate_parametric(): 
#================================================================
    "Check the propagate method of the fitResult object for evaluation of a parametric model"
    model = _getmodel_axis('parametric')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x)
    
    modeluq = fitResult.propagate(model, x, lb=np.zeros_like(x))

    assert_cis(modeluq)
#================================================================

def test_fit_propagate_nonparametric(): 
#================================================================
    "Check the propagate method of the fitResult object for evaluation of a nonparametric model"
    model = _getmodel_axis('nonparametric')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x)
    
    modeluq = fitResult.propagate(model, x, lb=np.zeros_like(x))

    assert_cis(modeluq)
#================================================================

def test_fit_propagate_nonparametric_vec(): 
#================================================================
    "Check the propagate method of the fitResult object for evaluation of a vectorized nonparametric model"
    model = _getmodel_axis('nonparametric_vec',vec=80)

    x = np.linspace(0,10,80)
    fitResult = fit(model,mock_data_fcn(x),x)
    
    modeluq = fitResult.propagate(model, x, lb=np.zeros_like(x))

    assert_cis(modeluq)
#================================================================

def test_fit_propagate_nonparametric_vec_normalized(): 
#================================================================
    "Check the propagate method of the fitResult object for evaluation of a vectorized nonparametric model"
    model = _getmodel_axis('nonparametric_vec_normalized',vec=80)

    x = np.linspace(0,10,80)
    fitResult = fit(model,mock_data_fcn(x),x)
    
    modeluq = fitResult.propagate(model, x, lb=np.zeros_like(x))

    assert_cis(modeluq)
#================================================================

def test_fit_propagate_semiparametric(): 
#================================================================
    "Check the propagate method of the fitResult object for evaluation of a semiparametric model"
    model = _getmodel_axis('semiparametric')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x)
    
    modeluq = fitResult.propagate(model, x, lb=np.zeros_like(x))

    assert_cis(modeluq)
#================================================================

def test_fit_propagate_semiparametric_vec(): 
#================================================================
    "Check the propagate method of the fitResult object for evaluation of a vectorized semiparametric model"
    model = _getmodel_axis('semiparametric_vec')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x)
    
    modeluq = fitResult.propagate(model, x, lb=np.zeros_like(x))

    assert_cis(modeluq)
#================================================================

def test_fit_propagate_callable(): 
#================================================================
    "Check the evaluate method of the fitResult object for evaluation of a callable function"
    model = _getmodel_axis('semiparametric_vec')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x)
    
    fcn = lambda mean1,mean2,std1,std2,amps: gauss2_design_axis(x,mean1,mean2,std1,std2)@amps
    modeluq = fitResult.propagate(fcn, lb=np.zeros_like(x))

    assert_cis(modeluq)
#================================================================

def test_fit_propagate_parametric_bootstrapped(): 
#================================================================
    "Check the propagate method of the fitResult object for evaluation of a parametric model"
    model = _getmodel_axis('parametric')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x,bootstrap=3)
    
    modeluq = fitResult.propagate(model, x, lb=np.zeros_like(x))

    assert_cis(modeluq)
#================================================================

def test_fit_propagate_nonparametric_bootstrapped(): 
#================================================================
    "Check the propagate method of the fitResult object for evaluation of a nonparametric model"
    model = _getmodel_axis('nonparametric')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x,bootstrap=3)
    
    modeluq = fitResult.propagate(model, x, lb=np.zeros_like(x))

    assert_cis(modeluq)
#================================================================

def test_fit_propagate_nonparametric_vec_bootstrapped(): 
#================================================================
    "Check the propagate method of the fitResult object for evaluation of a vectorized nonparametric model"
    model = _getmodel_axis('nonparametric_vec',vec=80)

    x = np.linspace(0,10,80)
    fitResult = fit(model,mock_data_fcn(x),x,bootstrap=3)
    
    modeluq = fitResult.propagate(model, x, lb=np.zeros_like(x))

    assert_cis(modeluq)
#================================================================

def test_fit_propagate_nonparametric_vec_normalized_bootstrapped(): 
#================================================================
    "Check the propagate method of the fitResult object for evaluation of a vectorized nonparametric model"
    model = _getmodel_axis('nonparametric_vec_normalized',vec=80)

    x = np.linspace(0,10,80)
    fitResult = fit(model,mock_data_fcn(x),x,bootstrap=3)
    
    modeluq = fitResult.propagate(model, x, lb=np.zeros_like(x))

    assert_cis(modeluq)
#================================================================

def test_fit_propagate_semiparametric_bootstrapped(): 
#================================================================
    "Check the propagate method of the fitResult object for evaluation of a semiparametric model"
    model = _getmodel_axis('semiparametric')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x,bootstrap=3)
    
    modeluq = fitResult.propagate(model, x, lb=np.zeros_like(x))

    assert_cis(modeluq)
#================================================================

def test_fit_propagate_semiparametric_vec_bootstrapped(): 
#================================================================
    "Check the propagate method of the fitResult object for evaluation of a vectorized semiparametric model"
    model = _getmodel_axis('semiparametric_vec')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x,bootstrap=3)
    
    modeluq = fitResult.propagate(model, x, lb=np.zeros_like(x))

    assert_cis(modeluq)
#================================================================

def test_fit_propagate_callable_bootstrapped(): 
#================================================================
    "Check the evaluate method of the fitResult object for evaluation of a callable function"
    model = _getmodel_axis('semiparametric_vec')

    x = np.linspace(0,10,200)
    fitResult = fit(model,mock_data_fcn(x),x,bootstrap=3)
    
    fcn = lambda mean1,mean2,std1,std2,amps: gauss2_design_axis(x,mean1,mean2,std1,std2)@amps
    modeluq = fitResult.propagate(fcn, lb=np.zeros_like(x))

    assert_cis(modeluq)
#================================================================

def test_pickle_fitresult():
# ======================================================================
    "Check that the fit results object can be pickled"

    fitResult = fit(mymodel,y)
    store_pickle(fitResult,'pickled_results')
    try:
        pickled_fitResult =read_pickle('pickled_results')
        os.remove("pickled_results.pkl") 
        assert np.allclose(pickled_fitResult.model.real,y.real) and np.allclose(pickled_fitResult.model.imag,y.imag)
    except Exception as exception:
        os.remove("pickled_results.pkl") 
        raise exception
# ======================================================================

def test_pickle_model():
# ======================================================================
    "Check that the model object can be pickled"

    store_pickle(mymodel,'pickled_model')
    try:
        pickled_model =read_pickle('pickled_model')
        fitResult = fit(pickled_model,y)
        os.remove("pickled_model.pkl") 
        assert np.allclose(fitResult.model.real,y.real) and np.allclose(fitResult.model.imag,y.imag)
    except Exception as exception:
        os.remove("pickled_model.pkl") 
        raise exception
# ======================================================================
