from collections import namedtuple
from deerlab.whitegaussnoise import whitegaussnoise
from deerlab.model import Model
from deerlab.fit import fit
import numpy as np 
import pytest
import os 
from deerlab.utils import store_pickle, read_pickle

# Fixtures 
# ----------------------------------------------------------------
x = np.linspace(0,10,80)

def gauss(axis,mean,std): 
    return np.exp(-(axis-mean)**2/std**2/2)
def gauss_fixed_axis(mean,std): 
    return gauss(x,mean,std)

def bigauss(axis,mean1,mean2,std1,std2,amp1,amp2):
    return amp1*gauss(axis,mean1,std1) + amp2*gauss(axis,mean2,std2)
def bigauss_fixed_axis(mean1,mean2,std1,std2,amp1,amp2):
    return bigauss(x,mean1,mean2,std1,std2,amp1,amp2)

def bigauss_design_matrix(axis,mean1,mean2,std1,std2):
    return np.atleast_2d([gauss(axis,mean1,std1), gauss(axis,mean2,std2)]).T
def bigauss_design_matrix_fixed_axis(mean1,mean2,std1,std2):
    return bigauss_design_matrix(x,mean1,mean2,std1,std2)

def identity_matrix(axis):
    return np.eye(len(axis))
def identity_matrix_fixed_axis():
    return identity_matrix(x)

def scaled_identity_matrix(axis,scale):
    return scale*np.eye(len(axis))
def scaled_identity_matrix_fixed_axis(scale):
    return scaled_identity_matrix(x,scale)

@pytest.fixture(scope='module')
def mock_x():
    return x

@pytest.fixture(scope='module')
def mock_data(mock_x):
    return bigauss(mock_x,mean1=3,mean2=4,std1=0.5,std2=0.2,amp1=0.5,amp2=0.6)

model_types = ['parametric','semiparametric','semiparametric_vec',
                'nonparametric','nonparametric_vec','nonparametric_vec_normalized']

def _generate_model(type,fixed_axis=False):
    if type=='parametric':
        if fixed_axis:
            model = Model(bigauss_fixed_axis)
        else:
            model = Model(bigauss,constants='axis')
        model.mean1.set(lb=0, ub=10, par0=2)
        model.mean2.set(lb=0, ub=10, par0=4)
        model.std1.set(lb=0.1, ub=5, par0=0.2)
        model.std2.set(lb=0.1, ub=5, par0=0.2)
        model.amp1.set(lb=0, ub=5, par0=1)
        model.amp2.set(lb=0, ub=5, par0=1)
    elif type=='semiparametric': 
        if fixed_axis:
            model = Model(bigauss_design_matrix_fixed_axis)
        else:
            model = Model(bigauss_design_matrix,constants='axis')
        model.mean1.set(lb=0, ub=10, par0=2)
        model.mean2.set(lb=0, ub=10, par0=4)
        model.std1.set(lb=0.1, ub=5, par0=0.2)
        model.std2.set(lb=0.1, ub=5, par0=0.2)
        model.addlinear('amp1',lb=0, ub=5)
        model.addlinear('amp2',lb=0, ub=5)
    elif type=='semiparametric_vec': 
        if fixed_axis:
            model = Model(bigauss_design_matrix_fixed_axis)
        else:
            model = Model(bigauss_design_matrix,constants='axis')
        model.mean1.set(lb=0, ub=10, par0=2)
        model.mean2.set(lb=0, ub=10, par0=4)
        model.std1.set(lb=0.01, ub=5, par0=0.2)
        model.std2.set(lb=0.01, ub=5, par0=0.2)
        model.addlinear('amps',lb=0, ub=5, vec=2)
    elif type=='nonparametric':
        if fixed_axis:
            model = Model(bigauss_design_matrix_fixed_axis(3,4,0.5,0.2))
        else:
            model = Model(lambda axis: bigauss_design_matrix(axis,3,4,0.5,0.2),constants='axis')        
        model.addlinear('amp1',lb=0)
        model.addlinear('amp2',lb=0)
    elif type=='nonparametric_vec': 
        if fixed_axis:
            model = Model(identity_matrix_fixed_axis)
        else:
            model = Model(identity_matrix,constants='axis')       
        model.addlinear('dist',lb=0,vec=len(x))
    elif type=='nonparametric_vec_normalized':
        if fixed_axis:
            model = Model(identity_matrix_fixed_axis)
        else:
            model = Model(identity_matrix,constants='axis')    
        model.addlinear('dist',lb=0,vec=len(x),normalization=lambda dist:dist/np.sum(dist))        
    return model
# ----------------------------------------------------------------

def test_construction_length():
# ================================================================
    "Check that the model is contructed correctly with the appropiate number of parameters"
    model = Model(gauss_fixed_axis)

    assert model.Nparam==2
# ================================================================

def test_construction_names():
# ================================================================
    "Check that the model is contructed correctly with the correct parameter names"
    model = Model(gauss_fixed_axis)

    assert 'mean' in model.__dict__ and 'std' in model.__dict__
# ================================================================

def test_parameters_set():
# ================================================================
    "Check that attributes of the parameters are editable"
    model = Model(gauss_fixed_axis)
    model.mean.set(lb=0,ub=10)

    assert getattr(model.mean,'lb')==0 and getattr(model.mean,'ub')==10 
# ================================================================

def test_parameters_setas():
# ================================================================
    "Check that attributes can be set to that of another parameter"
    modelref = Model(gauss_fixed_axis)
    modelref.mean.set(lb=0,ub=10)

    model = Model(gauss_fixed_axis)
    model.mean.setas(modelref.mean)

    assert getattr(model.mean,'lb')==getattr(modelref.mean,'lb') and getattr(model.mean,'ub')==getattr(modelref.mean,'ub') 
# ================================================================

def test_call_keywords():
# ================================================================
    "Check that calling the model with parameter returns the correct response"
    model = Model(gauss_fixed_axis)
    response  = model(mean=3,std=0.5)
    reference = gauss_fixed_axis(mean=3,std=0.5)

    assert np.allclose(response,reference)
# ================================================================

def test_call_positional():
# ================================================================
    "Check that calling the model with parameter returns the correct response"
    model = Model(gauss_fixed_axis)
    response  = model(3,0.5)
    reference = gauss_fixed_axis(3,0.5)
    
    assert np.allclose(response,reference)
# ================================================================

def test_call_mixed():
# ================================================================
    "Check that calling the model with parameter returns the correct response"
    model = Model(gauss_fixed_axis)
    response  = model(3,std=0.5)
    reference = gauss_fixed_axis(3,0.5)

    assert np.allclose(response,reference)
# ================================================================


def test_addlinear_length():
# ================================================================
    "Check that the model is contructed correctly with the appropiate number of parameters"
    model = Model(bigauss_design_matrix_fixed_axis)
    model.addlinear('amp1',lb=0)
    model.addlinear('amp2',lb=0)

    assert model.Nparam==6
# ================================================================

def test_addlinear_names():
# ================================================================
    "Check that linear parameters can be properly added"
    model = Model(bigauss_design_matrix_fixed_axis)
    model.addlinear('amp1',lb=0)
    model.addlinear('amp2',lb=0)

    assert 'amp1' in model.__dict__ and 'amp2' in model.__dict__
# ================================================================

def test_addlinear_set():
# ================================================================
    "Check that attributes of the linear parameters are editable"
    model = Model(bigauss_design_matrix_fixed_axis)
    model.addlinear('amp1')
    model.amp1.set(lb=0,ub=10)

    assert getattr(model.amp1,'lb')==0 and getattr(model.amp1,'ub')==10 
# ================================================================

def test_addlinear_normalization():
# ================================================================
    "Check that linear parameters can be properly added with a normalization function"
    model = Model(bigauss_design_matrix_fixed_axis)
    model.addlinear('amp1',lb=0, normalization=lambda x: 2*x)
    model.addlinear('amp2',lb=0, normalization=lambda x: 3*x)

    assert 'amp1' in model.__dict__ and 'amp2' in model.__dict__
# ================================================================

def test_addlinear_call_keywords():
# ================================================================
    "Check that calling the model with parameters returns the correct response"
    model = Model(bigauss_design_matrix_fixed_axis)
    model.addlinear('amp1')
    model.addlinear('amp2')
    response = model(mean1=3,mean2=4,std1=0.2,std2=0.3,amp1=0.5,amp2=0.4)
    reference = bigauss_fixed_axis(mean1=3,mean2=4,std1=0.2,std2=0.3,amp1=0.5,amp2=0.4)

    assert np.allclose(response,reference)
# ================================================================


def test_addlinear_call_positional():
# ================================================================
    "Check that calling the model with parameters returns the correct response"
    model = Model(bigauss_design_matrix_fixed_axis)
    model.addlinear('amp1')
    model.addlinear('amp2')
    response = model(3,4,0.2,0.3,0.5,0.4)
    reference = bigauss_fixed_axis(3,4,0.2,0.3,0.5,0.4)

    assert np.allclose(response,reference)
# ================================================================

def test_addlinear_call_mixed():
# ================================================================
    "Check that calling the model with parameters returns the correct response"
    model = Model(bigauss_design_matrix_fixed_axis)
    model.addlinear('amp1')
    model.addlinear('amp2')
    response = model(3,4,0.2,std2=0.3,amp1=0.5,amp2=0.4)
    reference = bigauss_fixed_axis(3,4,0.2,std2=0.3,amp1=0.5,amp2=0.4)

    assert np.allclose(response,reference)
# ================================================================

# ================================================================
def test_addlinear_vector_length():
    "Check that linear parameters can be defined as vectors"
    model = Model(identity_matrix_fixed_axis)
    model.addlinear('gaussian', vec=100)

    assert model.Nparam==100
# ================================================================

# ================================================================
def test_addlinear_vector_names():
    "Check that linear parameters can be defined as vectors"
    model = Model(identity_matrix_fixed_axis)
    model.addlinear('gaussian', vec=100)

    assert 'gaussian' in model.__dict__ 
# ================================================================

def test_addlinear_vector_set():
# ================================================================
    "Check that attributes of the vector linear parameters are editable"
    model = Model(identity_matrix_fixed_axis)
    model.addlinear('gaussian', vec=100)
    model.gaussian.set(lb=np.zeros(100))

    assert np.allclose(getattr(model.gaussian,'lb'),np.zeros(100)) 
# ================================================================

def test_addlinear_vector_call_keywords():
# ================================================================
    "Check that calling the model with scalar and vector parameters returns the correct response"
    model = Model(identity_matrix_fixed_axis)
    model.addlinear('gaussian', vec=len(x))
    reference = bigauss_fixed_axis(mean1=3,mean2=4,std1=0.2,std2=0.3,amp1=0.5,amp2=0.4)
    response = model(gaussian=reference)

    assert np.allclose(response,reference)
# ================================================================

def test_addlinear_vector_call_positional():
# ================================================================
    "Check that calling the model with scalar and vector parameters returns the correct response"
    model = Model(identity_matrix_fixed_axis)
    model.addlinear('gaussian', vec=len(x))
    reference = bigauss_fixed_axis(3,4,0.2,0.3,0.5,0.4)
    response = model(reference)

    assert np.allclose(response,reference)
# ================================================================

def test_mixed_vector_length():
# ================================================================
    "Check the definition of scalar nonlinear parameters and vector linear parameters"
    model = Model(scaled_identity_matrix_fixed_axis)
    model.addlinear('gaussian', vec=100)

    assert model.Nparam==101
# ================================================================

def test_mixed_vector_names():
# ================================================================
    "Check the definition of scalar nonlinear parameters and vector linear parameters"
    model = Model(scaled_identity_matrix_fixed_axis)
    model.addlinear('gaussian', vec=100)

    assert 'gaussian' in model.__dict__  and 'scale' in model.__dict__ 
# ================================================================

def test_mixed_vector_call_keywords():
# ================================================================
    "Check that calling the model with scalar and vector parameters returns the correct response"
    model = Model(scaled_identity_matrix_fixed_axis)
    model.addlinear('gaussian', vec=len(x))
    reference = 5*bigauss_fixed_axis(mean1=3,mean2=4,std1=0.2,std2=0.3,amp1=0.5,amp2=0.4)
    response = model(scale=5,gaussian=reference/5)

    assert np.allclose(response,reference)
# ================================================================

def test_mixed_vector_call_positional():
# ================================================================
    "Check that calling the model with scalar and vector parameters returns the correct response"
    model = Model(scaled_identity_matrix_fixed_axis)
    model.addlinear('gaussian', vec=len(x))
    reference = 5*bigauss_fixed_axis(3,4,0.2,0.3,0.5,0.4)
    response = model(5,reference/5)

    assert np.allclose(response,reference)
# ================================================================

def test_addnonlinear_length():
# ================================================================
    "Check that the model is contructed correctly with the appropiate number of parameters"
    model = Model(gauss_fixed_axis)
    model.addnonlinear('trivial1',lb=0)
    model.addnonlinear('trivial2',lb=0)

    assert model.Nparam==4
# ================================================================

def test_addnonlinear_names():
# ================================================================
    "Check that linear parameters can be properly added"
    model = Model(gauss_fixed_axis)
    model.addnonlinear('trivial1',lb=0)
    model.addnonlinear('trivial2',lb=0)

    assert hasattr(model,'trivial1') and hasattr(model,'trivial2')
# ================================================================

def test_addnonlinear_set():
# ================================================================
    "Check that attributes of the linear parameters are editable"
    model = Model(gauss_fixed_axis)
    model.addnonlinear('trivial1')
    model.trivial1.set(lb=0,ub=10)

    assert getattr(model.trivial1,'lb')==0 and getattr(model.trivial1,'ub')==10 
# ================================================================

def test_addnonlinear_call_keywords():
# ================================================================
    "Check that calling the model with parameters returns the correct response"
    model = Model(gauss_fixed_axis)
    model.addnonlinear('trivial1')
    model.addnonlinear('trivial2')
    response = model(mean=3,std=0.2,trivial1=1,trivial2=1)
    reference = gauss_fixed_axis(mean=3,std=0.2)

    assert np.allclose(response,reference)
# ================================================================

def test_addnonlinear_call_positional():
# ================================================================
    "Check that calling the model with parameters returns the correct response"
    model = Model(gauss_fixed_axis)
    model.addnonlinear('trivial1')
    model.addnonlinear('trivial2')
    response = model(3,0.2,1,1)
    reference = gauss_fixed_axis(3,0.2)

    assert np.allclose(response,reference)
# ================================================================

# ======================================================================
def test_preserve_original(mock_data): 
    "Check that the original model is not changed by the function"
    model = Model(gauss_fixed_axis)
    model.mean.par0 = 3
    model.std.par0 = 0.2
    
    _ = fit(model,mock_data)
    assert model._parameter_list() == ['mean','std']
# ======================================================================

def test_freeze():
# ================================================================
    "Check that a model parameter can be frozen to a fixed value"
    model = Model(gauss_fixed_axis)
    model.mean.freeze(3)

    assert model.mean.value==3 and model.mean.frozen==True
# ================================================================

def test_freeze_outofbounds():
# ================================================================
    "Check that a parameter cannot be frozen outside of the bounds"
    model = Model(gauss_fixed_axis)
    model.mean.set(lb=0,ub=10)
    with pytest.raises(ValueError):
        model.mean.freeze(-10)
# ================================================================

def test_freeze_vec_outofbounds():
# ================================================================
    "Check that a parameter cannot be frozen outside of the bounds"
    model = Model(bigauss_design_matrix_fixed_axis)
    model.addlinear('gaussian', vec=100, lb=0)
    with pytest.raises(ValueError):
        model.gaussian.freeze(np.full(100,-10))
# ================================================================

def test_unfreeze():
# ================================================================
    "Check that a model parameter can be frozen and then reversed"
    model = Model(gauss)
    model.mean.freeze(3)
    model.mean.unfreeze()

    assert model.mean.frozen==False and model.mean.value==None
# ================================================================

def test_model_with_constant_positional(): 
# ================================================================
    "Check that a model with axis can be defined and called"
    model = Model(gauss,constants='axis')

    x = np.linspace(0,10,300)
    reference = gauss(x,3,0.5)
    response = model(x,3,0.5)
        
    assert np.allclose(reference,response)
# ================================================================

def test_model_with_constant_keywords(): 
# ================================================================
    "Check that a model with axis can be defined and called via keywords"
    model = Model(gauss,constants='axis')

    x = np.linspace(0,10,300)
    reference = gauss(x,3,0.5)
    response = model(axis=x, mean=3, std=0.5)
        
    assert np.allclose(reference,response)
# ================================================================

def test_model_with_constant_mixed(): 
# ================================================================
    "Check that a model with axis can be defined and called via keywords"
    model = Model(gauss,constants='axis')

    x = np.linspace(0,10,300)
    reference = gauss(x,3,0.5)
    response = model(x, mean=3, std=0.5)
        
    assert np.allclose(reference,response)
# ================================================================

@pytest.mark.parametrize('model_type', model_types)
def test_fit_model(mock_data,model_type): 
# ================================================================
    "Check that a parametric model can be correctly fitted"
    model = _generate_model(model_type, fixed_axis=True)

    result = fit(model, mock_data)
    
    assert np.allclose(result.model, mock_data)
# ================================================================

@pytest.mark.parametrize('model_type, frozen_params', [
    ('parametric', [['mean1', 3],['mean2', 4]]),    
    ('semiparametric', [['amp1', 0.5],['amp2', 0.6]]),
    ('semiparametric', [['mean1', 3],['mean2', 4],['amp1', 0.5], ['amp2', 0.6]]),
])
def test_fit_model_with_frozen_parameters(mock_data, model_type, frozen_params): 
# ================================================================
    "Check that a parametric model can be correctly fitted"
    model = _generate_model(model_type, fixed_axis=True)

    for param,value in frozen_params:
        getattr(model,param).freeze(value)
    result = fit(model,mock_data)
    
    assert np.allclose(result.model,mock_data)
# ================================================================

@pytest.mark.parametrize('method', ['bootstrap','moment'])
@pytest.mark.parametrize('model_type', model_types)
def test_fit_model_confidence_intervals(mock_data,model_type,method): 
# ================================================================
    "Check the default confidence intervals of the fitted parameters"
    model = _generate_model(model_type, fixed_axis=True)

    if method=='bootstrap':
        results = fit(model,mock_data + whitegaussnoise(x,0.01,seed=1), bootstrap=3)
    else: 
        results = fit(model,mock_data + whitegaussnoise(x,0.01,seed=1))
    
    for param in model._parameter_list(): 
        parfit = getattr(results,param)
        ci_lower,ci_upper = getattr(results,f'{param}Uncert').ci(95).T
        if getattr(results,f'{param}Uncert').type=='bootstrap':
            assert np.allclose(np.round(parfit,10),np.round(getattr(results,f'{param}Uncert').median,10), atol=1e-3)
        assert np.all(np.round(parfit,12)<=np.round(ci_upper,12)) and np.all(np.round(parfit,12)>=np.round(ci_lower,12)) 
# ================================================================

# ================================================================
@pytest.mark.parametrize('model_type', model_types)
def test_fit_model_with_constant(mock_data,mock_x,model_type): 
    "Check that a model can be correctly fitted while specifying its constant"

    model = _generate_model(model_type)
    results = fit(model,mock_data,mock_x,noiselvl=1e-10)
    
    assert np.allclose(results.model,mock_data)
# ================================================================

# ================================================================
@pytest.mark.parametrize('model_type', model_types)
def test_fit_evaluate_model(mock_data,mock_x,model_type): 
    "Check that models can be evaluated from fit results"

    model = _generate_model(model_type, fixed_axis=False)
    results = fit(model,mock_data,mock_x,noiselvl=1e-10)
    response = results.evaluate(model,mock_x)
    if hasattr(results,'scale'):
        response *= results.scale
    
    assert np.allclose(response,mock_data)
# ================================================================

# ================================================================
@pytest.mark.parametrize('method', ['bootstrap','moment'])
@pytest.mark.parametrize('model_type', model_types)
def test_fit_propagate_through_model(mock_data,mock_x,model_type,method): 
    "Check that the uncertainty of fit results can be propagated through models"

    model = _generate_model(model_type, fixed_axis=False)
    
    if method=='bootstrap':
        results = fit(model,mock_data,mock_x, bootstrap=3)
    else: 
        results = fit(model,mock_data,mock_x)

    response_ci = results.propagate(model, mock_x, lb=np.zeros_like(mock_x)).ci(95)
    ci_lower = response_ci[:,0]
    ci_upper = response_ci[:,1]

    assert np.less_equal(ci_lower,ci_upper).all()
# ================================================================


@pytest.mark.parametrize('method', ['bootstrap','moment'])
@pytest.mark.parametrize('model_type', ['parametric'])
def test_fit_propagate_through_function(mock_data,mock_x,model_type,method):
    
    model = _generate_model(model_type, fixed_axis=False)
    
    if method=='bootstrap':
        results = fit(model,mock_data,mock_x, bootstrap=3)
    else: 
        results = fit(model,mock_data,mock_x)
    
    if model_type=='parametric':
        fun = lambda mean1,mean2,std1,std2,amp1,amp2: bigauss(mock_x,mean1,mean2,std1,std2,amp1,amp2)
        response_ci = results.propagate(fun, lb=np.zeros_like(mock_x)).ci(95)
    
    ci_lower = response_ci[:,0]
    ci_upper = response_ci[:,1]

    assert np.less_equal(ci_lower,ci_upper).all()

    

# ================================================================

def gauss_multiaxis(axis1,axis2,mean,std): 
    return np.exp(-(axis1-mean)**2/std**2/2)

def gauss2_multiaxis(axis1,axis2,mean1,mean2,std1,std2,amp1,amp2):
    return amp1*gauss(axis1,mean1,std1) + amp2*gauss(axis2,mean2,std2)

def test_model_with_multiple_constants(): 
# ================================================================
    "Check that a model with axis can be defined and called"
    model = Model(gauss_multiaxis,constants=['axis1','axis2'])

    x1 = np.linspace(0,5,300)
    x2 = np.linspace(5,10,300)
    reference = gauss_multiaxis(x1,x2,3,0.5)
    response = model(x1,x2,3,0.5)
        
    assert np.allclose(reference,response)
# ================================================================


def test_model_with_multiple_constants_fit(mock_data): 
# ================================================================
    "Check that a model with axis can be defined and called"
    model = Model(gauss2_multiaxis,constants=['axis1','axis2'])
    model.mean1.set(lb=0, ub=10, par0=2)
    model.mean2.set(lb=0, ub=10, par0=4)
    model.std1.set(lb=0.01, ub=5, par0=0.2)
    model.std2.set(lb=0.01, ub=5, par0=0.2)
    model.amp1.set(lb=0, ub=5, par0=1)
    model.amp2.set(lb=0, ub=5, par0=1)
    
    x1 = np.linspace(0,10,80)
    x2 = np.linspace(0,10,80)
    result = fit(model,mock_data,x1,x2)
    
    assert np.allclose(result.model,mock_data)        
# ================================================================

def test_model_constant_same_values_keywords(): 
# ================================================================
    "Check that a model with axis can be defined and called with same values as parameters"
    model = Model(gauss,constants='axis')

    reference = gauss(axis=3,mean=3,std=0.5)
    response = model(axis=3,mean=3,std=0.5)
        
    assert np.allclose(reference,response)
# ================================================================

def test_model_constant_same_values_positional(): 
# ================================================================
    "Check that a model with axis can be defined and called with same values as parameters"
    model = Model(gauss,constants='axis')

    reference = gauss(3,3,0.5)
    response = model(3,3,0.5)
        
    assert np.allclose(reference,response)
# ================================================================

# Fixtures 
# ---------------------------------------------------------------
@pytest.fixture(scope='module')
def complex_model():
    def _model(phase, center, std):
        y = gauss_fixed_axis(center, std)
        y = y*np.exp(-1j*phase)
        return y
    complex_model = Model(_model)
    complex_model.phase.set(par0=2*np.pi/5, lb=-np.pi, ub=np.pi)
    complex_model.center.set(par0=4, lb=1, ub=6)
    complex_model.std.set(par0=0.2, lb=0.05, ub=5)
    return complex_model 

@pytest.fixture(scope='module')
def complex_data(complex_model):
    return complex_model(phase=np.pi/5, center=3, std=0.5)     
# ---------------------------------------------------------------


def test_complex_model_complex_data(complex_model, complex_data):
# ======================================================================
    "Check the fit of a complex-valued model to complex-valued data"
    result = fit(complex_model, complex_data)

    assert np.allclose(result.model.real, complex_data.real) and np.allclose(result.model.imag, complex_data.imag)
# ======================================================================

def test_complex_model_real_data(complex_model, complex_data):
# ======================================================================
    "Check the fit of a complex-valued model to real-valued data"
    result = fit(complex_model, complex_data.real)

    assert np.allclose(result.model.real, complex_data.real) and np.allclose(result.model.imag, np.zeros_like(complex_data.real))
# ======================================================================

def test_real_model_complex_data(complex_model, complex_data):
# ======================================================================
    "Check the fit of a real-valued model to complex-valued data"
    complex_model.phase.freeze(0)
    result = fit(complex_model, complex_data.real + 1j*np.zeros_like(complex_data.imag))
    complex_model.phase.unfreeze()
    
    assert np.allclose(result.model.real, complex_data.real) and np.allclose(result.model.imag, np.zeros_like(complex_data.real))
# ======================================================================

def test_pickle_result(complex_model,complex_data):
# ======================================================================
    "Check that the fit results object can be pickled"

    result = fit(complex_model,complex_data)
    store_pickle(result,'pickled_results')
    try:
        pickled_result =read_pickle('pickled_results')
        os.remove("pickled_results.pkl") 
        assert np.allclose(pickled_result.model.real,complex_data.real) and np.allclose(pickled_result.model.imag,complex_data.imag)
    except Exception as exception:
        os.remove("pickled_results.pkl") 
        raise exception
# ======================================================================

def test_pickle_model(complex_model,complex_data):
# ======================================================================
    "Check that the model object can be pickled"

    store_pickle(complex_model,'pickled_model')
    try:
        pickled_model =read_pickle('pickled_model')
        result = fit(pickled_model,complex_data)
        os.remove("pickled_model.pkl") 
        assert np.allclose(result.model.real,complex_data.real) and np.allclose(result.model.imag,complex_data.imag)
    except Exception as exception:
        os.remove("pickled_model.pkl") 
        raise exception
# ======================================================================

# Fixtures
# ----------------------------------------------------------------------
@pytest.fixture(scope='module')
def mask(mock_x):
    # Construct mask 
    mask = np.ones_like(mock_x).astype(bool)   
    mask[(mock_x>3.8) & (mock_x<4.2)] = False 
    return mask 

@pytest.fixture(scope='module')
def mock_data_corrupted(mask,mock_data):
    corrupted_data = mock_data.copy()
    # Distort data outside of mask
    corrupted_data[~mask] = 0 
    return corrupted_data
# ----------------------------------------------------------------------

def test_fit_masking(mock_data_corrupted, mock_x, mock_data, mask): 
# ================================================================
    "Check that masking works"
    model = _generate_model('parametric', fixed_axis=False)

    results = fit(model, mock_data_corrupted,mock_x)
    results_masked = fit(model, mock_data_corrupted,mock_x,mask=mask)

    assert np.allclose(results_masked.model, mock_data) and not np.allclose(results.model, mock_data) 
# ================================================================

def test_masking_noiselvl(mock_data_corrupted, mock_x, mask):
# ======================================================================
    "Check that masking leads to correct noise level estimates"
    model = _generate_model('parametric', fixed_axis=False)

    noiselvl = 0.1
    noisy_data = mock_data_corrupted + whitegaussnoise(mock_x,noiselvl,seed=1)

    results = fit(model,noisy_data,mock_x)
    results_masked = fit(model,noisy_data,mock_x,mask=mask)

    assert results_masked.noiselvl!=results.noiselvl and abs(results_masked.noiselvl/noiselvl-1)<0.1
# ======================================================================

def test_masking_chi2red(mock_data_corrupted, mock_x, mask):
# ======================================================================
    "Check that masking is accounted for by the goodness-of-fit"
    model = _generate_model('parametric', fixed_axis=False)

    noiselvl = 0.02
    noisy_data = mock_data_corrupted + whitegaussnoise(x,noiselvl,seed=1)

    results = fit(model,noisy_data,mock_x, noiselvl=noiselvl)
    results_masked = fit(model,noisy_data,mock_x,mask=mask,noiselvl=noiselvl)

    assert results_masked.stats['chi2red']<results.stats['chi2red'] and abs(results_masked.stats['chi2red']-1)<0.2
# ======================================================================