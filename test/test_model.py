from collections import namedtuple
from deerlab.whitegaussnoise import whitegaussnoise
import numpy as np 
import pytest
from deerlab.model import Model

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

def test_fit_parametric(): 
#================================================================
    "Check that a parametric model can be correctly fitted"
    model = Model(gauss2)
    model.mean1.set(lb=0, ub=10, par0=2)
    model.mean2.set(lb=0, ub=10, par0=4)
    model.width1.set(lb=0.01, ub=5, par0=0.2)
    model.width2.set(lb=0.01, ub=5, par0=0.2)
    model.amp1.set(lb=0, ub=5, par0=1)
    model.amp2.set(lb=0, ub=5, par0=1)

    fit = model.fit(mock_data)
    
    assert np.allclose(fit.model,mock_data)
#================================================================

def test_fit_semiparametric(): 
#================================================================
    "Check that a semiparametric model can be correctly fitted"
    model = Model(gauss2_design)
    model.mean1.set(lb=0, ub=10, par0=2)
    model.mean2.set(lb=0, ub=10, par0=4)
    model.width1.set(lb=0.01, ub=5, par0=0.2)
    model.width2.set(lb=0.01, ub=5, par0=0.2)
    model.addlinear('amp1',lb=0, ub=5)
    model.addlinear('amp2',lb=0, ub=5)

    fit = model.fit(mock_data)
    
    assert np.allclose(fit.model,mock_data)
#================================================================

def test_fit_nonparametric(): 
#================================================================
    "Check that a semiparametric model can be correctly fitted"
    model = Model(gauss2_identity)
    model.addlinear('gaussian',vec=len(x),lb=np.zeros_like(x))

    fit = model.fit(mock_data,reg=True)
    
    assert np.allclose(fit.model,mock_data,atol=1e-3)
#================================================================