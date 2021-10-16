import numpy as np
import deerlab as dl 
from deerlab.model import fit,link,Model

# ======================================================================
def test_link_name(): 
    "Check that the parameters are linked to the proper name"
    model = dl.dd_gauss2

    linkedmodel = link(model,mean=['mean1','mean2'])

    assert hasattr(linkedmodel,'mean') and not hasattr(linkedmodel,'mean1') and not hasattr(linkedmodel,'mean2')
# ======================================================================

# ======================================================================
def test_link_Nparam(): 
    "Check that the parameters are linked resulting in proper number of parameters"
    model = dl.dd_gauss2

    linkedmodel = link(model,mean=['mean1','mean2'])

    assert linkedmodel.Nparam == model.Nparam - 1
# ======================================================================

# ======================================================================
def test_link_Nparam_list(): 
    "Check that the combined model has the right number of parameters"
    model = dl.dd_gauss2

    linkedmodel = link(model,mean=['mean1','mean2'])

    assert linkedmodel.Nparam == len(linkedmodel._parameter_list())
# ======================================================================

# ======================================================================
def test_link_multiple(): 
    "Check that multiple links can be established"
    model = dl.dd_gauss3

    linkedmodel = link(model,
            mean1=['mean1','mean2'],
            width2=['width2','width3'])

    assert linkedmodel.Nparam == model.Nparam - 2
# ======================================================================

# ======================================================================
def test_link_call(): 
    "Check that linked parameter models return the correct responses"
    model = dl.dd_gauss2

    linkedmodel = link(model,mean=['mean1','mean2'])

    x = np.linspace(0,10,400)
    ref = model(x,4,0.5,4,0.2,1,1)

    response = linkedmodel(x,4,0.5,0.2,1,1)

    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_link_fit(): 
    "Check that linked parameter models can be properly fitted"
    model = dl.dd_gauss2

    linkedmodel = link(model,
                mean=['mean1','mean2'],
                width=['width1','width2'])
    linkedmodel.mean.par0 = 3
    linkedmodel.width.par0 = 0.5

    x = np.linspace(0,10,400)
    ref = model(x,3,0.5,3,0.5,1,1)

    result = fit(model,ref,x)
    
    assert np.allclose(result.model,ref)
# ======================================================================

# ======================================================================
def test_link_linear_name(): 
    "Check that the parameters are linked to the proper name"
    model = dl.dd_gauss3

    linkedmodel = link(model,amp=['amp1','amp3'])

    assert hasattr(linkedmodel,'amp') and not hasattr(linkedmodel,'amp1') and not hasattr(linkedmodel,'amp3')
# ======================================================================

# ======================================================================
def test_link_linear_order(): 
    "Check that the parameters are linked to the proper name"
    model = dl.dd_gauss3

    linkedmodel1 = link(model,amp=['amp1','amp3'])
    linkedmodel2 = link(model,amp=['amp3','amp1'])

    for linkedmodel in [linkedmodel1,linkedmodel2]:
        assert hasattr(linkedmodel,'amp') and not hasattr(linkedmodel,'amp1') and not hasattr(linkedmodel,'amp3')
# ======================================================================

# ======================================================================
def test_link_linear_call(): 
    "Check that linked parameter models return the correct responses"
    model = dl.dd_gauss2

    linkedmodel = link(model,amp=['amp1','amp2'])

    x = np.linspace(0,10,400)
    ref = model(x,4,0.5,2,0.2,0.9,0.9)

    response = linkedmodel(x,4,0.5,2,0.2,0.9)

    assert np.allclose(response,ref)
# ======================================================================


double_vec = Model(lambda shift, scale: shift + scale*np.eye(80))
double_vec.addlinear('vec1',vec=40,lb=0,par0=0)
double_vec.addlinear('vec2',vec=40,lb=0,par0=0)

# ======================================================================
def test_vec_link_name(): 
    "Check that the parameters are linked to the proper name"
    model = double_vec

    linkedmodel = link(model,vec=['vec1','vec2'])

    assert hasattr(linkedmodel,'vec') and not hasattr(linkedmodel,'vec1') and not hasattr(linkedmodel,'vec2')
# ======================================================================

# ======================================================================
def test_vec_link_Nparam(): 
    "Check that the parameters are linked resulting in proper number of parameters"
    model = double_vec

    linkedmodel = link(model,vec=['vec1','vec2'])

    assert linkedmodel.Nparam == model.Nparam - 40
# ======================================================================

# ======================================================================
def test_vec_link_call(): 
    "Check that linked parameter models return the correct responses"
    model = double_vec

    linkedmodel = link(model,vec=['vec1','vec2'])

    x = np.linspace(0,10,40)
    ref = model(1,2,dl.dd_gauss(x,3,0.2),dl.dd_gauss(x,3,0.2))

    response = linkedmodel(1,2,dl.dd_gauss(x,3,0.2))

    assert np.allclose(response,ref)
# ======================================================================

# ======================================================================
def test_vec_link_fit(): 
    "Check that linked parameter models can be properly fitted"
    model = double_vec

    linkedmodel = link(model,vec=['vec1','vec2'])
    linkedmodel.shift.par0 = 1
    linkedmodel.scale.par0 = 2
    linkedmodel.shift.freeze(1)
    linkedmodel.scale.freeze(2)

    x = np.linspace(0,10,40)
    ref = model(1,2,dl.dd_gauss(x,3,0.6),dl.dd_gauss(x,3,0.6))

    result = fit(linkedmodel,ref,nonlin_tol=1e-3)
    
    assert np.allclose(result.model,ref,atol=1e-2)
# ======================================================================

