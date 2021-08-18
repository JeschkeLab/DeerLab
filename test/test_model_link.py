import numpy as np
import deerlab as dl 
from deerlab.model import fit,link

# ======================================================================
def test_link_name(): 
    "Check that the parameters are linked to the proper name"
    model = dl.dd_gauss2

    linkedmodel = link(model,mean=[model.mean1,model.mean2])

    assert hasattr(linkedmodel,'mean') and not hasattr(linkedmodel,'mean1') and not hasattr(linkedmodel,'mean2')
# ======================================================================

# ======================================================================
def test_link_Nparam(): 
    "Check that the parameters are linked resulting in proper number of parameters"
    model = dl.dd_gauss2

    linkedmodel = link(model,mean=[model.mean1,model.mean2])

    assert linkedmodel.Nparam == model.Nparam - 1
# ======================================================================

# ======================================================================
def test_link_Nparam_list(): 
    "Check that the combined model has the right number of parameters"
    model = dl.dd_gauss2

    linkedmodel = link(model,mean=[model.mean1,model.mean2])

    assert linkedmodel.Nparam == len(linkedmodel._parameter_list())
# ======================================================================

# ======================================================================
def test_link_multiple(): 
    "Check that multiple links can be established"
    model = dl.dd_gauss3

    linkedmodel = link(model,
            mean1=[model.mean1,model.mean2],
            width2=[model.width2,model.width3])

    assert linkedmodel.Nparam == model.Nparam - 2
# ======================================================================

# ======================================================================
def test_link_call(): 
    "Check that linked parameter models return the correct responses"
    model = dl.dd_gauss2

    linkedmodel = link(model,mean=[model.mean1,model.mean2])

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
                mean=[model.mean1,model.mean2],
                width=[model.width1,model.width2])
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

    linkedmodel = link(model,amp=[model.amp1,model.amp3])

    assert hasattr(linkedmodel,'amp') and not hasattr(linkedmodel,'amp1') and not hasattr(linkedmodel,'amp3')
# ======================================================================

# ======================================================================
def test_link_linear_order(): 
    "Check that the parameters are linked to the proper name"
    model = dl.dd_gauss3

    linkedmodel1 = link(model,amp=[model.amp1,model.amp3])
    linkedmodel2 = link(model,amp=[model.amp3,model.amp1])

    for linkedmodel in [linkedmodel1,linkedmodel2]:
        assert hasattr(linkedmodel,'amp') and not hasattr(linkedmodel,'amp1') and not hasattr(linkedmodel,'amp3')
# ======================================================================

# ======================================================================
def test_link_linear_call(): 
    "Check that linked parameter models return the correct responses"
    model = dl.dd_gauss2

    linkedmodel = link(model,amp=[model.amp1,model.amp2])

    x = np.linspace(0,10,400)
    ref = model(x,4,0.5,2,0.2,0.9,0.9)

    response = linkedmodel(x,4,0.5,2,0.2,0.9)

    assert np.allclose(response,ref)
# ======================================================================
