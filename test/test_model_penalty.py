import numpy as np
from deerlab import dd_gauss, whitegaussnoise
from deerlab.model import Penalty
from deerlab.fit import fit
from deerlab.utils import ovl
from copy import deepcopy 
import pytest 

# Fixtures 
# -----------------------------------------------------------------------
x = np.linspace(0,4,50)

available_selection_criteria = ['icc','aic','aicc','bic']

@pytest.fixture(scope='module')
def model():
    return deepcopy(dd_gauss)

@pytest.fixture(scope='module')
def mock_data():
    return dd_gauss(x,2,0.2) + whitegaussnoise(x,0.005,seed=1)

@pytest.fixture(scope='module')
def penalty_fcn(): 
    def _penalty_fcn(mean,std): 
        P = dd_gauss(x,mean,std)
        P = P/np.trapz(P,x)
        return np.sqrt(P*(x - np.trapz(P*x,x))**2*np.mean(np.diff(x)))
    return _penalty_fcn
# -----------------------------------------------------------------------


# ======================================================================
def test_type(penalty_fcn): 
    "Check that the output is the proper object type"
    penaltyobj = Penalty(penalty_fcn,'icc')
    assert isinstance(penaltyobj,Penalty)
# ======================================================================

# ======================================================================
def test_signature(penalty_fcn): 
    "Check that signature of the original function is properly taken"
    penaltyobj = Penalty(penalty_fcn,'icc')
    assert np.all(penaltyobj.signature==['mean','std'])
# ======================================================================

# ======================================================================
@pytest.mark.parametrize('selection',available_selection_criteria)
def test_selection_functionals(penalty_fcn,selection): 
    "Check that penalties with all functional can be constructed"
    penaltyobj = Penalty(penalty_fcn,selection)
    assert penaltyobj.selection==selection
# ======================================================================

# ======================================================================
def test_weight_boundaries(penalty_fcn): 
    "Check that the penalty weight boundaries can be adjusted"
    penaltyobj = Penalty(penalty_fcn,'icc')
    penaltyobj.weight.set(lb=1e-10,ub=1e1)
    assert penaltyobj.weight.lb==1e-10 and penaltyobj.weight.ub==1e1
# ======================================================================

# ======================================================================
def test_weight_freeze(penalty_fcn): 
    "Check that the penalty weight can be frozen to a value"
    penaltyobj = Penalty(penalty_fcn,'icc')
    penaltyobj.weight.freeze(0.5)
    assert penaltyobj.weight.frozen==True and penaltyobj.weight.value==0.5
# ======================================================================

# ======================================================================
@pytest.mark.parametrize('selection',available_selection_criteria)
def test_fit(penalty_fcn, model, mock_data, selection): 
    "Check fitting with a penalty with ICC-selected weight"
    penaltyobj = Penalty(penalty_fcn, selection)
    penaltyobj.weight.set(lb=1e-6,ub=1e1)
    result = fit(model,mock_data,x,penalties=penaltyobj)
    assert ovl(result.model,mock_data)>0.975
# ======================================================================

# ======================================================================
@pytest.mark.parametrize('case',['bounded','frozen'])
def test_fit_with_penalty_weight(penalty_fcn, model, mock_data, case): 
    "Check fitting with a penalty with bounded weight"    
    penaltyobj = Penalty(penalty_fcn,'icc')
    if case == 'bounded':
        penaltyobj.weight.set(lb=1e-10,ub=1e1)
    elif case == 'frozen':
        penaltyobj.weight.freeze(0.00001)
    result = fit(model,mock_data,x,penalties=penaltyobj)
    assert ovl(result.model,mock_data)>0.975
# ======================================================================

# ======================================================================
def test_fit_with_multiple_penalties(penalty_fcn, model, mock_data): 
    "Check fitting with multiple penalties"
    penaltyobj = Penalty(penalty_fcn,'icc')
    penaltyobj2 = deepcopy(penaltyobj)
    penaltyobj.weight.freeze(0.00001)
    penaltyobj2.weight.freeze(0.00001)
    result = fit(model,mock_data,x,penalties=[penaltyobj,penaltyobj2])
    assert ovl(result.model,mock_data)>0.975
# ======================================================================
