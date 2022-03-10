import numpy as np
from deerlab import dd_gauss, whitegaussnoise
from deerlab.model import Penalty,fit
from deerlab.utils import ovl
from copy import deepcopy 

x = np.linspace(0,4,50)
mock_data = dd_gauss(x,2,0.2) + whitegaussnoise(x,0.005,seed=1)
def penalty_fcn(mean,std): 
    P = dd_gauss(x,mean,std)
    P = P/np.trapz(P,x)
    return np.sqrt(P*(x - np.trapz(P*x,x))**2*np.mean(np.diff(x)))

# ======================================================================
def test_type(): 
    "Check that the output is the proper object type"
    
    penaltyobj = Penalty(penalty_fcn,'icc')

    assert isinstance(penaltyobj,Penalty)
# ======================================================================

# ======================================================================
def test_signature(): 
    "Check that signature of the original function is properly taken"
    
    penaltyobj = Penalty(penalty_fcn,'icc')

    assert np.all(penaltyobj.signature==['mean','std'])
# ======================================================================

# ======================================================================
def test_functionals_icc(): 
    "Check that penalties with all functional can be constructed"
    
    penaltyobj = Penalty(penalty_fcn,'icc')

    assert penaltyobj.selection=='icc'
# ======================================================================


# ======================================================================
def test_functionals_aic(): 
    "Check that penalties with all functional can be constructed"
    
    penaltyobj = Penalty(penalty_fcn,'aic')

    assert penaltyobj.selection=='aic'
# ======================================================================


# ======================================================================
def test_functionals_bic(): 
    "Check that penalties with all functional can be constructed"
    
    penaltyobj = Penalty(penalty_fcn,'bic')

    assert penaltyobj.selection=='bic'
# ======================================================================


# ======================================================================
def test_functionals_aicc(): 
    "Check that penalties with all functional can be constructed"
    
    penaltyobj = Penalty(penalty_fcn,'aicc')

    assert penaltyobj.selection=='aicc'
# ======================================================================

# ======================================================================
def test_weight_boundaries(): 
    "Check that the penalty weight boundaries can be adjusted"

    penaltyobj = Penalty(penalty_fcn,'icc')
    penaltyobj.weight.set(lb=1e-10,ub=1e1)

    assert penaltyobj.weight.lb==1e-10 and penaltyobj.weight.ub==1e1
# ======================================================================

# ======================================================================
def test_weight_freeze(): 
    "Check that the penalty weight can be frozen to a value"

    penaltyobj = Penalty(penalty_fcn,'icc')
    penaltyobj.weight.freeze(0.5)

    assert penaltyobj.weight.frozen==True and penaltyobj.weight.value==0.5
# ======================================================================

# ======================================================================
def test_fit_icc(): 
    "Check fitting with a penalty with ICC-selected weight"

    model = deepcopy(dd_gauss)
    penaltyobj = Penalty(penalty_fcn,'icc')
    penaltyobj.weight.set(lb=1e-6,ub=1e1)

    result = fit(model,mock_data,x,penalties=penaltyobj)

    assert ovl(result.model,mock_data)>0.975
# ======================================================================

# ======================================================================
def test_fit_aic(): 
    "Check fitting with a penalty with AIC-selected weight"

    model = deepcopy(dd_gauss)
    penaltyobj = Penalty(penalty_fcn,'aic')
    penaltyobj.weight.set(lb=1e-6,ub=1e1)

    result = fit(model,mock_data,x,penalties=penaltyobj)

    assert ovl(result.model,mock_data)>0.975
# ======================================================================

# ======================================================================
def test_fit_bic(): 
    "Check fitting with a penalty with AIC-selected weight"
 
    model = deepcopy(dd_gauss)
    penaltyobj = Penalty(penalty_fcn,'bic')
    penaltyobj.weight.set(lb=1e-6,ub=1e1)

    result = fit(model,mock_data,x,penalties=penaltyobj)

    assert ovl(result.model,mock_data)>0.975
# ======================================================================

# ======================================================================
def test_fit_aicc(): 
    "Check fitting with a penalty with AICc-selected weight"

    model = deepcopy(dd_gauss)
    penaltyobj = Penalty(penalty_fcn,'aicc')
    penaltyobj.weight.set(lb=1e-6,ub=1e1)

    result = fit(model,mock_data,x,penalties=penaltyobj)

    assert ovl(result.model,mock_data)>0.975
# ======================================================================

# ======================================================================
def test_fit_weight_bounded(): 
    "Check fitting with a penalty with bounded weight"
    
    model = deepcopy(dd_gauss)
    penaltyobj = Penalty(penalty_fcn,'aicc')
    penaltyobj.weight.set(lb=1e-10,ub=1e1)

    result = fit(model,mock_data,x,penalties=penaltyobj)

    assert ovl(result.model,mock_data)>0.975
# ======================================================================

# ======================================================================
def test_fit_weight_unbounded(): 
    "Check fitting with a penalty with unbounded weight"
    
    model = deepcopy(dd_gauss)
    penaltyobj = Penalty(penalty_fcn,'icc')

    result = fit(model,mock_data,x,penalties=penaltyobj)

    assert not ovl(result.model,mock_data)>0.975
# ======================================================================

# ======================================================================
def test_fit_weight_frozen(): 
    "Check fitting with a penalty with frozen weight"
    
    model = deepcopy(dd_gauss)
    penaltyobj = Penalty(penalty_fcn,'icc')
    penaltyobj.weight.freeze(0.00001)

    result = fit(model,mock_data,x,penalties=penaltyobj)

    assert ovl(result.model,mock_data)>0.975
# ======================================================================

# ======================================================================
def test_fit_multiple_penalties(): 
    "Check fitting with multiple penalties"
    
    model = deepcopy(dd_gauss)
    penaltyobj = Penalty(penalty_fcn,'icc')
    penaltyobj2 = deepcopy(penaltyobj)
    penaltyobj.weight.freeze(0.00001)
    penaltyobj2.weight.freeze(0.00001)

    result = fit(model,mock_data,x,penalties=[penaltyobj,penaltyobj2])

    assert ovl(result.model,mock_data)>0.975
# ======================================================================
