import numpy as np
from bg_models import *

def assert_bgmodel(model):
    "Check basic functionality of a background model"

    t = np.linspace(-5,5,500)

    # Extract model information
    info = model()
    par0 = info['Start']
    lower = info['Lower']
    upper = info['Upper']

    # Calculate under different conditions
    B1 = model(t,par0)
    B2 = model(t.T,par0)
    B3 = model(t,lower)
    B4 = model(t,upper)
    B5 = model(t,par0,1)
    B6 = model(t,2*par0,0.5)

    # Assert
    passed = np.zeros(5, dtype=bool)
    passed[0] = all(B1 == B2)
    passed[1] = all(~np.isnan(B1)) and all(~np.isnan(B2)) and all(~np.isnan(B3)) and all(~np.isnan(B4))
    passed[2] = all(B5 == B6)
    errors = []
    if not passed[0]:
        errors.append("Dimensionality is not correct")
    if not passed[1]:
        errors.append("Some conditions have returned NaN values")
    if not passed[2]:
        errors.append("Specifying pathway amplitude leads to incorrect result")

    # assert no error message has been registered, else print messages
    assert not errors, "Errors occured:\n{}".format("\n".join(errors))
 

def test_bg_exp():
    assert_bgmodel(bg_exp)

def test_bg_hom3d():
    assert_bgmodel(bg_hom3d)
