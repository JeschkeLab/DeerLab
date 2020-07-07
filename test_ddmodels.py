import numpy as np
from dd_models import *

def assert_ddmodel(model):
    "Check dimensionality, non-negativity, and values at boundaries of a distance distribution model"

    r = np.linspace(0,20,500)
    rnus = np.sqrt(np.linspace(1.5,6^2,500))

    # Extract model information
    info = model()
    par0 = info['Start']
    lower = info['Lower']
    upper = info['Upper']

    # Calculate under different conditions
    P1 = model(r,par0)
    P2 = model(r.T,par0)
    P3 = model(r,lower)
    P4 = model(r,upper)
    P5 = model(rnus,par0)

    # Assert
    passed = np.zeros(5, dtype=bool)
    passed[0] = all(P1 == P2)
    passed[1] = all(P1 >= 0)
    passed[2] = all(P1 >= 0) and all(P2 >= 0)
    passed[3] = all(~np.isnan(P1)) and all(~np.isnan(P2)) and all(~np.isnan(P3)) and all(~np.isnan(P4))
    passed[4] = np.round(np.trapz(P5,rnus),4) == 1
    errors = []
    if not passed[0]:
        errors.append("Dimensionality is not correct")
    if not passed[1]:
        errors.append("Start values violate non-negativity constraint")
    if not passed[2]:
        errors.append("Boundary values violate non-negativity constraint")
    if not passed[3]:
        errors.append("Some conditions have returned NaN values")
    if not passed[4]:
        errors.append("Non-uniform trapezoidal integration failed")

    # assert no error message has been registered, else print messages
    assert not errors, "Errors occured:\n{}".format("\n".join(errors))
 

def test_dd_gauss():
    assert_ddmodel(dd_gauss)

def test_dd_gauss2():
    assert_ddmodel(dd_gauss2)

def test_dd_gauss3():
    assert_ddmodel(dd_gauss3)