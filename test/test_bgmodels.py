import numpy as np
from deerlab.bg_models import *

def assert_bgmodel(model,Bref):
    "Check the correct behaviour of the core functionality of a background model"

    t = np.linspace(-5,5,500)

    # Extract model information
    info = model()
    par0 = info['Start']
    lower = info['Lower']
    upper = info['Upper']
    paramnames = info['Parameters']
    units = info['Units']

    # Calculate under different conditions
    B1 = model(t,par0)
    B2 = model(t.T,par0)
    B3 = model(t,lower)
    B4 = model(t,upper)
    B5 = model(t,par0,0)
    B6 = model(2.5,par0)
    print(B6)
    # Assert
    passed = np.zeros(5, dtype=bool)
    passed[0] = all(B1 == B2)
    passed[1] = all(~np.isnan(B1)) and all(~np.isnan(B2)) and all(~np.isnan(B3)) and all(~np.isnan(B4))
    passed[2] = all(B5 == 1)
    passed[3] = abs(B6 - Bref) < 1e-8
    passed[4] = len(paramnames) == len(par0) and len(units) == len(par0)
    errors = []
    if not passed[0]:
        errors.append("Dimensionality is not correct")
    if not passed[1]:
        errors.append("Some conditions have returned NaN values")
    if not passed[2]:
        errors.append("Specifying pathway amplitude leads to incorrect result")
    if not passed[3]:
        errors.append("Absolute values of returned by model do no match the reference.")
    if not passed[4]:
        errors.append("The number of parameter names and units are not equal.")
    # assert no error message has been registered, else print messages
    assert not errors, "Errors occured:\n{}".format("\n".join(errors))
 

def test_bg_exp():
    Bref = 0.416862019678508 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(bg_exp,Bref)

def test_bg_hom3d():
    Bref = 0.882785339742350 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(bg_hom3d,Bref)

def test_bg_hom3dex():
    Bref = 0.882896490000000 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(bg_hom3dex,Bref)

def test_bg_strexp():
    Bref = 0.535261428518990 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(bg_strexp,Bref)

def test_bg_prodstrexp():
    Bref = 0.286504796860190 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(bg_prodstrexp,Bref)

def test_bg_sumstrexp():
    Bref = 0.535261428518990 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(bg_sumstrexp,Bref)

def test_bg_poly1():
    Bref = -1.500000000000000 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(bg_poly1,Bref)

def test_bg_poly2():
    Bref = -7.750000000000000 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(bg_poly2,Bref)

def test_bg_poly3():
    Bref = -23.37500000000000 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(bg_poly3,Bref)