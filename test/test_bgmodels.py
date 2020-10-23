import numpy as np
import deerlab as dl

def assert_bgmodel(model,Bref,physical=False):
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
    B5 = model(2.5,par0)
    if physical:
        B6 = model(t,par0,0)
    
    # Assert
    passed = np.zeros(5, dtype=bool)
    passed[0] = all(B1 == B2)
    passed[1] = all(~np.isnan(B1)) and all(~np.isnan(B2)) and all(~np.isnan(B3)) and all(~np.isnan(B4))
    if physical:
        passed[2] = all(B6 == 1)
    else:
        passed[2] = True
    passed[3] = abs(B5 - Bref) < 1e-8
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
    assert not errors, "{}".format(", ".join(errors))
 

def test_bg_hom3d():
    Bref = 0.882785339742350 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(dl.bg_hom3d,Bref,physical=True)

def test_bg_homfractal():
    Bref = 0.882785339742350 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(dl.bg_homfractal,Bref,physical=True)

def test_bg_hom3dex():
    Bref = 0.882896490000000 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(dl.bg_hom3dex,Bref,physical=True)

def test_bg_exp():
    Bref = 0.416862019678508 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(dl.bg_exp,Bref)

def test_bg_strexp():
    Bref = 0.535261428518990 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(dl.bg_strexp,Bref)

def test_bg_prodstrexp():
    Bref = 0.286504796860190 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(dl.bg_prodstrexp,Bref)

def test_bg_sumstrexp():
    Bref = 0.535261428518990 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(dl.bg_sumstrexp,Bref)

def test_bg_poly1():
    Bref = -1.500000000000000 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(dl.bg_poly1,Bref)

def test_bg_poly2():
    Bref = -7.750000000000000 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(dl.bg_poly2,Bref)

def test_bg_poly3():
    Bref = -23.37500000000000 # Reference from DeerLab 0.9.2 on MATLAB
    assert_bgmodel(dl.bg_poly3,Bref)
