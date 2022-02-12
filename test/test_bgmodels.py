import numpy as np
import deerlab as dl


 #---------------------------------------------------------------------------------------
def assert_bgmodel(model,Bref):
    "Check the correct behaviour of the core functionality of a background model"

    t = np.linspace(-5,5,500)

    # Extract model information
    meta = model.getmetadata()
    par0 = meta['par0']
    lower = meta['lb']
    upper = meta['ub']
    paramnames = meta['names']
    units = meta['units']

    # Calculate under different conditions
    B1 = model(t,*par0)
    B2 = model(t.T,*par0)
    B3 = model(t,*lower)
    B4 = model(t,*upper)
    B5 = model(2.5,*par0)
    
    # Assert
    assert all(B1 == B2)
    assert all(~np.isnan(B1)) and all(~np.isnan(B2)) and all(~np.isnan(B3)) and all(~np.isnan(B4))
    assert (abs(B5.real - Bref.real) < 1e-5) and (abs(B5.imag - Bref.imag) < 1e-4)
    assert len(paramnames) == len(par0) and len(units) == len(par0)
 #---------------------------------------------------------------------------------------

def test_bg_hom3d():
    Bref = 0.8827946221563685 # Reference from numerical integration
    assert_bgmodel(dl.bg_hom3d,Bref)

def test_bg_hom3d_phase():
    Bref = 0.9998641918107823-0.01648022859583373j # Reference from numerical integration
    assert_bgmodel(dl.bg_hom3d_phase,Bref)

def test_bg_homfractal():
    Bref = 0.9987108379313128 # Reference numerical integration
    assert_bgmodel(dl.bg_homfractal,Bref)

def test_bg_homfractal_phase():
    Bref = 0.999999988356966-0.00015259773164288322j # Reference numerical integration
    assert_bgmodel(dl.bg_homfractal_phase,Bref)

def test_bg_hom3dex():
    Bref = 0.8828975362813715 # Reference from numerical integration
    assert_bgmodel(dl.bg_hom3dex,Bref)

def test_bg_hom3dex_phase():
    Bref = 0.9998641899814665-0.016480339581022234j # Reference from numerical integration
    assert_bgmodel(dl.bg_hom3dex_phase,Bref)

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
