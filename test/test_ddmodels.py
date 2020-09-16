import numpy as np
import deerlab as dl

def assert_ddmodel(model):
#============================================================================
    "Check dimensionality, non-negativity, and values at boundaries of a distance distribution model"

    r = np.linspace(0,20,500)
    rnus = np.sqrt(np.linspace(1.5,6**2,500))

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
    passed[4] = np.round(np.trapz(P5,rnus),2) == 1
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
#=============================================================================

def test_dd_gauss():
    assert_ddmodel(dl.dd_gauss)

def test_dd_gauss2():
    assert_ddmodel(dl.dd_gauss2)

def test_dd_gauss3():
    assert_ddmodel(dl.dd_gauss3)

def test_dd_gengauss():
    assert_ddmodel(dl.dd_gengauss)
    
def test_dd_skewgauss():
    assert_ddmodel(dl.dd_skewgauss)

def test_dd_rice():
    assert_ddmodel(dl.dd_rice)

def test_dd_rice2():
    assert_ddmodel(dl.dd_rice2)

def test_dd_rice3():
    assert_ddmodel(dl.dd_rice3)

def test_dd_randcoil():
    assert_ddmodel(dl.dd_randcoil)

def test_dd_circle():
    assert_ddmodel(dl.dd_circle)

def test_dd_cos():
    assert_ddmodel(dl.dd_cos)

def test_dd_sphere():
    assert_ddmodel(dl.dd_sphere)

def test_dd_shell():
    assert_ddmodel(dl.dd_shell)

def test_dd_shellshell():
    assert_ddmodel(dl.dd_shellshell)

def test_dd_shellsphere():
    assert_ddmodel(dl.dd_shellsphere)

def test_dd_shellvoidsphere():
    assert_ddmodel(dl.dd_shellvoidsphere)

def test_dd_shellvoidshell():
    assert_ddmodel(dl.dd_shellvoidshell)

def test_dd_spherepoint():
    assert_ddmodel(dl.dd_spherepoint)

def test_dd_spheresurf():
    assert_ddmodel(dl.dd_spheresurf)

def test_dd_triangle():
    assert_ddmodel(dl.dd_triangle)

def test_dd_uniform():
    assert_ddmodel(dl.dd_uniform)

def test_dd_wormchain():
    assert_ddmodel(dl.dd_wormchain)

def test_dd_wormgauss():
    assert_ddmodel(dl.dd_wormgauss)