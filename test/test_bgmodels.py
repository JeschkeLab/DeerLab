import numpy as np
import deerlab as dl


def assert_bgmodel_value(Bmodel, Bref):
    "Check value at par0 against provided reference value"
    meta = Bmodel.getmetadata()
    par0 = meta["par0"]
    B = Bmodel(2.5, *par0)
    assert (abs(B.real - Bref.real) < 1e-5) and (abs(B.imag - Bref.imag) < 1e-4)


def assert_bgmodel_behavior(model):
    "Check the correct behavior of the core functionality of a background model"

    t = np.linspace(-5, 5, 100)

    # Extract model information
    meta = model.getmetadata()
    par0 = meta["par0"]
    lb = meta["lb"]
    ub = meta["ub"]
    names = meta["names"]
    units = meta["units"]

    # Calculate under different conditions
    B_par0 = model(t, *par0)
    B_par0_T = model(t.T, *par0)
    B_lb = model(t, *lb)
    B_ub = model(t, *ub)

    # Assert
    assert all(B_par0 == B_par0_T)
    assert all(~np.isnan(B_par0))
    assert all(~np.isnan(B_lb))
    assert all(~np.isnan(B_ub))
    assert len(names) == len(par0) and len(units) == len(par0)


def test_bgmodel_behavior():
    models = [
        dl.bg_hom3d,
        dl.bg_homfractal,
        dl.bg_hom3dex,
        dl.bg_hom3dex_phase,
        dl.bg_hom3d_phase,
        dl.bg_homfractal_phase,
        dl.bg_exp,
        dl.bg_strexp,
        dl.bg_prodstrexp,
        dl.bg_sumstrexp,
        dl.bg_poly1,
        dl.bg_poly2,
        dl.bg_poly3,
    ]
    for model in models:
        assert_bgmodel_behavior(model)


def test_bg_hom3d_value():
    t = 1.2
    conc = 143
    lam = 0.47
    B = dl.bg_hom3d(t, conc, lam)
    Bref = 0.92269778443975   # Value from numerical integration
    assert abs(B - Bref) < 1e-4


def test_bg_hom3d_phase_value():
    t = 1.2
    conc = 143
    lam = 0.47
    B = dl.bg_hom3d_phase(t, conc, lam)
    Bref = 0.99994353+0.01062727j   # Value from numerical integration
    assert np.abs(B - Bref) < 1e-4

def test_bg_homfractal_value():
    conc = 1e-3
    dim = 2.2
    lam = 0.47
    t = 0.1
    B = dl.bg_homfractal(t, conc, dim, lam)
    Bref = 0.9443678720378192    # Value from numerical integration
    assert abs(B - Bref) < 1e-4


def test_bg_homfractal_phase():
    conc = 1e-3
    dim = 2.2
    lam = 0.47
    t = 0.1
    B = dl.bg_homfractal_phase(t, conc, dim, lam)
    Bref = 0.9999771819749+0.006755407429425621j    # Value from numerical integration
    assert abs(B - Bref) < 1e-4


def test_bg_hom3dex():
    t = 1.2
    Rex = 3
    conc = 143
    lam = 0.47
    B = dl.bg_hom3dex(t, conc, Rex, lam)
    Bref = 0.926971275363862   # Value from numerical integration
    assert abs(B - Bref) < 1e-4


def test_bg_hom3dex_phase():
    t = 1.2
    Rex = 3
    conc = 143
    lam = 0.47
    B = dl.bg_hom3dex_phase(t, conc, Rex, lam)
    Bref = 0.999943779807736+0.010603642007255394j   # Value from numerical integration
    assert abs(B - Bref) < 1e-4


def test_bg_exp_value():
    k = 0.5656
    t = 2.454
    B = dl.bg_exp(t, k)
    Bref = np.exp(-k * abs(t))
    assert abs(B - Bref) < 1e-10


def test_bg_strexp_value():
    k = 1.3434
    t = 1.445
    xi = 2.7
    B = dl.bg_strexp(t, k, xi)
    Bref = np.exp(-k * abs(t) ** xi)
    assert abs(B - Bref) < 1e-10


def test_bg_prodstrexp_value():
    k1 = 1.3434
    k2 = 0.96767
    t = 1.445
    xi1 = 2.7
    xi2 = 2.4
    B = dl.bg_prodstrexp(t, k1, xi1, k2, xi2)
    Bref = np.exp(-k1 * abs(t) ** xi1) * np.exp(-k2 * abs(t) ** xi2)
    assert abs(B - Bref) < 1e-10


def test_bg_sumstrexp():
    k1 = 1.3434
    k2 = 0.96767
    w = 0.734
    t = 1.445
    xi1 = 2.7
    xi2 = 2.4
    B = dl.bg_sumstrexp(t, k1, xi1, w, k2, xi2)
    Bref = w*np.exp(-k1 * abs(t) ** xi1) + (1-w)*np.exp(-k2 * abs(t) ** xi2)
    assert abs(B-Bref) <  1e-10


def test_bg_poly1_value():
    # Test bg_poly1 model against explicit evaluation
    t = 0.7
    p0 = 1.2
    p1 = -0.23
    B = dl.bg_poly1(t, p0, p1)
    Bref = np.polyval([p1, p0], abs(t))
    assert abs(B - Bref) < 1e-10


def test_bg_poly2_value():
    # Test bg_poly2 model against explicit evaluation
    t = 0.7
    p0 = 1.2
    p1 = -0.23
    p2 = 0.0298
    B = dl.bg_poly2(t, p0, p1, p2)
    Bref = np.polyval([p2, p1, p0], abs(t))
    assert abs(B - Bref) < 1e-10


def test_bg_poly3_value():
    # Test bg_poly3 model against explicit evaluation
    t = 0.7
    p0 = 1.2
    p1 = -0.23
    p2 = 0.0298
    p3 = -0.143
    B = dl.bg_poly3(t, p0, p1, p2, p3)
    Bref = np.polyval([p3, p2, p1, p0], abs(t))
    assert abs(B - Bref) < 1e-10


def test_bg_hom3d_fractal():
    # Make sure bg_homfractal with dim=3 gives the same result as bg_hom3d
    conc = 100
    dim = 3
    lam = 0.4
    t = 0.789
    B1 = dl.bg_hom3d(t, conc, lam)
    B2 = dl.bg_homfractal(t, conc, dim, lam)
    assert abs(B1 - B2) < 1e-6
