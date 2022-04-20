import numpy as np
from deerlab import distancerange

def test_distancerange_outsize():
    "Verify that the output has the correct size"

    t = np.linspace(-1,5,300)
    rminmax = distancerange(t)
    assert len(rminmax)==2

    nr = 200
    r = distancerange(t,nr)
    assert len(r)==nr


def test_distancerange_rmin():
    "Verify that the minimum distance is correct"

    t, dt = np.linspace(-1, 5, 300, retstep=True)
    rmin, _ = distancerange(t)

    D = 52.04  # MHz nm^3
    nu_max = 1/dt/2/2*0.85  # Nyquist freq, with buffer
    rmin_ref = (D/nu_max)**(1/3)

    assert np.max(np.abs(rmin - rmin_ref)) < 1e-10


def test_distancerange_rmax():
    "Verify that the maximum distance is correct"

    t, dt = np.linspace(-1, 5, 300, retstep=True)
    _, rmax = distancerange(t)

    trange = max(t) - min(t)
    D = 52.04  # MHz nm^3
    Tmax = trange*2  # maximum period length
    rmax_ref = (D*Tmax)**(1/3)

    assert np.max(np.abs(rmax - rmax_ref)) < 1e-10
