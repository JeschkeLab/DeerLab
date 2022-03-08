import numpy as np
from deerlab import time2dist

def test_time2dist_outsize():
    "Verify that the output has the correct size"

    t = np.linspace(-1,5,300)
    rminmax = time2dist(t)
    assert len(rminmax)==2

    nr = 200
    r = time2dist(t,nr)
    assert len(r)==nr


def test_time2dist_rmin():
    "Verify that the minimum distance is correct"

    t, dt = np.linspace(-1, 5, 300, retstep=True)
    rmin, _ = time2dist(t)

    D = 52.04  # MHz nm^3
    nu_max = 1/dt/2/2*0.85  # Nyquist freq, with buffer
    rmin_ref = (D/nu_max)**(1/3)

    assert np.max(np.abs(rmin - rmin_ref)) < 1e-10


def test_time2dist_rmax():
    "Verify that the maximum distance is correct"

    t, dt = np.linspace(-1, 5, 300, retstep=True)
    _, rmax = time2dist(t)

    trange = max(t) - min(t)
    D = 52.04  # MHz nm^3
    Tmax = trange*2  # maximum period length
    rmax_ref = (D*Tmax)**(1/3)

    assert np.max(np.abs(rmax - rmax_ref)) < 1e-10
