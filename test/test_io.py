import numpy as np
import pytest
import os
import tempfile
import deerlab as dl
from deerlab import dd_gauss, bg_hom3d, FitResult
from deerlab.dipolarmodel import dipolarmodel
from deerlab.fit import fit
from deerlab.io import save, load, json_dumps, json_loads
from matplotlib.figure import Figure
import h5py
import io

from test_uqresult import uncertainty_quantification_simulation

# Shared test data
t = np.linspace(-0.5, 5, 100)
r = np.linspace(2, 5, 50)

@pytest.fixture(scope='module')
def fitresult():
    Bfcn = lambda t, lam: bg_hom3d(t, 50, lam)
    Pr = dd_gauss(r, 3, 0.2)
    V = 1e5 * dl.dipolarkernel(t, r, mod=0.3, bg=Bfcn) @ Pr
    model = dipolarmodel(t, r, dd_gauss, bg_hom3d, npathways=1)
    return fit(model, V, ftol=1e-4)

# ======================================================================
def test_fitresult_to_dict(fitresult):
    "Check that FitResult can be converted to a dictionary"

    d = fitresult.to_dict()

    assert isinstance(d, dict)
    assert 'model' in d
# ======================================================================

def test_fitresult_save_filebuffer(fitresult):
    "Check that FitResult can be saved to an HDF5 fileIO buffer"


    buf = io.BytesIO()
    print(buf.getbuffer().nbytes)
    save(buf, fitresult, format='hdf5')
    assert buf.getbuffer().nbytes > 0
    

# ======================================================================

def test_fitresult_save_hdf5(fitresult):
    "Check that FitResult can be saved to an HDF5 file"

    with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
        fname = f.name
    try:
        save(fname, fitresult)
        assert os.path.exists(fname)
        # Check that the file can be opened and has the expected attributes
        with h5py.File(fname, 'r') as file:
            assert file.attrs['format'] == 'deerlab'
            assert file.attrs['version'] == dl.__VERSION__
            assert file.attrs['object_class'] == 'FitResult'
            assert 'model' in file.keys()
    finally:
        os.remove(fname)


# ======================================================================
def test_fitresult_save_load_hdf5(fitresult):
    "Check that FitResult can be saved and loaded from an HDF5 file"

    with tempfile.NamedTemporaryFile(suffix='.hdf5', delete=False) as f:
        fname = f.name
    try:
        save(fname, fitresult)
        loaded = load(fname)
        assert isinstance(loaded, FitResult)
        assert isinstance(loaded.modelUncert,dl.UQResult)
        assert np.allclose(loaded.model, fitresult.model)
        assert np.allclose(loaded.regparam, fitresult.regparam)
        assert isinstance(loaded.plot(), Figure)
        
    finally:
        os.remove(fname)

# ======================================================================
def test_fitresult_to_from_json_string(fitresult):
    "Check that FitResult can be converted to a JSON string"

    json_str = json_dumps(fitresult)
    assert isinstance(json_str, str)
    assert len(json_str) > 0
    fitresult_loaded = json_loads(json_str)
    assert isinstance(fitresult_loaded, FitResult)
    assert np.allclose(fitresult_loaded.model, fitresult.model)
    assert isinstance(fitresult_loaded.modelUncert,dl.UQResult)
    assert np.allclose(fitresult_loaded.regparam, fitresult.regparam)

# ======================================================================
def test_fitresult_save_load_json(fitresult):
    "Check that FitResult can be saved and loaded from a JSON file"

    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        fname = f.name
    try:
        save(fname, fitresult)
        loaded = load(fname)
        assert isinstance(loaded, FitResult)
        assert np.allclose(loaded.model, fitresult.model)
        assert isinstance(loaded.modelUncert,dl.UQResult)
        assert np.allclose(loaded.regparam, fitresult.regparam)
        assert isinstance(loaded.plot(), Figure)
        

    finally:
        os.remove(fname)


# ======================================================================
def test_fitresult_save_load_toml(fitresult):
    "Check that FitResult can be saved and loaded from a TOML file"

    with tempfile.NamedTemporaryFile(suffix='.toml', delete=False) as f:
        fname = f.name
    try:
        save(fname, fitresult)
        loaded = load(fname)
        assert isinstance(loaded, FitResult)
        assert np.allclose(loaded.model, fitresult.model)
        assert isinstance(loaded.modelUncert,dl.UQResult)
        assert np.allclose(loaded.regparam, fitresult.regparam)
        assert isinstance(loaded.plot(), Figure)

    finally:
        os.remove(fname)


# ======================================================================

@pytest.mark.parametrize('method', ['moment', 'bootstrap', 'profile'])
def test_UQResult_saving(uncertainty_quantification_simulation, method):
    "Check that UQResult can be saved and loaded from a JSON string"

    # Retrieve the results of the mock simulation
    uq_objects, references = uncertainty_quantification_simulation
    uq = uq_objects[method]

    json_str = json_dumps(uq)
    assert isinstance(json_str, str)
    assert len(json_str) > 0
    UQresult_loaded = json_loads(json_str)
    assert isinstance(UQresult_loaded, dl.UQResult)
    assert np.allclose(UQresult_loaded.mean, uq.mean)
    assert UQresult_loaded.type == uq.type
    assert np.allclose(UQresult_loaded.std, uq.std)
    assert np.allclose(UQresult_loaded.lb, uq.lb)
    assert np.allclose(UQresult_loaded.ub, uq.ub)