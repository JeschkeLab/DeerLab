# __init__.py
from . import dd_models as _dd_models_mod
from . import bg_models as _bg_models_mod

# Define __getattr__ early so submodules that do `from deerlab import bg_*`
# during their own import (e.g. dipolarmodel) can resolve names via this hook.
def __getattr__(name):
    if name in _dd_models_mod.__all__:
        return _dd_models_mod.__getattr__(name)
    if name in _bg_models_mod.__all__:
        return _bg_models_mod.__getattr__(name)
    raise AttributeError(f"module 'deerlab' has no attribute {name!r}")

from .model import Model, Penalty, Parameter, link, lincombine, merge, relate
from .deerload import deerload
from .selregparam import selregparam
from .dipolarkernel import dipolarkernel
from .dipolarbackground import dipolarbackground,dipolarbackgroundmodel
from .dipolarmodel import dipolarmodel,ExperimentInfo, dipolarpenalty, ex_4pdeer,ex_3pdeer,ex_rev5pdeer,ex_fwd5pdeer,ex_ridme,ex_sifter,ex_dqc
from .solvers import snlls, fnnls, cvxnnls
from .regoperator import regoperator
from .correctphase import correctphase
from .correctzerotime import correctzerotime
from .noiselevel import noiselevel
from .whitegaussnoise import whitegaussnoise
from .bootstrap_analysis import bootstrap_analysis
from .fftspec import fftspec
from .distancerange import distancerange
from .classes import UQResult
from .fit import fit
from .fitresult import FitResult
from .diststats import diststats
from .profile_analysis import profile_analysis
from .utils import *
