# __init__.py
from importlib.metadata import version as get_version

try:
    __VERSION__ = get_version('DeerLab')
except Exception:
    __VERSION__ = 'unknown'

from .dd_models import *
from .bg_models import *
from .model import Model, Penalty, Parameter, link, lincombine, merge, relate
from .deerload import deerload
from .selregparam import selregparam
from .dipolarkernel import dipolarkernel
from .dipolarbackground import dipolarbackground
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
from .io import save, load, json_dumps, json_loads