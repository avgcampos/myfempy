from __future__ import absolute_import

from .setup.fea import newAnalysis
from .core.solver.staticlinear import StaticLinear
from .core.solver.staticlineariterative import StaticLinearIterative
from .core.solver.modallinear import ModalLinear
from .core.solver.harmoniclinear import HarmonicLinear
from .__about__ import __version__

__all__ = ["newAnalysis",
           "StaticLinear",
           "StaticLinearIterative",
           "ModalLinear",
           "HarmonicLinear",
           "__version__"]

