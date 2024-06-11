from __future__ import absolute_import

from .__about__ import __version__
from .core.solver.harmoniclinear import HarmonicLinear
from .core.solver.modallinear import ModalLinear
from .core.solver.staticlinear import StaticLinear
from .core.solver.staticlineariterative import StaticLinearIterative
from .setup.fea import newAnalysis

__all__ = [
    "newAnalysis",
    "StaticLinear",
    "StaticLinearIterative",
    "ModalLinear",
    "HarmonicLinear",
    "__version__",
]
