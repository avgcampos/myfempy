from __future__ import absolute_import

from .setup.feaanalysis import newAnalysis
from .core.solver.staticlinear import StaticLinear
from .core.solver.modallinear import ModalLinear
from .core.solver.harmoniclinear import HarmonicLinear
from .__about__ import __version__

__all__ = ["newAnalysis",
           "StaticLinear",
           "ModalLinear",
           "HarmonicForced",
           "__version__"]

