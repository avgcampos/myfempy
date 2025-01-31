from __future__ import absolute_import

from .core.solver.cyclicsymm import StaticLinearCyclicSymm
from .core.solver.dynsteadystatelinear import DynamicSteadyStateLinear
from .core.solver.dynmodal import DynamicModalLinear
from .core.solver.steadystatelineariterative import SteadyStateLinearIterative
from .core.solver.steadystatelinear import SteadyStateLinear
from .setup.fea import newAnalysis

from .utils.utils import get_version
__version__ = get_version()

__all__ = [
    "__version__",
    "newAnalysis",
    "SteadyStateLinear",
    "SteadyStateLinearIterative",
    "StaticLinearCyclicSymm",
    "DynamicModalLinear",
    "DynamicSteadyStateLinear",
]
