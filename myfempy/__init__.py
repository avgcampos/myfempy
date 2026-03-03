from __future__ import absolute_import
# API
from .api.main import newAnalysis
# SOLVERS
from .core.solver.steadystatelinear import SteadyStateLinear
from .core.solver.steadystatelineariterative import SteadyStateLinearIterative
from .core.solver.dyneigen import DynamicEigenLinear
from .core.solver.dynharmonicresponse import DynamicHarmonicResponseLinear
from .core.solver.cyclicsymm import StaticLinearCyclicSymmPlane
from .core.solver.homogenplanefullcell import HomogenizationPlane
# from .core.solver.homogenplaneinfperiodic import HomogenPlaneInfPeriodic
from .core.solver.phonocrystalinplane import PhononicCrystalInPlane
# API CLASS
from .core.elements.element import Element
from .core.geometry.geometry import Geometry
from .core.shapes.shape import Shape
from .core.mesh.mesh import Mesh
from .core.material.material import Material
# API CLASS SET MATERIAL
from .core.material.planestress import PlaneStress
from .core.material.planestrain import PlaneStrain
from .core.material.uniaxialstress import UniAxialStress
from .core.material.solidelastic import SolidElastic
from .core.material.heatplane import HeatPlane
from .core.material.heatsolid import HeatSolid
# VERSION
from .utils.utils import get_version
__version__ = get_version()

__all__ = [
    "__version__",
    "newAnalysis",
    "SteadyStateLinear",
    "SteadyStateLinearIterative",
    "StaticLinearCyclicSymmPlane",
    "DynamicEigenLinear",
    "DynamicHarmonicResponseLinear",
    "PhononicCrystalInPlane",
    "HomogenizationPlane",
    # "HomogenPlaneInfPeriodic"
    "Mesh",
    "Shape",
    "Element",
    "Geometry",
    "Material",
    "PlaneStress",
    "PlaneStrain",
    "UniAxialStress",
    "SolidElastic",
    "HeatPlane",
    "HeatSolid"
]
