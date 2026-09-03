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
from .core.solver.homogenplaneinfperiodic import HomogenizationPlaneBCPeriodic
from .core.solver.phonocrystalinplane import PhononicCrystalPlaneBCPeriodic
# API CLASS
from .core.elements.element import Element
from .core.shapes.shape import Shape
from .core.mesh.mesh import Mesh
from .core.material.material import Material
from .core.geometry.geometry import Geometry
from .core.physic.physics import Physics
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
__author__ = 'Antonio Vinicius Garcia Campos'
__license__ = "GPLv3"
__email__ = "antviniciuscampos@gmail.com"

__docformat__ = "google"
__doc__ = """
myfempy
=======

myfempy is a python package based on finite element method to multiphysics
analysis. The code is open source and *intended for educational and scientific
purposes only, not recommended to commercial use. The name myfempy is an acronym
for MultiphYsics Finite Elements Module to PYthon. You can help us by contributing
with the main project, send us a mensage on https://github.com/avgcampos/myfempy/discussions/10
If you use myfempy in your research, the  developers would be grateful if you 
could cite in your work.
						
Web Documentation:
------------------
For examples and a description of the API, please visit the online documentation.: https://myfempy.readthedocs.io

Example:
--------
from myfempy import newAnalysis, SteadyStateLinear
FEA = newAnalysis(SteadyStateLinear)
FEA.Model(modeldata:dict)
FEA.Physic(physicdata:dict)
FEA.Solve(solverset:dict)
FEA.PostProcess(postprocset:dict)

"""

__all__ = [
    "__version__",
    "get_about"
    "newAnalysis",
    "SteadyStateLinear",
    "SteadyStateLinearIterative",
    "StaticLinearCyclicSymmPlane",
    "DynamicEigenLinear",
    "DynamicHarmonicResponseLinear",
    "HomogenizationPlane",
    "HomogenizationPlaneBCPeriodic",
    "PhononicCrystalPlaneBCPeriodic",
    "Mesh",
    "Shape",
    "Element",
    "Geometry",
    "Material",
    "Physics",
    "PlaneStress",
    "PlaneStrain",
    "UniAxialStress",
    "SolidElastic",
    "HeatPlane",
    "HeatSolid"
]
