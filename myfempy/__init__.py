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

__docformat__ = "google"

__doc__ = """

==========================================================================
                            __                                
         _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
        | '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
        | | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
        |_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
                    |___/                       |_|     |___/ 
        myfempy -- MultiphYsics Finite Element Module to PYthon    
                    COMPUTATIONAL ANALYSIS PROGRAM                   
        Copyright (C) 2022-2026 Antonio Vinicius Garcia Campos        
==========================================================================
This Python file is part of myfempy project.

myfempy is a python package based on finite element method to multiphysics
analysis. The code is open source and *intended for educational and scientific
purposes only, not recommended to commercial use. The name myfempy is an acronym
for MultiphYsics Finite Elements Module to PYthon. You can help us by contributing
with the main project, send us a mensage on https://github.com/avgcampos/myfempy/discussions/10
If you use myfempy in your research, the  developers would be grateful if you 
could cite in your work.
																		
The code is written by Antonio Vinicius Garcia Campos.                                  
																		
A github repository, with the most up to date version of the code,      
can be found here: https://github.com/avgcampos/myfempy.                 
																		
The code is open source and intended for educational and scientific     
purposes only. If you use myfempy in your research, the developers      
would be grateful if you could cite this. The myfempy project is published
under the GPLv3, see the myfempy LICENSE on
https://github.com/avgcampos/myfempy/blob/main/LICENSE.
																		
Disclaimer:                                                             
The authors reserve all rights but do not guarantee that the code is    
free from errors. Furthermore, the authors shall not be liable in any   
event caused by the use of the program.

"""

__all__ = [
    "__version__",
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
