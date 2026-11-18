from __future__ import annotations

from numpy import empty, float64, linspace, pi, unique, zeros, array
from numpy.linalg import solve
from scipy.sparse.linalg import spsolve

from myfempy.core.solver.assemblerSR import AssemblerSR
from myfempy.core.solver.assemblerMP import AssemblerMP
from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps

__docformat__ = "google"

__doc__ = """
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


class DynamicHarmonicResponseLinear(Solver):
    """
    Dynamic Harmonic Response Forced System Steady State Linear Solver Class <ConcreteClassService>
    """

    def getMatrixAssembler(
        Model,
        inci=None,
        coord=None,
        tabmat=None,
        tabgeo=None,
        intgauss=None,
        MP=None,
        max_workers=None,
    ):
        matrix = dict()
        if MP:
            matrix["stiffness"] = AssemblerSR.getGlobalMatrixAssemblerMP(
                Model,
                Model.element.getStifLinearMat,
                Model.shape.getIntNumK,
                inci,
                coord,
                tabmat,
                tabgeo,
                intgauss,
                max_workers,
            )
            matrix["mass"] = AssemblerMP.getGlobalMatrixAssembler(
                Model,
                Model.element.getMassConsistentMat,
                Model.shape.getIntNumM,
                inci,
                coord,
                tabmat,
                tabgeo,
                intgauss,
                max_workers,
            )
        else:
            matrix["stiffness"] = AssemblerSR.getGlobalMatrixAssembler(
                Model,
                Model.element.getStifLinearMat,
                Model.shape.getIntNumK,
                inci,
                coord,
                tabmat,
                tabgeo,
                intgauss,
            )
            matrix["mass"] = AssemblerSR.getGlobalMatrixAssembler(
                Model,
                Model.element.getMassConsistentMat,
                Model.shape.getIntNumM,
                inci,
                coord,
                tabmat,
                tabgeo,
                intgauss,
            )
        return matrix

    def getLoadAssembler(loadaply, nodetot, nodedof):
        return AssemblerSR.getLoadAssembler(loadaply, nodetot, nodedof)

    def getConstrains(constrains, nodetot, nodedof):
        return AssemblerSR.getConstrains(constrains, nodetot, nodedof)

    def getDirichletNH(constrains, nodetot, nodedof):
        return empty(
            (nodedof * nodetot, len(unique(constrains[:, 3][constrains[:, 3] != 0]))),
            dtype=float64,
        )

    def runSolve(Model, Physic, assembly, constrainsdof, solverset):
        fulldofs = Model.modelinfo["fulldofs"]

        solution = dict()
        stiffness = assembly["stiffness"]
        mass = assembly["mass"]
        forcelist = assembly["loads"]

        freedof = constrainsdof["freedof"]

        twopi = 2 * pi
        img = 1j
        freqStart = twopi * solverset["STEPSET"]["start"]
        freqEnd = twopi * solverset["STEPSET"]["end"]
        freqStep = setSteps(solverset["STEPSET"])
        w_range = linspace(freqStart, freqEnd, freqStep)
        
        try:
            fdam = solverset["DAMP"]
        except:
            fdam = 1E-30

        U = zeros((fulldofs, freqStep), dtype=float64)

        alpha, beta = DynamicHarmonicResponseLinear.__rayleigh_from_zeta(freqStart, freqEnd, fdam, fdam)

        sK = stiffness[:, freedof][freedof, :]
        sM = mass[:, freedof][freedof, :]
        sC = alpha * sM + beta * sK
        sF = forcelist[freedof, 0]
        for ww in range(freqStep):
            Zw = sK - (w_range[ww]**2) * sM + (img * w_range[ww]) * sC
            try:
                U[freedof, ww] = spsolve(A=Zw, b=sF)
            except:
                raise RuntimeError(f"spsolve não convergiu no passo {ww})")
        solution["U"] = U
        solution["FREQ"] = (w_range / (twopi)).reshape(-1, 1)
        return solution

    def __rayleigh_from_zeta(w1, w2, z1, z2):
        """alpha, beta para zeta = z1 em w1 e z2 em w2 (resolve o sistema 2x2)."""
        A = 0.5*array([[1/(w1+1E-30), w1], [1/w2, w2]])
        return solve(A, [z1, z2])