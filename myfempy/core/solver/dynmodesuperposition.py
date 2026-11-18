from __future__ import annotations

from numpy import empty, float64, linspace, pi, unique, zeros, sqrt, nonzero, zeros_like, broadcast_to, maximum
from scipy.sparse.linalg import spsolve, eigsh

from myfempy.core.solver.assemblerSR import AssemblerSR
from myfempy.core.solver.assemblerMP import AssemblerMP
# from myfempy.core.solver.dyneigen import DynamicEigenLinear
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


class DynamicModeSuperposition(Solver):
    """
    Dynamic Mode Superposition Response System Steady State Solver Class <ConcreteClassService>
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
        return zeros(
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
        modeEnd = solverset["STEPSET"]["modeEnd"]

        try:
            fdam = solverset["DAMP"]
        except:
            fdam = 1E-30

        U = zeros((fulldofs, freqStep), dtype=float64)

        sK = stiffness[:, freedof][freedof, :]
        sM = mass[:, freedof][freedof, :]
        sF = forcelist[freedof, 0]
        f_in_dof = nonzero(sF)

        WN_rad, PHI = eigsh(
            A=sK,
            M=sM,
            k=modeEnd,
            sigma=1E-9,
            which="LM",
        )
        WF_rad = linspace(freqStart, freqEnd, freqStep)
        WNmin = maximum(WN_rad, 1E-9)
        WR = sqrt(WNmin)
        zeta = broadcast_to(fdam, WN_rad.shape)
        PHI_f_in = PHI[f_in_dof[0][0], :]
        U0 = spsolve(sK, sF)
        R = zeros_like(U0)
        if solverset.get("RFSC"):
            R = U0 - PHI @ (PHI_f_in/WN_rad)
        
        for ww in range(freqStep):
            try:
                aux = PHI_f_in/(WN_rad - WF_rad[ww]**2 + 2*img*zeta*WR*WF_rad[ww])
                U[freedof, ww] = PHI @ aux + R
            except:
                raise RuntimeError(f"solver não convergiu no passo {ww})")
        solution["U"] = U
        solution["FREQ"] = (WF_rad / (twopi)).reshape(-1, 1)
        return solution
