from __future__ import annotations


from numpy import (array, arange, concatenate, empty, float64, newaxis, pi, sqrt,
                   unique, zeros)
from scipy.sparse.linalg import eigsh

from myfempy.core.solver.assemblerfull import AssemblerFULL
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


class DynamicEigenLinear(Solver):
    """
    Dynamic Eigen (modal problem) Linear Solver Class <ConcreteClassService>
    """
    def getMatrixAssembler(Model, inci = None, coord = None, tabmat = None, tabgeo = None, intgauss = None, MP = None):
        matrix = dict()
        if MP:
            matrix["stiffness"] = AssemblerFULL.getGlobalMatrixAssemblerMP(
                        Model, Model.element.getStifLinearMat, Model.shape.getIntNumK, inci, coord, tabmat, tabgeo, intgauss,

                    )
            matrix["mass"] = AssemblerFULL.getGlobalMatrixAssemblerMP(
                Model, Model.element.getMassConsistentMat, Model.shape.getIntNumM, inci, coord, tabmat, tabgeo, intgauss,
            )
        else:
            matrix["stiffness"] = AssemblerFULL.getGlobalMatrixAssembler(
                Model, Model.element.getStifLinearMat, Model.shape.getIntNumK, inci, coord, tabmat, tabgeo, intgauss,

            )
            matrix["mass"] = AssemblerFULL.getGlobalMatrixAssembler(
                Model, Model.element.getMassConsistentMat, Model.shape.getIntNumM, inci, coord, tabmat, tabgeo, intgauss,
            )
        return matrix

    def getLoadAssembler(loadaply, nodetot, nodedof):
        return empty(
            (
                nodedof * nodetot,
                len(unique(loadaply[:, 3])),
            )
        )

    def getConstrains(constrains, nodetot, nodedof):
        return AssemblerFULL.getConstrains(constrains, nodetot, nodedof)

    def getDirichletNH(constrains, nodetot, nodedof):
        return empty(
            (nodedof * nodetot, len(unique(constrains[:, 3][constrains[:, 3] != 0]))),
            dtype=float64,
        )

    def runSolve(Model, Physic, assembly, constrainsdof, solverset):
        fulldofs = Model.modelinfo["fulldofs"]
        solution = dict()
        modeEnd = setSteps(solverset["STEPSET"])
        stiffness = assembly["stiffness"]
        mass = assembly["mass"]
        # forcelist = assembly["loads"]
        U = zeros((fulldofs, modeEnd), dtype=float64)
        freedof = constrainsdof["freedof"]
        try:
            W, U[freedof, :] = eigsh(
                A=stiffness[:, freedof][freedof, :],
                M=mass[:, freedof][freedof, :],
                k=modeEnd,
                sigma=1,
                which="LM",
                maxiter=1000,
            )
        except:
            pass
        Wlist = arange(0, modeEnd + 1)
        Wrad = sqrt(W)
        Whz = Wrad / (2 * pi)
        w_range = concatenate(
            (Wlist[1:, newaxis], Wrad[:, newaxis], Whz[:, newaxis]), axis=1
        )
        solution["U"] = U
        solution["FREQ"] = w_range
        return solution
