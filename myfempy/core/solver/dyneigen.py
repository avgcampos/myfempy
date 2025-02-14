from __future__ import annotations


from numpy import (arange, concatenate, empty, float64, newaxis, pi, sqrt,
                   unique, zeros)
from scipy.sparse.linalg import eigsh

from myfempy.core.solver.assemblerfull import AssemblerFULL
from myfempy.core.solver.assemblerfull_parallel import AssemblerFULLPOOL
from myfempy.core.solver.assemblersymm import AssemblerSYMM
from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps


class DynamicEigenLinear(Solver):
    """
    Dynamic Eigen (modal problem) Linear Solver Class <ConcreteClassService>
    """

    def getMatrixAssembler(
        Model, inci, coord, tabmat, tabgeo, intgauss, SYMM=None, MP=None
    ):
        matrix = dict()
        if SYMM:
            matrix["stiffness"] = AssemblerSYMM.getLinearStiffnessGlobalMatrixAssembler(
                Model,
                inci,
                coord,
                tabmat,
                tabgeo,
                intgauss,
                type_assembler="linear_stiffness",
                MP=MP,
            )
            matrix["mass"] = AssemblerSYMM.getMassConsistentGlobalMatrixAssembler(
                Model,
                inci,
                coord,
                tabmat,
                tabgeo,
                intgauss,
                type_assembler="mass_consistent",
                MP=MP,
            )
        else:
            if MP:
                matrix["stiffness"] = AssemblerFULLPOOL.getLinearStiffnessGlobalMatrixAssembler(
                    Model,
                    inci,
                    coord,
                    tabmat,
                    tabgeo,
                    intgauss,
                    type_assembler="linear_stiffness",
                    MP=MP,
                )
                matrix["mass"] = AssemblerFULLPOOL.getMassConsistentGlobalMatrixAssembler(
                    Model,
                    inci,
                    coord,
                    tabmat,
                    tabgeo,
                    intgauss,
                    type_assembler="linear_stiffness",
                    MP=MP,
                )
            else:
                matrix["stiffness"] = AssemblerFULL.getLinearStiffnessGlobalMatrixAssembler(
                    Model,
                    inci,
                    coord,
                    tabmat,
                    tabgeo,
                    intgauss,
                    type_assembler="linear_stiffness",
                    MP=MP,
                )
                matrix["mass"] = AssemblerFULL.getMassConsistentGlobalMatrixAssembler(
                    Model,
                    inci,
                    coord,
                    tabmat,
                    tabgeo,
                    intgauss,
                    type_assembler="linear_stiffness",
                    MP=MP,
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

    def runSolve(assembly, constrainsdof, modelinfo, solverset):
        fulldofs = modelinfo["fulldofs"]
        solution = dict()
        modeEnd = setSteps(solverset["STEPSET"])
        stiffness = assembly["stiffness"]
        mass = assembly["mass"]
        forcelist = assembly["loads"]
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
