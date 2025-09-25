from __future__ import annotations


from numpy import empty, float64, linspace, pi, unique, zeros
from scipy.sparse.linalg import minres, spsolve

from myfempy.core.solver.assemblerfull import AssemblerFULL
from myfempy.core.solver.assemblerfull_parallel import AssemblerFULLPOOL
from myfempy.core.solver.assemblersymm import AssemblerSYMM
# from myfempy.core.alglin import linsolve_spsolve
from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps


class DynamicHarmonicResponseLinear(Solver):
    """
    Dynamic Harmonic Response Forced System Steady State Linear Solver Class <ConcreteClassService>
    """

    def getMatrixAssembler(Model, inci = None, coord = None, tabmat = None, tabgeo = None, intgauss = None, SYMM=None, MP=None):
        matrix = dict()
        if SYMM:
            matrix["stiffness"] = AssemblerSYMM.getLinearStiffnessGlobalMatrixAssembler(
                Model, inci, coord, tabmat, tabgeo, intgauss,
            )
            matrix["mass"] = AssemblerSYMM.getLinearStiffnessGlobalMatrixAssembler(
                Model, inci, coord, tabmat, tabgeo, intgauss,
            )
        else:
            if MP:
                matrix["stiffness"] = AssemblerFULLPOOL.getLinearStiffnessGlobalMatrixAssembler(
                    Model, inci, coord, tabmat, tabgeo, intgauss,
                    MP=MP,
                )
                matrix["mass"] = AssemblerFULLPOOL.getMassConsistentGlobalMatrixAssembler(
                    Model, inci, coord, tabmat, tabgeo, intgauss,
                    MP=MP,
                )
            else:
                matrix["stiffness"] = AssemblerFULL.getLinearStiffnessGlobalMatrixAssembler(
                    Model, inci, coord, tabmat, tabgeo, intgauss,
                )
                matrix["mass"] = AssemblerFULL.getMassConsistentGlobalMatrixAssembler(
                    Model, inci, coord, tabmat, tabgeo, intgauss,
                )
        return matrix

    def getLoadAssembler(loadaply, nodetot, nodedof):
        return AssemblerFULL.getLoadAssembler(loadaply, nodetot, nodedof)

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
        stiffness = assembly["stiffness"]
        mass = assembly["mass"]
        forcelist = assembly["loads"]

        freedof = constrainsdof["freedof"]

        twopi = 2 * pi
        freqStart = (twopi) * solverset["STEPSET"]["start"]
        freqEnd = (twopi) * solverset["STEPSET"]["end"]
        freqStep = setSteps(solverset["STEPSET"])
        w_range = linspace(freqStart, freqEnd, freqStep)

        U = zeros((fulldofs, freqStep), dtype=float64)
        U0 = U[freedof, 0]

        sA = stiffness[:, freedof][freedof, :]
        sM = mass[:, freedof][freedof, :]
        for ww in range(freqStep):
            Wn = w_range[ww]
            Dw = sA - (Wn**2) * sM
            try:
                U[freedof, ww], info = minres(
                    A=Dw, b=forcelist[freedof, 0], x0=U0, tol=1e-10, maxiter=1000
                )
            except:
                raise info
        solution["U"] = U
        solution["FREQ"] = w_range / (twopi)
        return solution
