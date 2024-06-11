from __future__ import annotations

from os import environ

environ["OMP_NUM_THREADS"] = "3"

from numpy import empty, float64, linspace, pi, unique, zeros
from scipy.sparse.linalg import minres, spsolve

from myfempy.core.solver.assemblerfull import AssemblerFULL
from myfempy.core.solver.assemblersymm import AssemblerSYMM
# from myfempy.core.alglin import linsolve_spsolve
from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps


class HarmonicLinear(Solver):

    """
    Harmonic Forced System Linear Solver Class <ConcreteClassService>
    """

    def getMatrixAssembler(
        Model, inci, coord, tabmat, tabgeo, intgauss, SYMM=None, MP=None
    ):
        matrix = dict()
        if SYMM:
            matrix["stiffness"] = AssemblerSYMM.getMatrixAssembler(
                Model,
                inci,
                coord,
                tabmat,
                tabgeo,
                intgauss,
                type_assembler="linear_stiffness",
                MP=MP,
            )
            matrix["mass"] = AssemblerSYMM.getMatrixAssembler(
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
            matrix["stiffness"] = AssemblerFULL.getMatrixAssembler(
                Model,
                inci,
                coord,
                tabmat,
                tabgeo,
                intgauss,
                type_assembler="linear_stiffness",
                MP=MP,
            )
            matrix["mass"] = AssemblerFULL.getMatrixAssembler(
                Model,
                inci,
                coord,
                tabmat,
                tabgeo,
                intgauss,
                type_assembler="mass_consistent",
                MP=MP,
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

    def runSolve(assembly, constrainsdof, fulldofs, solverset):
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
