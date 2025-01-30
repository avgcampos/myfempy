from __future__ import annotations

# from os import environ
# environ['OMP_NUM_THREADS'] = '3'
from numpy import dot, float64, zeros
from scipy.sparse.linalg import minres

from myfempy.core.solver.assemblerfull import AssemblerFULL
from myfempy.core.solver.assemblersymm import AssemblerSYMM
from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps


class SteadyStateLinearIterative(Solver):

    """
    Steady State Linear Iterative Solver Class <ConcreteClassService>
    """

    def getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, SYMM, MP):
        
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
        return matrix

    def getLoadAssembler(loadaply, nodetot, nodedof):
        return AssemblerFULL.getLoadAssembler(loadaply, nodetot, nodedof)

    def getConstrains(constrains, nodetot, nodedof):
        return AssemblerFULL.getConstrains(constrains, nodetot, nodedof)

    def getDirichletNH(constrains, nodetot, nodedof):
        return AssemblerFULL.getDirichletNH(constrains, nodetot, nodedof)

    def runSolve(assembly, constrainsdof, modelinfo, solverset):
        fulldofs = modelinfo["fulldofs"]

        solution = dict()
        nsteps = setSteps(solverset["STEPSET"])

        stiffness = assembly["stiffness"]
        forcelist = assembly["loads"]

        U0 = zeros((fulldofs), dtype=float64)  # empty((fulldofs, 1))
        U1 = zeros((fulldofs), dtype=float64)  # empty((fulldofs, 1))
        U = zeros((fulldofs, nsteps), dtype=float64)  # empty((fulldofs, nsteps))
        Uc = assembly["bcdirnh"]

        freedof = constrainsdof["freedof"]
        constdof = constrainsdof["constdof"]

        for step in range(nsteps):
            forcelist[freedof, step] = forcelist[freedof, step] - dot(
                stiffness[:, constdof][freedof, :].toarray(), Uc[constdof, step]
            )
            try:
                U1[freedof], info = minres(
                    A=stiffness[:, freedof][freedof, :],
                    b=forcelist[freedof, step],
                    tol=1e-10,
                    maxiter=1000,
                )
            except:
                raise info
            U1[constdof] = Uc[constdof, step]
            U1[:] += U0[:]
            U[:, step] = U1
            U0[:] = U1[:]
        solution["U"] = U
        return solution
