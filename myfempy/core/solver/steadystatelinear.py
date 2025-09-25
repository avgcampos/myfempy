from __future__ import annotations


from numpy import dot, float64, zeros
from scipy.sparse.linalg import spsolve
import scipy.sparse as sp

from myfempy.core.solver.assemblerfull import AssemblerFULL
from myfempy.core.solver.assemblerfull_parallel import AssemblerFULLPOOL
from myfempy.core.solver.assemblersymm import AssemblerSYMM
# from myfempy.core.alglin import linsolve_spsolve
from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps

class SteadyStateLinear(Solver):
    """
    Steady State Linear Solver Class <ConcreteClassService>
    """

    # @profile
    def getMatrixAssembler(Model, inci = None, coord = None, tabmat = None, tabgeo = None, intgauss = None, SYMM=None, MP=None):
       
        matrix = dict()
        
        if SYMM:
            assembler = AssemblerSYMM.getLinearStiffnessGlobalMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss)
        else:
            if MP:
                assembler = AssemblerFULLPOOL.getLinearStiffnessGlobalMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, MP)
            else:
                assembler = AssemblerFULL.getLinearStiffnessGlobalMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss)
        
        matrix["stiffness"] = assembler #sp.load_npz('sparse_matrix.npz')
        
        return matrix

    def getLoadAssembler(loadaply, nodetot, nodedof):
        return AssemblerFULL.getLoadAssembler(loadaply, nodetot, nodedof)

    def getConstrains(constrains, nodetot, nodedof):
        return AssemblerFULL.getConstrains(constrains, nodetot, nodedof)

    def getDirichletNH(constrains, nodetot, nodedof):
        return AssemblerFULL.getDirichletNH(constrains, nodetot, nodedof)

    def runSolve(Model, Physic, assembly, constrainsdof, solverset):
        fulldofs = Model.modelinfo["fulldofs"]

        solution = dict()
        nsteps = setSteps(solverset["STEPSET"])

        stiffness = assembly["stiffness"]
        forcelist = assembly["loads"]

        U0 = zeros((fulldofs), dtype=float64)           # empty((fulldofs, 1))
        U1 = zeros((fulldofs), dtype=float64)           # empty((fulldofs, 1))
        U = zeros((fulldofs, nsteps), dtype=float64)    # empty((fulldofs, nsteps))
        Uc = assembly["bcdirnh"]

        freedof = constrainsdof["freedof"]
        constdof = constrainsdof["constdof"]

        for step in range(nsteps):
            forcelist[freedof, step] = forcelist[freedof, step] - dot(
                stiffness[:, constdof][freedof, :].toarray(), Uc[constdof, step]
            )
            try:
                U1[freedof] = spsolve(
                    stiffness[:, freedof][freedof, :], forcelist[freedof, step]
                )
            except:
                pass
            U1[constdof] = Uc[constdof, step]
            U1[:] += U0[:]
            U[:, step] = U1
            U0[:] = U1[:]
        solution["U"] = U
        return solution
