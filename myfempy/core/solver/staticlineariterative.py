from __future__ import annotations

from os import environ
environ['OMP_NUM_THREADS'] = '3'

from numpy import zeros, float64
from scipy.sparse.linalg import spsolve, minres

from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps
from myfempy.expe.asmb_cython.import_assembler_cython2py import getMatrixAssemblerSYMM
from myfempy.core.solver.assemblersymm import AssemblerSYMM

class StaticLinearIterative(Solver):

    """
     Static Linear Solver Iterative Class <ConcreteClassService>
    """
    
    def getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss):  
        matrix = dict()
        # matrix['stiffness'] = AssemblerSYMM.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler = 'linear_stiffness')
        matrix['stiffness'] = getMatrixAssemblerSYMM(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler = 'linear_stiffness')
        return matrix

    def getLoadAssembler(loadaply, nodetot, nodedof):
        return AssemblerSYMM.getLoadAssembler(loadaply, nodetot, nodedof)

    def getConstrains(constrains, nodetot, nodedof):
        return AssemblerSYMM.getConstrains(constrains, nodetot, nodedof)
        
    def runSolve(fulldofs, assembly, forcelist, freedof, solverset):
        solution = dict()
        
        nsteps = setSteps(solverset["STEPSET"])
        U0 = zeros((fulldofs, 1), dtype=float64)     #empty((fulldofs, 1))
        U1 = zeros((fulldofs, 1), dtype=float64)     #empty((fulldofs, 1))
        U = zeros((fulldofs, nsteps), dtype=float64) #empty((fulldofs, nsteps))
        stiffness = assembly['stiffness']
        sA = stiffness[:, freedof][freedof, :]
        sU0 = U0[freedof, 0]
        for step in range(nsteps):
            try:
                U1[freedof, 0], info = minres(A=sA, b=forcelist[freedof, step].toarray(), x0=sU0, tol=1E-10, maxiter=999)
            except:
                raise info          
            U1[freedof, 0] += U0[freedof, 0]
            U[freedof, step] = U1[freedof, 0]
            U0[freedof, 0] = U1[freedof, 0]
        solution['U'] = U
        solution['INFO'] = info
        return solution