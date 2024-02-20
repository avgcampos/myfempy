from __future__ import annotations

from numpy import zeros, float64
from scipy import empty, sparse

# from myfempy.core.alglin import linsolve_spsolve
from myfempy.core.solver.solver import Solver

from myfempy.experimental.asmb_cython.import_assembler_cython2py import getMatrixAssemblerSYMM

class StaticLinear(Solver):

    """
     Static Linear Solver Class <ConcreteClassService>
    """
    
    def getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss):  
        matrix = dict()
        # matrix['stiffness'] = Solver.getMatrixAssemblerFULL(Model, inci, coord, tabmat, tabgeo, intgauss,  type_assembler = 'linear_stiffness')
        # matrix['stiffness'] = Solver.getMatrixAssemblerSYMM(Model, inci, coord, tabmat, tabgeo, intgauss,  type_assembler = 'linear_stiffness')
        matrix['stiffness'] = getMatrixAssemblerSYMM(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler = 'linear_stiffness')
        return matrix

    def getLoadAssembler(loadaply, nodetot, nodedof):
        return Solver.getLoadAssembler(loadaply, nodetot, nodedof)
    
    def addMatrix(A, A_add, loc):
        return Solver.addMatrix(A, A_add, loc)
    
    def getConstrains(constrains, nodetot, nodedof):
        return Solver.getConstrains(constrains, nodetot, nodedof)
    
    def setSteps(steps):
        return Solver.setSteps(steps)
        

    def runSolve(fulldofs, assembly, forcelist, freedof, solverset):    
        solution = dict()    
        stiffness = assembly['stiffness']
        nsteps = StaticLinear.setSteps(solverset["STEPSET"])
        U0 = zeros((fulldofs, 1))     #empty((fulldofs, 1))
        U1 = zeros((fulldofs, 1))     #empty((fulldofs, 1))
        U = zeros((fulldofs, nsteps)) #empty((fulldofs, nsteps))
        
        # U0 = sparse.csc_matrix((fulldofs, 1), dtype=float64)
        # U1 = sparse.csc_matrix((fulldofs, 1), dtype=float64)   
        # U = sparse.csc_matrix((fulldofs, nsteps), dtype=float64)
                        
        for step in range(nsteps):
            U1[freedof, 0] = Solver.getLinSysSolve(stiffness[:, freedof][freedof, :], forcelist[freedof, step])                    
            U1[freedof, 0] += U0[freedof, 0]
            U[freedof, step] = U1[freedof, 0]
            U0[freedof, 0] = U1[freedof, 0]        
        solution['U'] = U
        return solution