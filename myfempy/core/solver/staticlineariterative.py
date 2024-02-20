from __future__ import annotations

import numpy as np
import scipy as sp

# from myfempy.core.alglin import linsolve_gmres
from myfempy.core.solver import Solver


class StaticLinearIterative(Solver):

    """
     Static Linear Solver Iterative Class <ConcreteClassService>
    """
    
    def getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss):  
        matrix = dict()
        matrix['stiffness'] = Solver.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss,  type_assembler = 'linear_stiffness')
        return matrix

    def getLoadAssembler(loadaply, nodetot, nodedof):
        return Solver.getLoadAssembler(loadaply, nodetot, nodedof)
    
    def addMatrix(A, A_add, loc):
        return Solver.addMatrix(A, A_add, loc)
    
    def getConstrains(constrains, nodetot, nodedof):
        return Solver.getConstrains(constrains, nodetot, nodedof)
    
    def setSteps(steps):
        return Solver.setSteps(steps)
        

    def Solve(fulldofs, assembly, forcelist, freedof, solverset):        
        stiffness = assembly['stiffness']
        nsteps = StaticLinearIterative.setSteps(solverset["STEPSET"])
        U0 = np.zeros((fulldofs, 1))     #empty((fulldofs, 1))
        U1 = np.zeros((fulldofs, 1))     #empty((fulldofs, 1))
        U = np.zeros((fulldofs, nsteps)) #empty((fulldofs, nsteps))
        Nshape = np.size(freedof)
        sA = stiffness[:, freedof][freedof, :]
        def A_ilu(x):
            return sp.sparse.linalg.spsolve(sA, x)
        M = sp.sparse.linalg.LinearOperator((Nshape, Nshape), A_ilu)
        solution = dict()
        for step in range(nsteps):
            try:
                U1[freedof, 0], info = Solver.getGenMinResSolve(sA, forcelist[freedof, step].toarray(), U0[freedof, 0], M)
            except:
                raise info          
            U1[freedof, 0] += U0[freedof, 0]
            U[freedof, step] = U1[freedof, 0]
            U0[freedof, 0] = U1[freedof, 0]
        solution['U'] = U
        solution['INFO'] = info
        return solution