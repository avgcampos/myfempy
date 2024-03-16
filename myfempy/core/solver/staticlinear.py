from __future__ import annotations

from os import environ
environ['OMP_NUM_THREADS'] = '3'

from numpy import zeros, float64
from scipy.sparse.linalg import spsolve 

# from myfempy.core.alglin import linsolve_spsolve
from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps
from myfempy.core.solver.assemblersymm import AssemblerSYMM
from myfempy.core.solver.assemblerfull import AssemblerFULL

class StaticLinear(Solver):

    """
     Static Linear Solver Class <ConcreteClassService>
    """
    
    def getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, SYMM=None, MP=None):  
        matrix = dict()
        if SYMM:
            matrix['stiffness'] = AssemblerSYMM.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler = 'linear_stiffness', MP=MP)
        else:
            matrix['stiffness'] = AssemblerFULL.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler = 'linear_stiffness', MP=MP)
        return matrix

    def getLoadAssembler(loadaply, nodetot, nodedof):
        return AssemblerSYMM.getLoadAssembler(loadaply, nodetot, nodedof)

    def getConstrains(constrains, nodetot, nodedof):
        return AssemblerSYMM.getConstrains(constrains, nodetot, nodedof)
       
    def runSolve(fulldofs, assembly, forcelist, freedof, solverset):    
        solution = dict()    
        stiffness = assembly['stiffness']
        nsteps = setSteps(solverset["STEPSET"])
        U0 = zeros((fulldofs, 1), dtype=float64)     #empty((fulldofs, 1))
        U1 = zeros((fulldofs, 1), dtype=float64)     #empty((fulldofs, 1))
        U = zeros((fulldofs, nsteps), dtype=float64) #empty((fulldofs, nsteps))
        sA = stiffness[:, freedof][freedof, :]             
        for step in range(nsteps):
            U1[freedof, 0] = spsolve(sA, forcelist[freedof, step])                  
            U1[freedof, 0] += U0[freedof, 0]
            U[freedof, step] = U1[freedof, 0]
            U0[freedof, 0] = U1[freedof, 0]        
        solution['U'] = U
        return solution