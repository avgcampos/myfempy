from __future__ import annotations

from os import environ
environ['OMP_NUM_THREADS'] = '3'

from numpy import zeros, empty, unique, arange, sqrt, concatenate, pi, newaxis, float64
from scipy.sparse.linalg import eigsh

from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps
from myfempy.expe.asmb_cython.import_assembler_cython2py import getMatrixAssemblerSYMM
from myfempy.core.solver.assemblersymm import AssemblerSYMM

class ModalLinear(Solver):
    
    """
     Modal(eig problem) Linear Solver Class <ConcreteClassService>
    """
    
    def getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss):  
        matrix = dict()
        matrix['stiffness'] = AssemblerSYMM.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler = 'linear_stiffness')
        # matrix['stiffness'] = getMatrixAssemblerSYMM(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler = 'linear_stiffness')
        matrix['mass'] = AssemblerSYMM.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss,  type_assembler = 'mass_consistent')
        return matrix 
    
    def getLoadAssembler(loadaply, nodetot, nodedof):
        return empty((nodedof *nodetot, len(unique(loadaply[:, 3])),))
          
    def getConstrains(constrains, nodetot, nodedof):
        return AssemblerSYMM.getConstrains(constrains, nodetot, nodedof)
   
    def runSolve(fulldofs, assembly, forcelist, freedof, solverset):      
        solution = dict()
        stiffness = assembly['stiffness']
        mass = assembly['mass']
        modeEnd = setSteps(solverset["STEPSET"])
        U = zeros((fulldofs, modeEnd), dtype=float64)
        # W, U[freedof, :] = Solver.getEigHerSysSolve(stiffness[:, freedof][freedof, :], mass[:, freedof][freedof, :], modeEnd)
        W, U[freedof, :] = eigsh(A=stiffness[:, freedof][freedof, :], M=mass[:, freedof][freedof, :], k=modeEnd, sigma=1, which="LM", maxiter=1000)
        Wlist = arange(0, modeEnd + 1)
        Wrad = sqrt(W)
        Whz = Wrad / (2 * pi)
        w_range = concatenate((Wlist[1:, newaxis], Wrad[:, newaxis], Whz[:, newaxis]), axis=1)
        solution['U'] = U
        solution['FREQ'] = w_range
        return solution