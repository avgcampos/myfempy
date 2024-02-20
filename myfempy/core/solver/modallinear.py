from __future__ import annotations

from scipy import empty
import numpy as np

# from myfempy.core.alglin import eigsolve_eigsh
from myfempy.core.solver.solver import Solver


class ModalLinear(Solver):
    
    """
     Modal(eig problem) Linear Solver Class <ConcreteClassService>
    """
    
    def getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss):
        matrix = dict()
        matrix['stiffness'] = Solver.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss,  type_assembler = 'linear_stiffness')
        matrix['mass'] = Solver.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss,  type_assembler = 'mass_consistent')
        return matrix
    
    
    def getLoadAssembler(loadaply, nodetot, nodedof):
        forcevec = empty((nodedof *nodetot, len(np.unique(loadaply[:, 3])),))
        return forcevec
          
    def getConstrains(constrains, nodetot, nodedof):
        return Solver.getConstrains(constrains, nodetot, nodedof)
    
    def setSteps(steps):
        return Solver.setSteps(steps)          

    def Solve(fulldofs, assembly, forcelist, freedof, solverset):      
        solution = dict()
        stiffness = assembly['stiffness']
        mass = assembly['mass']
        modeEnd = ModalLinear.setSteps(solverset["STEPSET"])
        U = np.zeros((fulldofs, modeEnd))
        W, U[freedof, :] = Solver.getEigHerSysSolve(stiffness[:, freedof][freedof, :], mass[:, freedof][freedof, :], modeEnd)
        Wlist = np.arange(0, modeEnd + 1)
        Wrad = np.sqrt(W)
        Whz = Wrad / (2 * np.pi)
        w_range = np.concatenate((Wlist[1:, np.newaxis], Wrad[:, np.newaxis], Whz[:, np.newaxis]), axis=1)
        solution['U'] = U
        solution['FREQ'] = w_range
        return solution