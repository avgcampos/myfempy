from __future__ import annotations

import time

import numpy as np
# import jax.numpy as np
import scipy.sparse.linalg as spla

from myfempy.core.alglin import eigsolve_eigsh
from myfempy.core.solver import Solver


class ModalLinear(Solver):
    '''Modal(eig problem) Linear Solver Class <ConcreteClassService>'''
    
    def getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss):
                
        matrix = dict()
        startstep = time.time()
        matrix['stiffness'] = Solver.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss,  type_assembler = 'linear_stiffness')
        matrix['mass'] = Solver.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss,  type_assembler = 'mass_consistent')
        endstep = time.time()
        print("\nGLOBAL ASSEMBLY TIME ", "\ TIME SPEND: ", endstep - startstep, " SEC")
        return matrix
    
    
    def getLoadAssembler(loadaply, nodetot, nodedof):
        forcevec = np.zeros((nodedof *nodetot, len(np.unique(loadaply[:, 3])),))
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
        startstep = time.time()
        W, U[freedof, :] = eigsolve_eigsh(stiffness[:, freedof][freedof, :], mass[:, freedof][freedof, :], modeEnd)
        endstep = time.time()
        print("\nSTEP --: SUCCESSFUL CONVERGED\n")
        print("\ TIME SPEND: ", endstep - startstep, " SEC")
        
        Wlist = np.arange(0, modeEnd + 1)
        Wrad = np.sqrt(W)
        Whz = Wrad / (2 * np.pi)
        
        w_range = np.concatenate(
            (Wlist[1:, np.newaxis], Wrad[:, np.newaxis], Whz[:, np.newaxis]), axis=1
        )
        
        solution['U'] = U
        solution['FREQ'] = w_range
        
        return solution