from __future__ import annotations

import numpy as np
from scipy import empty

from myfempy.core.solver.scipy_solve import eigsolve_eigsh
from myfempy.core.solver.solver import Solver


class HarmonicModal(Solver):
   
    """
    HarmonicModal Harmonic Forced System Modal Solver Class <ConcreteClassService>
    """

    def getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss):
        matrix = dict()
        matrix['stiffness'] = Solver.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss,  type_assembler = 'linear_stiffness')
        matrix['mass'] = Solver.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss,  type_assembler = 'mass_consistent')
        return matrix
    
    
    def getLoadAssembler(loadaply, nodetot, nodedof):
        return Solver.getLoadAssembler(loadaply, nodetot, nodedof)  
          
    def getConstrains(constrains, nodetot, nodedof):
        return Solver.getConstrains(constrains, nodetot, nodedof)
    
    def setSteps(steps):
        return Solver.setSteps(steps)      
    
    # def Solve(fulldofs, assembly, forcelist, freedof, solverset):