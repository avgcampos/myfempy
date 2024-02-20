from __future__ import annotations

import numpy as np
from scipy import empty

# from myfempy.core.alglin import linsolve_spsolve
from myfempy.core.solver.solver import Solver


class HarmonicLinear(Solver):
    
    """
    Harmonic Forced System Linear Solver Class <ConcreteClassService>
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

    def Solve(fulldofs, assembly, forcelist, freedof, solverset):
        solution = dict()
        stiffness = assembly['stiffness']
        mass = assembly['mass']
        twopi = 2 * np.pi
        freqStart = (twopi) * solverset["STEPSET"]["start"]
        freqEnd = (twopi) * solverset["STEPSET"]["end"]
        freqStep = HarmonicLinear.setSteps(solverset["STEPSET"])
        w_range = np.linspace(freqStart, freqEnd, freqStep)
        U = np.zeros((fulldofs, freqStep))
        for ww in range(freqStep):
            Wn = w_range[ww]
            Dw = (stiffness[:, freedof][freedof, :]) - (Wn**2) * (mass[:, freedof][freedof, :])
            U[freedof, ww] = Solver.getLinSysSolve(Dw, forcelist[freedof, :])
        solution['U'] = U
        solution['FREQ'] = w_range / (twopi)
        return solution