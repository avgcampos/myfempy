from __future__ import annotations

import time

import numpy as np
# import jax.numpy as np
import scipy.sparse.linalg as spla

from myfempy.core.alglin import linsolve_direct
from myfempy.core.solver import Solver


class HarmonicForced(Solver):
    '''Harmonic Forced System Linear Solver Class <ConcreteClassService>'''
    
    def getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss):
                
        matrix = dict()
        startstep = time.time()
        matrix['stiffness'] = Solver.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss,  type_assembler = 'linear_stiffness')
        matrix['mass'] = Solver.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss,  type_assembler = 'mass_consistent')
        endstep = time.time()
        print("\nGLOBAL ASSEMBLY TIME ", "\ TIME SPEND: ", endstep - startstep, " SEC")
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
        freqStep = HarmonicForced.setSteps(solverset["STEPSET"])
        w_range = np.linspace(freqStart, freqEnd, freqStep)

        U = np.zeros((fulldofs, freqStep))
        startstep = time.time()
        for ww in range(freqStep):
            Wn = w_range[ww]
            Dw = (stiffness[:, freedof][freedof, :]) - (Wn**2) * (mass[:, freedof][freedof, :])
            U[freedof, ww] = linsolve_direct(Dw, forcelist[freedof, :])
            # U[freedof, ww], info = spla.bicgstab(A=Dw, b=forcelist[freedof, :].toarray(), tol=solverset['TOL'])
            # if info > 0:
            #     pass
            # elif info < 0:
            #     print('ILLEGAL INPUT OR BREAKDOWN')
            # else:
            #     pass
        endstep = time.time()
        print("\nSTEP --: SUCCESSFUL CONVERGED\n")
        print("\ TIME SPEND: ", endstep - startstep, " SEC")
                
        solution['U'] = U
        solution['FREQ'] = w_range / (twopi)
        
        return solution