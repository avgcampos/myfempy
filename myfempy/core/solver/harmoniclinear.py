from __future__ import annotations

from os import environ
environ['OMP_NUM_THREADS'] = '3'

from numpy import zeros, linspace, pi, float64
from scipy.sparse.linalg import spsolve, minres 

# from myfempy.core.alglin import linsolve_spsolve
from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps
from myfempy.core.solver.assemblersymm import AssemblerSYMM
from myfempy.core.solver.assemblerfull import AssemblerFULL


class HarmonicLinear(Solver):
    
    """
    Harmonic Forced System Linear Solver Class <ConcreteClassService>
    """

    def getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, SYMM=None, MP=None):  
        matrix = dict()
        if SYMM:
            matrix['stiffness'] = AssemblerSYMM.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler = 'linear_stiffness', MP=MP)
            matrix['mass'] = AssemblerSYMM.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss,  type_assembler = 'mass_consistent', MP=MP)
        else:
            matrix['stiffness'] = AssemblerFULL.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler = 'linear_stiffness', MP=MP)
            matrix['mass'] = AssemblerFULL.getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss,  type_assembler = 'mass_consistent', MP=MP)
        return matrix  
        
    def getLoadAssembler(loadaply, nodetot, nodedof):
        return AssemblerSYMM.getLoadAssembler(loadaply, nodetot, nodedof)

    def getConstrains(constrains, nodetot, nodedof):
        return AssemblerSYMM.getConstrains(constrains, nodetot, nodedof)
        

    def runSolve(fulldofs, assembly, forcelist, freedof, solverset):
        solution = dict()
        stiffness = assembly['stiffness']
        mass = assembly['mass']
        twopi = 2 * pi
        freqStart = (twopi) * solverset["STEPSET"]["start"]
        freqEnd = (twopi) * solverset["STEPSET"]["end"]
        freqStep = setSteps(solverset["STEPSET"])
        w_range = linspace(freqStart, freqEnd, freqStep)
        U = zeros((fulldofs, freqStep), dtype=float64)
        U0 = U[freedof, 0]
        sA = stiffness[:, freedof][freedof, :]
        sM = mass[:, freedof][freedof, :]
        for ww in range(freqStep):
            Wn = w_range[ww]
            Dw = sA - (Wn**2) * sM
            # U[freedof, ww] = Solver.getLinSysSolve(Dw, forcelist[freedof, :])
            try:
                U[freedof, ww], info = minres(A=Dw, b=forcelist[freedof, :].toarray(), x0=U0, tol=1E-10, maxiter=1000)
            except:
                raise info    
        solution['U'] = U
        solution['FREQ'] = w_range / (twopi)
        return solution