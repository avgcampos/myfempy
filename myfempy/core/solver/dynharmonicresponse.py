from __future__ import annotations


from numpy import empty, float64, linspace, pi, unique, zeros
from scipy.sparse.linalg import minres, spsolve

from myfempy.core.solver.assemblerfull import AssemblerFULL
from myfempy.core.solver.assemblerfull_parallel import AssemblerFULLPOOL
from myfempy.core.solver.assemblersymm import AssemblerSYMM
# from myfempy.core.alglin import linsolve_spsolve
from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps


__docformat__ = "google"

__doc__ = """

==========================================================================
                            __                                
         _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
        | '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
        | | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
        |_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
                    |___/                       |_|     |___/ 
        myfempy -- MultiphYsics Finite Element Module to PYthon    
                    COMPUTATIONAL ANALYSIS PROGRAM                   
        Copyright (C) 2022-2026 Antonio Vinicius Garcia Campos        
==========================================================================
This Python file is part of myfempy project.

myfempy is a python package based on finite element method to multiphysics
analysis. The code is open source and *intended for educational and scientific
purposes only, not recommended to commercial use. The name myfempy is an acronym
for MultiphYsics Finite Elements Module to PYthon. You can help us by contributing
with the main project, send us a mensage on https://github.com/avgcampos/myfempy/discussions/10
If you use myfempy in your research, the  developers would be grateful if you 
could cite in your work.
																		
The code is written by Antonio Vinicius Garcia Campos.                                  
																		
A github repository, with the most up to date version of the code,      
can be found here: https://github.com/avgcampos/myfempy.                 
																		
The code is open source and intended for educational and scientific     
purposes only. If you use myfempy in your research, the developers      
would be grateful if you could cite this. The myfempy project is published
under the GPLv3, see the myfempy LICENSE on
https://github.com/avgcampos/myfempy/blob/main/LICENSE.
																		
Disclaimer:                                                             
The authors reserve all rights but do not guarantee that the code is    
free from errors. Furthermore, the authors shall not be liable in any   
event caused by the use of the program.

"""


class DynamicHarmonicResponseLinear(Solver):
    """
    Dynamic Harmonic Response Forced System Steady State Linear Solver Class <ConcreteClassService>
    """
    def getMatrixAssembler(Model, inci = None, coord = None, tabmat = None, tabgeo = None, intgauss = None, SYMM=None, MP=None):
        matrix = dict()
        if SYMM:
            matrix["stiffness"] = AssemblerSYMM.getLinearStiffnessGlobalMatrixAssembler(
                Model, inci, coord, tabmat, tabgeo, intgauss,
            )
            matrix["mass"] = AssemblerSYMM.getLinearStiffnessGlobalMatrixAssembler(
                Model, inci, coord, tabmat, tabgeo, intgauss,
            )
        else:
            if MP:
                matrix["stiffness"] = AssemblerFULLPOOL.getLinearStiffnessGlobalMatrixAssembler(
                    Model, inci, coord, tabmat, tabgeo, intgauss,
                    MP=MP,
                )
                matrix["mass"] = AssemblerFULLPOOL.getMassConsistentGlobalMatrixAssembler(
                    Model, inci, coord, tabmat, tabgeo, intgauss,
                    MP=MP,
                )
            else:
                matrix["stiffness"] = AssemblerFULL.getLinearStiffnessGlobalMatrixAssembler(
                    Model, inci, coord, tabmat, tabgeo, intgauss,
                )
                matrix["mass"] = AssemblerFULL.getMassConsistentGlobalMatrixAssembler(
                    Model, inci, coord, tabmat, tabgeo, intgauss,
                )
        return matrix

    def getLoadAssembler(loadaply, nodetot, nodedof):
        return AssemblerFULL.getLoadAssembler(loadaply, nodetot, nodedof)

    def getConstrains(constrains, nodetot, nodedof):
        return AssemblerFULL.getConstrains(constrains, nodetot, nodedof)

    def getDirichletNH(constrains, nodetot, nodedof):
        return empty(
            (nodedof * nodetot, len(unique(constrains[:, 3][constrains[:, 3] != 0]))),
            dtype=float64,
        )

    def runSolve(Model, Physic, assembly, constrainsdof, solverset):
        fulldofs = Model.modelinfo["fulldofs"]

        solution = dict()
        stiffness = assembly["stiffness"]
        mass = assembly["mass"]
        forcelist = assembly["loads"]

        freedof = constrainsdof["freedof"]

        twopi = 2 * pi
        freqStart = twopi * solverset["STEPSET"]["start"]
        freqEnd = twopi * solverset["STEPSET"]["end"]
        freqStep = setSteps(solverset["STEPSET"])
        w_range = linspace(freqStart, freqEnd, freqStep)

        U = zeros((fulldofs, freqStep), dtype=float64)
        U0 = U[freedof, 0]

        sA = stiffness[:, freedof][freedof, :]
        sM = mass[:, freedof][freedof, :]
        for ww in range(freqStep):
            Wn = w_range[ww]
            Dw = sA - (Wn**2) * sM
            try:
                U[freedof, ww], info = minres(
                    A=Dw, b=forcelist[freedof, 0], x0=U0, tol=1e-10, maxiter=1000
                )
            except:
                raise info
        solution["U"] = U
        solution["FREQ"] = w_range / (twopi)
        return solution
