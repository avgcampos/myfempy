from __future__ import annotations


from numpy import (array, arange, concatenate, empty, float64, newaxis, pi, sqrt,
                   unique, zeros)
from scipy.sparse.linalg import eigsh

from myfempy.core.solver.assemblerfull import AssemblerFULL
from myfempy.core.solver.assemblerfull_parallel import AssemblerFULLPOOL
from myfempy.core.solver.assemblersymm import AssemblerSYMM
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


class DynamicEigenLinear(Solver):
    """
    Dynamic Eigen (modal problem) Linear Solver Class <ConcreteClassService>
    """
    def getMatrixAssembler(Model, inci = None, coord = None, tabmat = None, tabgeo = None, intgauss = None, SYMM=None, MP=None):
        matrix = dict()
        if SYMM:
            matrix["stiffness"] = AssemblerSYMM.getLinearStiffnessGlobalMatrixAssembler(
                Model, inci, coord, tabmat, tabgeo, intgauss,
            )
            matrix["mass"] = AssemblerSYMM.getMassConsistentGlobalMatrixAssembler(
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
        return empty(
            (
                nodedof * nodetot,
                len(unique(loadaply[:, 3])),
            )
        )

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
        modeEnd = setSteps(solverset["STEPSET"])
        stiffness = assembly["stiffness"]
        mass = assembly["mass"]
        # forcelist = assembly["loads"]
        U = zeros((fulldofs, modeEnd), dtype=float64)
        freedof = constrainsdof["freedof"]
        try:
            W, U[freedof, :] = eigsh(
                A=stiffness[:, freedof][freedof, :],
                M=mass[:, freedof][freedof, :],
                k=modeEnd,
                sigma=1,
                which="LM",
                maxiter=1000,
            )
        except:
            pass
        Wlist = arange(0, modeEnd + 1)
        Wrad = sqrt(W)
        Whz = Wrad / (2 * pi)
        w_range = concatenate(
            (Wlist[1:, newaxis], Wrad[:, newaxis], Whz[:, newaxis]), axis=1
        )
        solution["U"] = U
        solution["FREQ"] = w_range
        return solution
