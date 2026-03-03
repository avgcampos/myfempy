from __future__ import annotations


from numpy import dot, float64, zeros
from scipy.sparse.linalg import spsolve
import scipy.sparse as sp

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


class SteadyStateLinear(Solver):
    """
    Steady State Linear Solver Class <ConcreteClassService>
    """
    def getMatrixAssembler(Model, inci = None, coord = None, tabmat = None, tabgeo = None, intgauss = None, SYMM=None, MP=None):
       
        matrix = dict()
        
        if SYMM:
            assembler = AssemblerSYMM.getLinearStiffnessGlobalMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss)
        else:
            if MP:
                assembler = AssemblerFULLPOOL.getLinearStiffnessGlobalMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, MP)
            else:
                assembler = AssemblerFULL.getLinearStiffnessGlobalMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss)
        
        matrix["stiffness"] = assembler #sp.load_npz('sparse_matrix.npz')
        
        return matrix

    def getLoadAssembler(loadaply, nodetot, nodedof):
        return AssemblerFULL.getLoadAssembler(loadaply, nodetot, nodedof)

    def getConstrains(constrains, nodetot, nodedof):
        return AssemblerFULL.getConstrains(constrains, nodetot, nodedof)

    def getDirichletNH(constrains, nodetot, nodedof):
        return AssemblerFULL.getDirichletNH(constrains, nodetot, nodedof)

    def runSolve(Model, Physic, assembly, constrainsdof, solverset):
        fulldofs = Model.modelinfo["fulldofs"]

        solution = dict()
        nsteps = setSteps(solverset["STEPSET"])

        stiffness = assembly["stiffness"]
        forcelist = assembly["loads"]

        U0 = zeros((fulldofs), dtype=float64)          
        U1 = zeros((fulldofs), dtype=float64)           
        U = zeros((fulldofs, nsteps), dtype=float64)
        Uc = assembly["bcdirnh"]

        freedof = constrainsdof["freedof"]
        constdof = constrainsdof["constdof"]

        for step in range(nsteps):
            forcelist[freedof, step] = forcelist[freedof, step] - dot(
                stiffness[:, constdof][freedof, :].toarray(), Uc[constdof, step]
            )
            try:
                U1[freedof] = spsolve(
                    stiffness[:, freedof][freedof, :], forcelist[freedof, step]
                )
            except:
                pass
            U1[constdof] = Uc[constdof, step]
            U1[:] += U0[:]
            U[:, step] = U1
            U0[:] = U1[:]
        solution["U"] = U
        return solution
