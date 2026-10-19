from __future__ import annotations

from numpy import (concatenate, flatnonzero, isin, where, zeros, dot)
from scipy.sparse import block_diag, identity, vstack
from scipy.sparse.linalg import minres

from myfempy.core.solver.assemblerSR import AssemblerSR
from myfempy.core.solver.assemblerMP import AssemblerMP
from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps

__docformat__ = "google"

__doc__ = """
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


class StaticLinearCyclicSymmPlane(Solver):
    """
    Static Linear Cyclic Symmetry Plane Solver Class <ConcreteClassService>
    """

    def getMatrixAssembler(
        Model,
        inci=None,
        coord=None,
        tabmat=None,
        tabgeo=None,
        intgauss=None,
        MP=None,
        max_workers=None,
    ):
        matrix = dict()
        if MP:
            matrix["stiffness"] = AssemblerMP.getGlobalMatrixAssembler(
                Model,
                Model.element.getStifLinearMat,
                Model.shape.getIntNumK,
                inci,
                coord,
                tabmat,
                tabgeo,
                intgauss,
                max_workers,
            )
        else:
            matrix["stiffness"] = AssemblerSR.getGlobalMatrixAssembler(
                Model,
                Model.element.getStifLinearMat,
                Model.shape.getIntNumK,
                inci,
                coord,
                tabmat,
                tabgeo,
                intgauss,
            )
        return matrix

    def getLoadAssembler(loadaply, nodetot, nodedof):
        return AssemblerSR.getLoadAssembler(loadaply, nodetot, nodedof)

    def getRotationMatrix(node_list, coord, ndof):
        return AssemblerSR.getRotationMatrix(node_list, coord, ndof)

    def getConstrains(constrains, nodetot, nodedof):
        cs_left = where(constrains[:, 1] == 11)
        cs_right = where(constrains[:, 1] == 12)

        cs_left_constrain = constrains[cs_left[0], :]
        cs_right_constrain = constrains[cs_right[0], :]

        __, fixed_left_dof, __ = AssemblerSR.getConstrains(
            cs_left_constrain, nodetot, nodedof
        )

        __, fixed_right_dof, __ = AssemblerSR.getConstrains(
            cs_right_constrain, nodetot, nodedof
        )

        fixed_nodes = where(constrains[:, 1] == 0)
        fixed_constrain = constrains[fixed_nodes[0], :]

        freedof, fixedof_constrain, constdof = AssemblerSR.getConstrains(
            fixed_constrain, nodetot, nodedof
        )

        testl = isin(freedof, fixed_left_dof, assume_unique=True)
        testr = isin(freedof, fixed_right_dof, assume_unique=True)

        freedof = freedof[testl == testr]

        fixedof = [fixed_left_dof, fixed_right_dof, fixedof_constrain]

        return freedof, fixedof, constdof

    def getDirichletNH(constrains, nodetot, nodedof):
        return AssemblerSR.getDirichletNH(constrains, nodetot, nodedof)

    def runSolve(Model, Physic, assembly, constrainsdof, solverset):
        nsteps = setSteps(solverset["STEPSET"])
        leftdof, rightdof, fixeddof = constrainsdof["fixedof"]
        interdof = constrainsdof["freedof"]
        ni, nl = len(interdof), len(leftdof)

        solution = dict()

        # ordem dos GDL: internos | esquerda | direita
        order = concatenate((interdof, leftdof, rightdof))

        # T = rotação (cs local -> global) @ condensação periódica (direita = esquerda)
        RM_left = AssemblerSR.getRotationMatrix(Physic.csleft, Model.coord, nl)
        RM_right = AssemblerSR.getRotationMatrix(Physic.csright, Model.coord, len(rightdof))

        rot = block_diag((identity(ni), RM_left, RM_right), format="csc")
        I = identity(ni + nl, format="csc")
        Rcs = rot @ vstack([I, I[ni:]], format="csc")

        # sistema reduzido, apenas GDL livres
        freedof = flatnonzero(~isin(order[: ni + nl], fixeddof))
        K = (Rcs.T @ assembly["stiffness"][order, :][:, order] @ Rcs).tocsr()
        F = Rcs.T @ assembly["loads"][order, :]  # esparsa @ densa -> ndarray
        Kff = K[freedof][:, freedof]
        Fff = F[freedof]

        u = zeros(ni + nl)
        U = zeros((Model.modelinfo["fulldofs"], nsteps))
        for step in range(nsteps):
            du, info = minres(Kff, Fff[:, step], rtol=1e-10, maxiter=1000)
            if info != 0:
                raise RuntimeError(f"minres não convergiu no passo {step} (info={info})")
            u[freedof] += du  # deslocamento acumulado entre passos
            U[order, step] = Rcs @ u

        solution["U"] = U
        return solution