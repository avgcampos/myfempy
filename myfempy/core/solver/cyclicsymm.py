from __future__ import annotations

from os import environ

environ["OMP_NUM_THREADS"] = "3"

from numpy import (arange, array, concatenate, dot, float64, in1d, int16,
                   setdiff1d, where, zeros)
from scipy.sparse import csc_matrix, eye, hstack, vstack
from scipy.sparse.linalg import minres, spsolve

from myfempy.core.solver.assemblerfull import AssemblerFULL
from myfempy.core.solver.assemblersymm import AssemblerSYMM
# from myfempy.core.alglin import linsolve_spsolve
from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps


class StaticLinearCyclicSymm(Solver):

    """
    Static Linear Cyclic Symmetry Solver Class <ConcreteClassService>
    """

    def getMatrixAssembler(
        Model, inci, coord, tabmat, tabgeo, intgauss, SYMM=None, MP=None
    ):
        matrix = dict()

        if SYMM:
            matrix["stiffness"] = AssemblerSYMM.getLinearStiffnessGlobalMatrixAssembler(
                Model,
                inci,
                coord,
                tabmat,
                tabgeo,
                intgauss,
                type_assembler="linear_stiffness",
                MP=MP,
            )
        else:
            matrix["stiffness"] = AssemblerFULL.getLinearStiffnessGlobalMatrixAssembler(
                Model,
                inci,
                coord,
                tabmat,
                tabgeo,
                intgauss,
                type_assembler="linear_stiffness",
                MP=MP,
            )
        return matrix

    def getLoadAssembler(loadaply, nodetot, nodedof):
        return AssemblerFULL.getLoadAssembler(loadaply, nodetot, nodedof)

    def getRotationMatrix(node_list, coord, ndof):
        return AssemblerFULL.getRotationMatrix(node_list, coord, ndof)

    def getConstrains(constrains, nodetot, nodedof):
        cs_left = where(constrains[:, 1] == 11)
        cs_right = where(constrains[:, 1] == 12)

        cs_left_constrain = constrains[cs_left[0], :]
        cs_right_constrain = constrains[cs_right[0], :]

        __, fixed_left_dof, __ = AssemblerFULL.getConstrains(
            cs_left_constrain, nodetot, nodedof
        )

        __, fixed_right_dof, __ = AssemblerFULL.getConstrains(
            cs_right_constrain, nodetot, nodedof
        )

        fixed_nodes = where(constrains[:, 1] == 0)
        fixed_constrain = constrains[fixed_nodes[0], :]

        freedof, fixedof_constrain, constdof = AssemblerFULL.getConstrains(
            fixed_constrain, nodetot, nodedof
        )

        # testl = in1d(fixedof_constrain, fixed_left_dof, assume_unique=True, invert = True)
        # testr = in1d(fixedof_constrain, fixed_right_dof, assume_unique=True, invert = True)

        # fixedof_constrain = fixedof_constrain[testl == testr]

        fixedof = [fixed_left_dof, fixed_right_dof, fixedof_constrain]

        return freedof, fixedof, constdof

    def getDirichletNH(constrains, nodetot, nodedof):
        return AssemblerFULL.getDirichletNH(constrains, nodetot, nodedof)

    def runSolve(assembly, constrainsdof, modelinfo, solverset):
        fulldofs = modelinfo["fulldofs"]

        solution = dict()
        nsteps = setSteps(solverset["STEPSET"])

        leftdof = constrainsdof["fixedof"][0]
        rightdof = constrainsdof["fixedof"][1]
        fixeddof = constrainsdof["fixedof"][2]
        freedof = constrainsdof["freedof"]
        constdof = constrainsdof["constdof"]

        testl = in1d(freedof, leftdof, assume_unique=True)
        testr = in1d(freedof, rightdof, assume_unique=True)

        interdof = freedof[testl == testr]

        stiffness = assembly["stiffness"]
        forcelist = assembly["loads"]

        RM_left = AssemblerFULL.getRotationMatrix(
            modelinfo["csleft"], modelinfo["coord"], leftdof.shape[0]
        )
        RM_right = AssemblerFULL.getRotationMatrix(
            modelinfo["csright"], modelinfo["coord"], rightdof.shape[0]
        )

        FG_cell = vstack(
            [forcelist[interdof, :], forcelist[leftdof, :], forcelist[rightdof, :]]
        )

        KG_cell = vstack(
            [
                hstack(
                    [
                        stiffness[interdof, :][:, interdof],
                        stiffness[interdof, :][:, leftdof],
                        stiffness[interdof, :][:, rightdof],
                    ]
                ),
                hstack(
                    [
                        stiffness[leftdof, :][:, interdof],
                        stiffness[leftdof, :][:, leftdof],
                        stiffness[leftdof, :][:, rightdof],
                    ]
                ),
                hstack(
                    [
                        stiffness[rightdof, :][:, interdof],
                        stiffness[rightdof, :][:, leftdof],
                        stiffness[rightdof, :][:, rightdof],
                    ]
                ),
            ]
        )

        Iii = eye(interdof.shape[0], interdof.shape[0], dtype=int16)
        Ill = eye(leftdof.shape[0], leftdof.shape[0], dtype=int16)
        Irr = eye(leftdof.shape[0], leftdof.shape[0], dtype=int16)
        Zil = csc_matrix((interdof.shape[0], leftdof.shape[0]), dtype=int16)
        Zll = csc_matrix((leftdof.shape[0], leftdof.shape[0]), dtype=int16)

        Rcs = vstack(
            [
                hstack([Iii, Zil, Zil]),
                hstack([Zil.transpose(), RM_left, Zll]),
                hstack([Zil.transpose(), Zll.transpose(), RM_right]),
            ]
        )

        KG_cell = dot(dot(Rcs.transpose(), KG_cell), Rcs)
        FG_cell = dot(Rcs.transpose(), FG_cell)

        #  CONDENSACAO ESTATICA
        MATCON = vstack(
            [
                hstack([Iii, Zil]),
                hstack([Zil.transpose(), Ill]),
                hstack([Zil.transpose(), Irr]),
            ]
        )

        KG_con_cs = dot(dot(MATCON.transpose(), KG_cell), MATCON)
        FG_con_cs = dot(MATCON.transpose(), FG_cell)

        fulldof_con_cs = concatenate((interdof, leftdof), axis=0)

        fixed_list_con_cs = setdiff1d(fixeddof, rightdof)

        fixed_list_full = setdiff1d(fixed_list_con_cs, leftdof)

        leftdof_con_cs = where(
            in1d(fulldof_con_cs, leftdof, assume_unique=True) == True
        )[0]
        interdof_con_cs = where(
            in1d(fulldof_con_cs, interdof, assume_unique=True) == True
        )[0]
        free_list_full = concatenate((interdof, leftdof, rightdof), axis=0)
        rightdof_con_cs = where(
            in1d(free_list_full, rightdof, assume_unique=True) == True
        )[0]

        fixedof_con_cs = where(
            in1d(fulldof_con_cs, fixed_list_con_cs, assume_unique=True) == True
        )[0]
        freedof_con_cs = where(
            in1d(fulldof_con_cs, fixed_list_con_cs, assume_unique=True) == False
        )[0]

        U0 = zeros((fulldof_con_cs.shape[0]), dtype=float64)  # empty((fulldofs, 1))
        U1 = zeros((fulldof_con_cs.shape[0]), dtype=float64)  # empty((fulldofs, 1))
        U = zeros((fulldofs, nsteps), dtype=float64)  # empty((fulldofs, nsteps))

        # Uc = assembly["bcdirnh"]
        for step in range(nsteps):
            # FG_con_cs[freedof_con_cs, step] = FG_con_cs[freedof_con_cs, step] - dot(
            #     KG_con_cs[:, constdof][freedof_con_cs, :].toarray(), Uc[constdof, step]
            # )
            try:
                U1[freedof_con_cs], info = minres(
                    A=KG_con_cs[:, freedof_con_cs][freedof_con_cs, :],
                    b=FG_con_cs[freedof_con_cs, step].toarray(),
                    tol=1e-10,
                    maxiter=1000,
                )
            except:
                raise info
            # U1[constdof] = Uc[constdof, step]

            U1[:] += U0[:]

            U_exp = dot(MATCON.toarray(), U1)

            Uxy = concatenate(
                (
                    U_exp[interdof_con_cs],
                    dot(RM_left.toarray(), U_exp[leftdof_con_cs]),
                    dot(RM_right.toarray(), U_exp[rightdof_con_cs]),
                ),
                axis=0,
            )

            U[free_list_full, step] = Uxy
            U0[:] = U1[:]
        solution["U"] = U
        return solution
