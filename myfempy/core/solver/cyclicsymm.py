from __future__ import annotations

from numpy import (arange, array, concatenate, dot, float64, in1d, int16,
                   setdiff1d, where, zeros, sqrt)

from scipy.sparse import csc_matrix, lil_matrix, eye, hstack, vstack
from scipy.sparse.linalg import minres, spsolve

from myfempy.core.solver.assemblerfull import AssemblerFULL
from myfempy.core.solver.assemblerfull_parallel import AssemblerFULLPOOL
from myfempy.core.solver.assemblersymm import AssemblerSYMM
from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps


class StaticLinearCyclicSymmPlane(Solver):
    """
    Static Linear Cyclic Symmetry Plane Solver Class <ConcreteClassService>
    """

    def getMatrixAssembler(
        Model, SYMM=None, MP=None
    ):
        matrix = dict()

        if SYMM:
            matrix["stiffness"] = AssemblerSYMM.getLinearStiffnessGlobalMatrixAssembler(
                Model,
            )
        else:
            if MP:
                matrix["stiffness"] = AssemblerFULLPOOL.getLinearStiffnessGlobalMatrixAssembler(
                    Model,
                    MP=MP,
                )
            else:
                matrix["stiffness"] = AssemblerFULL.getLinearStiffnessGlobalMatrixAssembler(
                    Model,
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

    def runSolve(Model, Physic, assembly, constrainsdof, solverset):
        fulldofs = Model.modelinfo["fulldofs"]

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

        RM_left = StaticLinearCyclicSymmPlane.getRotationMatrix2D(Physic.csleft, Model.coord, leftdof.shape[0])
        RM_right = StaticLinearCyclicSymmPlane.getRotationMatrix2D(Physic.csright, Model.coord, rightdof.shape[0])

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
        Irr = eye(rightdof.shape[0], rightdof.shape[0], dtype=int16)
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

    
    # https://en.wikipedia.org/wiki/Rotation_matrix
    def getRotationMatrix2D(node_list, coord, ndof):
        # Initialize RM as a sparse matrix
        RM = lil_matrix((ndof, ndof))
        for n in range(node_list.shape[0]):
            nol = int(node_list[n] - 1)
            RonX = coord[nol, 1]  # - Og[0]
            RonY = coord[nol, 2]  # - Og[1]
            Ron = sqrt(RonX**2 + RonY**2)
            S_the = RonY / Ron
            C_the = RonX / Ron
            # Assign values to the sparse matrix
            RM[2 * n, 2 * n] = C_the
            RM[2 * n, 2 * n + 1] = -S_the
            RM[2 * n + 1, 2 * n] = S_the
            RM[2 * n + 1, 2 * n + 1] = C_the

        # Convert to CSR format for more efficient arithmetic and matrix-vector operations
        return RM.tocsr()

    # def getRotationMatrix3D(node_list, coord, ndof):
    #     """
    #     https://en.wikipedia.org/wiki/Rotation_matrix

    #     Arguments:
    #         node_list -- _description_
    #         coord -- _description_
    #         ndof -- _description_

    #     Returns:
    #         RM
    #     """
    #     # Inicialize a RM como uma matriz esparsa
    #     RM = lil_matrix((ndof, ndof))
    #     for n in range(node_list.shape[0]):
    #         nol = int(node_list[n] - 1)
    #         RonX = coord[nol, 1]
    #         RonY = coord[nol, 2]
    #         RonZ = coord[nol, 3]
    #         Ron = np.sqrt(RonX**2 + RonY**2 + RonZ**2)
    #         S_phi = RonY / Ron
    #         C_phi = RonX / Ron
    #         S_theta = RonZ / Ron
    #         C_theta = np.sqrt(RonX**2 + RonY**2) / Ron

    #         # Atribuir valores à matriz esparsa
    #         RM[3 * n, 3 * n] = C_theta
    #         RM[3 * n, 3 * n + 1] = -S_phi * S_theta
    #         RM[3 * n, 3 * n + 2] = C_phi * S_theta
    #         RM[3 * n + 1, 3 * n] = S_phi
    #         RM[3 * n + 1, 3 * n + 1] = C_phi * C_theta
    #         RM[3 * n + 1, 3 * n + 2] = -S_theta
    #         RM[3 * n + 2, 3 * n] = -S_theta
    #         RM[3 * n + 2, 3 * n + 1] = C_theta
    #         RM[3 * n + 2, 3 * n + 2] = C_theta * C_phi

    #     # Converter para o formato CSR para operações aritméticas e de matriz-vetor mais eficientes
    #     return RM.tocsr() 