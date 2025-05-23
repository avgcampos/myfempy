from __future__ import annotations

from numpy import (arange, array, empty, concatenate, unique, dot, float64, in1d, int16, int64,
                   sort, pi, where, zeros, sqrt, real, linspace, ceil, exp)

from scipy.sparse import csc_matrix, lil_matrix, eye, hstack, vstack
from scipy.sparse.linalg import minres, spsolve

from myfempy.core.solver.assemblerfull import AssemblerFULL
from myfempy.core.solver.assemblerfull_parallel import AssemblerFULLPOOL
from myfempy.core.solver.assemblersymm import AssemblerSYMM
from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps


class HomogenizationPlane(Solver):
    """
    Homogenization Plane Solver Class <ConcreteClassService>
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

    def getConstrains(constrains, nodetot, nodedof):
        pc_left = where(constrains[:, 1] == 13)
        pc_right = where(constrains[:, 1] == 14)
        pc_bottom = where(constrains[:, 1] == 15)
        pc_top = where(constrains[:, 1] == 16)
        pc_bottom_left = where(constrains[:, 1] == 17)
        pc_bottom_right = where(constrains[:, 1] == 18)
        pc_top_left = where(constrains[:, 1] == 19)
        pc_top_right = where(constrains[:, 1] == 20)

        pc_left_constrain = constrains[pc_left[0], :]
        pc_right_constrain = constrains[pc_right[0], :]
        pc_bottom_constrain = constrains[pc_bottom[0], :]
        pc_top_constrain = constrains[pc_top[0], :]
        pc_bottom_left_constrain = constrains[pc_bottom_left[0], :]
        pc_bottom_right_constrain = constrains[pc_bottom_right[0], :]
        pc_top_left_constrain = constrains[pc_top_left[0], :]
        pc_top_right_constrain = constrains[pc_top_right[0], :]

        testl2bl = in1d(pc_left_constrain[:,0], pc_bottom_left_constrain[:,0], assume_unique=True, invert=True)
        pc_left_constrain = pc_left_constrain[testl2bl,:]
        
        testl2tl = in1d(pc_left_constrain[:,0], pc_top_left_constrain[:,0], assume_unique=True, invert=True)
        pc_left_constrain = pc_left_constrain[testl2tl,:]
                        
        testb2bl = in1d(pc_bottom_constrain[:,0], pc_bottom_left_constrain[:,0], assume_unique=True, invert=True)
        pc_bottom_constrain = pc_bottom_constrain[testb2bl,:]

        testb2br = in1d(pc_bottom_constrain[:,0], pc_bottom_right_constrain[:,0], assume_unique=True, invert=True)
        pc_bottom_constrain = pc_bottom_constrain[testb2br,:]
        
        testt2tl = in1d(pc_top_constrain[:,0], pc_top_left_constrain[:,0], assume_unique=True, invert=True)
        pc_top_constrain = pc_top_constrain[testt2tl,:]
        
        testt2tr = in1d(pc_top_constrain[:,0], pc_top_right_constrain[:,0], assume_unique=True, invert=True)
        pc_top_constrain = pc_top_constrain[testt2tr,:]

        testr2tr = in1d(pc_right_constrain[:,0], pc_top_right_constrain[:,0], assume_unique=True, invert=True)
        pc_right_constrain = pc_right_constrain[testr2tr,:]

        testr2br = in1d(pc_right_constrain[:,0], pc_bottom_right_constrain[:,0], assume_unique=True, invert=True)
        pc_right_constrain = pc_right_constrain[testr2br,:]

        __, pc_left_dof, __ = AssemblerFULL.getConstrains(pc_left_constrain, nodetot, nodedof)
        __, pc_right_dof, __ = AssemblerFULL.getConstrains(pc_right_constrain, nodetot, nodedof)
        __, pc_bottom_dof, __ = AssemblerFULL.getConstrains(pc_bottom_constrain, nodetot, nodedof)
        __, pc_top_dof, __ = AssemblerFULL.getConstrains(pc_top_constrain, nodetot, nodedof)
        __, pc_bottom_left_dof, __ = AssemblerFULL.getConstrains(pc_bottom_left_constrain, nodetot, nodedof)
        __, pc_bottom_right_dof, __ = AssemblerFULL.getConstrains(pc_bottom_right_constrain, nodetot, nodedof)
        __, pc_top_left_dof, __ = AssemblerFULL.getConstrains(pc_top_left_constrain, nodetot, nodedof)
        __, pc_top_right_dof, __ = AssemblerFULL.getConstrains(pc_top_right_constrain, nodetot, nodedof)

        full_dofs = arange(0, nodedof * nodetot, 1, int)

        pc_dofs = concatenate((pc_left_dof, pc_bottom_dof, pc_bottom_left_dof, pc_right_dof, pc_top_dof, pc_bottom_right_dof, pc_top_left_dof, pc_top_right_dof), axis=0)
        testpc2i = in1d(full_dofs, pc_dofs, assume_unique=True, invert=True)
        freedof = full_dofs[testpc2i]

        constdof = [pc_left_dof, pc_bottom_dof, pc_bottom_left_dof, pc_right_dof, pc_top_dof, pc_bottom_right_dof, pc_top_left_dof, pc_top_right_dof]
        fixedof = [pc_bottom_left_dof]

        return freedof, fixedof, constdof

    def getDirichletNH(constrains, nodetot, nodedof):
        return AssemblerFULL.getDirichletNH(constrains, nodetot, nodedof)

    def runSolve(Model, Physic, assembly, constrainsdof, solverset):
        fulldofs = Model.modelinfo["fulldofs"]

        solution = dict()
        nsteps = setSteps(solverset["STEPSET"])
        
        fixeddof = constrainsdof["fixedof"]
        idof = constrainsdof["freedof"]
        ldof = constrainsdof["constdof"][0]
        bdof = constrainsdof["constdof"][1]
        bldof = constrainsdof["constdof"][2]
        rdof = constrainsdof["constdof"][3]
        tdof = constrainsdof["constdof"][4]
        brdof = constrainsdof["constdof"][5]
        tldof = constrainsdof["constdof"][6]
        trdof = constrainsdof["constdof"][7]

        # fulldof_pc = concatenate((idof, ldof, bdof, bldof, rdof, tdof, brdof, tldof, trdof), axis=0)
        fulldof_pc = concatenate((idof, ldof, bdof, bldof), axis=0)
        freedof_pc = where(in1d(fulldof_pc, fixeddof, assume_unique=True) == False)[0]
        
        Kg_fem = assembly["stiffness"]
        Fg_fem = assembly["loads"]

        # CONDENSACAO PERIODICA INF
        KG_cell = vstack([
                        hstack([Kg_fem[idof, :][:, idof],   Kg_fem[idof, :][:, ldof],   Kg_fem[idof, :][:, bdof],   Kg_fem[idof, :][:, bldof],  Kg_fem[idof, :][:, rdof],   Kg_fem[idof, :][:, tdof],   Kg_fem[idof, :][:, brdof],  Kg_fem[idof, :][:, tldof],  Kg_fem[idof, :][:, trdof]]),
                        hstack([Kg_fem[ldof, :][:, idof],   Kg_fem[ldof, :][:, ldof],   Kg_fem[ldof, :][:, bdof],   Kg_fem[ldof, :][:, bldof],  Kg_fem[ldof, :][:, rdof],   Kg_fem[ldof, :][:, tdof],   Kg_fem[ldof, :][:, brdof],  Kg_fem[ldof, :][:, tldof],  Kg_fem[ldof, :][:, trdof]]),
                        hstack([Kg_fem[bdof, :][:, idof],   Kg_fem[bdof, :][:, ldof],   Kg_fem[bdof, :][:, bdof],   Kg_fem[bdof, :][:, bldof],  Kg_fem[bdof, :][:, rdof],   Kg_fem[bdof, :][:, tdof],   Kg_fem[bdof, :][:, brdof],  Kg_fem[bdof, :][:, tldof],  Kg_fem[bdof, :][:, trdof]]),
                        hstack([Kg_fem[bldof, :][:, idof],  Kg_fem[bldof, :][:, ldof],  Kg_fem[bldof, :][:, bdof],  Kg_fem[bldof, :][:, bldof], Kg_fem[bldof, :][:, rdof],  Kg_fem[bldof, :][:, tdof],  Kg_fem[bldof, :][:, brdof], Kg_fem[bldof, :][:, tldof], Kg_fem[bldof, :][:, trdof]]),
                        hstack([Kg_fem[rdof, :][:, idof],   Kg_fem[rdof, :][:, ldof],   Kg_fem[rdof, :][:, bdof],   Kg_fem[rdof, :][:, bldof],  Kg_fem[rdof, :][:, rdof],   Kg_fem[rdof, :][:, tdof],   Kg_fem[rdof, :][:, brdof],  Kg_fem[rdof, :][:, tldof],  Kg_fem[rdof, :][:, trdof]]),
                        hstack([Kg_fem[tdof, :][:, idof],   Kg_fem[tdof, :][:, ldof],   Kg_fem[tdof, :][:, bdof],   Kg_fem[tdof, :][:, bldof],  Kg_fem[tdof, :][:, rdof],   Kg_fem[tdof, :][:, tdof],   Kg_fem[tdof, :][:, brdof],  Kg_fem[tdof, :][:, tldof],  Kg_fem[tdof, :][:, trdof]]),
                        hstack([Kg_fem[brdof, :][:, idof],  Kg_fem[brdof, :][:, ldof],  Kg_fem[brdof, :][:, bdof],  Kg_fem[brdof, :][:, bldof], Kg_fem[brdof, :][:, rdof],  Kg_fem[brdof, :][:, tdof],  Kg_fem[brdof, :][:, brdof], Kg_fem[brdof, :][:, tldof], Kg_fem[brdof, :][:, trdof]]),
                        hstack([Kg_fem[tldof, :][:, idof],  Kg_fem[tldof, :][:, ldof],  Kg_fem[tldof, :][:, bdof],  Kg_fem[tldof, :][:, bldof], Kg_fem[tldof, :][:, rdof],  Kg_fem[tldof, :][:, tdof],  Kg_fem[tldof, :][:, brdof], Kg_fem[tldof, :][:, tldof], Kg_fem[tldof, :][:, trdof]]),
                        hstack([Kg_fem[trdof, :][:, idof],  Kg_fem[trdof, :][:, ldof],  Kg_fem[trdof, :][:, bdof],  Kg_fem[trdof, :][:, bldof], Kg_fem[trdof, :][:, rdof],  Kg_fem[trdof, :][:, tdof],  Kg_fem[trdof, :][:, brdof], Kg_fem[trdof, :][:, tldof], Kg_fem[trdof, :][:, trdof]]),
                        ])
        
        FG_cell = vstack([Fg_fem[idof, :], Fg_fem[ldof, :], Fg_fem[bdof, :], Fg_fem[bldof, :], Fg_fem[rdof, :], Fg_fem[tdof, :], Fg_fem[brdof, :], Fg_fem[tldof, :], Fg_fem[trdof, :]])

        Iii = eye(idof.shape[0], idof.shape[0], dtype=int64)
        Zil = csc_matrix((idof.shape[0], ldof.shape[0]), dtype=int64)
        Zib = csc_matrix((idof.shape[0], bdof.shape[0]), dtype=int64)
        Zibl = csc_matrix((idof.shape[0], bldof.shape[0]), dtype=int64)
        
        Zli = csc_matrix((ldof.shape[0], idof.shape[0]), dtype=int64)
        Ill = eye(ldof.shape[0], ldof.shape[0], dtype=int64)
        Zlb = csc_matrix((ldof.shape[0], bdof.shape[0]), dtype=int64)
        Zlbl = csc_matrix((ldof.shape[0], bldof.shape[0]), dtype=int64)

        Zbi = csc_matrix((bdof.shape[0], idof.shape[0]), dtype=int64)
        Zbl = csc_matrix((bdof.shape[0], ldof.shape[0]), dtype=int64)
        Ibb = eye(bdof.shape[0], bdof.shape[0], dtype=int64)
        Zbbl = csc_matrix((bdof.shape[0], bldof.shape[0]), dtype=int64)

        Zbli = csc_matrix((bldof.shape[0], idof.shape[0]), dtype=int64)
        Zbll = csc_matrix((bldof.shape[0], ldof.shape[0]), dtype=int64)
        Zblb = csc_matrix((bldof.shape[0], bdof.shape[0]), dtype=int64)
        Iblbl = eye(bldof.shape[0], bldof.shape[0], dtype=int64)

        Zri = csc_matrix((rdof.shape[0], idof.shape[0]), dtype=int64)
        Irl = eye(rdof.shape[0], ldof.shape[0], dtype=int64)
        Zrb = csc_matrix((rdof.shape[0], bdof.shape[0]), dtype=int64)
        Zrbl = csc_matrix((rdof.shape[0], bldof.shape[0]), dtype=int64)

        Zti = csc_matrix((tdof.shape[0], idof.shape[0]), dtype=int64)
        Ztl = csc_matrix((tdof.shape[0], ldof.shape[0]), dtype=int64)
        Itb = eye(tdof.shape[0], bdof.shape[0], dtype=int64)
        Ztbl = csc_matrix((tdof.shape[0], bldof.shape[0]), dtype=int64)

        Zbri = csc_matrix((brdof.shape[0], idof.shape[0]), dtype=int64)
        Zbrl = csc_matrix((brdof.shape[0], ldof.shape[0]), dtype=int64)
        Zbrb = csc_matrix((brdof.shape[0], bdof.shape[0]), dtype=int64)
        Ibrbl = eye(brdof.shape[0], bldof.shape[0], dtype=int64)

        Ztli = csc_matrix((tldof.shape[0], idof.shape[0]), dtype=int64)
        Ztll = csc_matrix((tldof.shape[0], ldof.shape[0]), dtype=int64)
        Ztlb = csc_matrix((tldof.shape[0], bdof.shape[0]), dtype=int64)
        Itlbl = eye(tldof.shape[0], bldof.shape[0], dtype=int64)

        Ztri = csc_matrix((trdof.shape[0], idof.shape[0]), dtype=int64)
        Ztrl = csc_matrix((trdof.shape[0], ldof.shape[0]), dtype=int64)
        Ztrb = csc_matrix((trdof.shape[0], bdof.shape[0]), dtype=int64)
        Itrbl = eye(trdof.shape[0], bldof.shape[0], dtype=int64)

        MAT_R = vstack([ 
            hstack([Iii, Zil, Zib, Zibl]),
            hstack([Zli, Ill, Zlb, Zlbl]),
            hstack([Zbi, Zbl, Ibb, Zbbl]),
            hstack([Zbli, Zbll, Zblb, Iblbl]),
            hstack([Zri, Irl, Zrb, Zrbl]),
            hstack([Zti, Ztl, Itb, Ztbl]),
            hstack([Zbri, Zbrl, Zbrb, Ibrbl]),
            hstack([Ztli, Ztll, Ztlb, Itlbl]),
            hstack([Ztri, Ztrl, Ztrb, Itrbl]),
        ])

        # Static Condensation !!!
        KG_cell_SC = (dot(dot(MAT_R.transpose(), KG_cell), MAT_R)).tocsc()
        FG_cell_SC = (dot(MAT_R.transpose(), FG_cell)).tocsc()

        U0 = zeros((fulldof_pc.shape[0], Fg_fem.shape[1]), dtype=float64)  # empty((fulldofs, nsteps))
        try:
            X = spsolve(A=KG_cell_SC[:, freedof_pc][freedof_pc, :], b=FG_cell_SC[freedof_pc, :])
        except:
            raise 'erro'
        U0[freedof_pc, :] = X.toarray()
        U_FULL_PC = dot(MAT_R.toarray(), U0)

        # # Uc = assembly["bcdirnh"]
        # for step in range(nsteps):
        #     # FG_con_cs[freedof_con_cs, step] = FG_con_cs[freedof_con_cs, step] - dot(
        #     #     KG_con_cs[:, constdof][freedof_con_cs, :].toarray(), Uc[constdof, step]
        #     # )
        #     try:
        #         U1[freedof_con_cs], info = minres(
        #             A=KG_con_cs[:, freedof_con_cs][freedof_con_cs, :],
        #             b=FG_con_cs[freedof_con_cs, step].toarray(),
        #             tol=1e-10,
        #             maxiter=1000,
        #         )
        #     except:
        #         raise info
        #     # U1[constdof] = Uc[constdof, step]

        #     U1[:] += U0[:]

        #     U_exp = dot(MATCON.toarray(), U1)

        #     Uxy = concatenate(
        #         (
        #             U_exp[interdof_con_cs],
        #             dot(RM_left.toarray(), U_exp[leftdof_con_cs]),
        #             dot(RM_right.toarray(), U_exp[rightdof_con_cs]),
        #         ),
        #         axis=0,
        #     )

        #     U[free_list_full, step] = Uxy
        #     U0[:] = U1[:]

        solution["U"] = U_FULL_PC
        
        return solution