from __future__ import annotations

from numpy import (arange, array, empty, concatenate, unique, dot, float64, in1d, int16, int64,
                   sort, pi, where, zeros, sqrt, real, linspace, ceil, exp)

from scipy.sparse import csc_matrix, lil_matrix, eye, hstack, vstack
from scipy.sparse.linalg import eigsh

from myfempy.core.solver.assemblerfull import AssemblerFULL
from myfempy.core.solver.assemblerfull_parallel import AssemblerFULLPOOL
from myfempy.core.solver.assemblersymm import AssemblerSYMM
from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps


class PhononicCrystalInPlane(Solver):
    """
    Phononic Crystal In-Plane Solver Class <ConcreteClassService>
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
            matrix["mass"] = AssemblerSYMM.getMassConsistentGlobalMatrixAssembler(
                Model,
                inci,
                coord,
                tabmat,
                tabgeo,
                intgauss,
                type_assembler="mass_consistent",
                MP=MP,
            )
        else:
            if MP:
                matrix["stiffness"] = AssemblerFULLPOOL.getLinearStiffnessGlobalMatrixAssembler(
                    Model,
                    inci,
                    coord,
                    tabmat,
                    tabgeo,
                    intgauss,
                    type_assembler="linear_stiffness",
                    MP=MP,
                )
                matrix["mass"] = AssemblerFULLPOOL.getMassConsistentGlobalMatrixAssembler(
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
                matrix["mass"] = AssemblerFULL.getMassConsistentGlobalMatrixAssembler(
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
        return empty(
            (
                nodedof * nodetot,
                len(unique(loadaply[:, 3])),
            )
        )
        

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
        fixedof = []

        return freedof, fixedof, constdof

    def getDirichletNH(constrains, nodetot, nodedof):
        return AssemblerFULL.getDirichletNH(constrains, nodetot, nodedof)

    def runSolve(assembly, constrainsdof, modelinfo, solverset):
        fulldofs = modelinfo["fulldofs"]

        solution = dict()
        modeEnd = setSteps(solverset["STEPSET"])
        cont_co = solverset["IBZ"]

        mu, tot_steps = PhononicCrystalInPlane.__setIBZ(cont_co)

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
        
        Kg_fem = assembly["stiffness"]
        Mg_fem = assembly["mass"]

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
        
        MG_cell =vstack([
                        hstack([Mg_fem[idof, :][:, idof],   Mg_fem[idof, :][:, ldof],   Mg_fem[idof, :][:, bdof],   Mg_fem[idof, :][:, bldof],  Mg_fem[idof, :][:, rdof],   Mg_fem[idof, :][:, tdof],   Mg_fem[idof, :][:, brdof],  Mg_fem[idof, :][:, tldof],  Mg_fem[idof, :][:, trdof]]),
                        hstack([Mg_fem[ldof, :][:, idof],   Mg_fem[ldof, :][:, ldof],   Mg_fem[ldof, :][:, bdof],   Mg_fem[ldof, :][:, bldof],  Mg_fem[ldof, :][:, rdof],   Mg_fem[ldof, :][:, tdof],   Mg_fem[ldof, :][:, brdof],  Mg_fem[ldof, :][:, tldof],  Mg_fem[ldof, :][:, trdof]]),
                        hstack([Mg_fem[bdof, :][:, idof],   Mg_fem[bdof, :][:, ldof],   Mg_fem[bdof, :][:, bdof],   Mg_fem[bdof, :][:, bldof],  Mg_fem[bdof, :][:, rdof],   Mg_fem[bdof, :][:, tdof],   Mg_fem[bdof, :][:, brdof],  Mg_fem[bdof, :][:, tldof],  Mg_fem[bdof, :][:, trdof]]),
                        hstack([Mg_fem[bldof, :][:, idof],  Mg_fem[bldof, :][:, ldof],  Mg_fem[bldof, :][:, bdof],  Mg_fem[bldof, :][:, bldof], Mg_fem[bldof, :][:, rdof],  Mg_fem[bldof, :][:, tdof],  Mg_fem[bldof, :][:, brdof], Mg_fem[bldof, :][:, tldof], Mg_fem[bldof, :][:, trdof]]),
                        hstack([Mg_fem[rdof, :][:, idof],   Mg_fem[rdof, :][:, ldof],   Mg_fem[rdof, :][:, bdof],   Mg_fem[rdof, :][:, bldof],  Mg_fem[rdof, :][:, rdof],   Mg_fem[rdof, :][:, tdof],   Mg_fem[rdof, :][:, brdof],  Mg_fem[rdof, :][:, tldof],  Mg_fem[rdof, :][:, trdof]]),
                        hstack([Mg_fem[tdof, :][:, idof],   Mg_fem[tdof, :][:, ldof],   Mg_fem[tdof, :][:, bdof],   Mg_fem[tdof, :][:, bldof],  Mg_fem[tdof, :][:, rdof],   Mg_fem[tdof, :][:, tdof],   Mg_fem[tdof, :][:, brdof],  Mg_fem[tdof, :][:, tldof],  Mg_fem[tdof, :][:, trdof]]),
                        hstack([Mg_fem[brdof, :][:, idof],  Mg_fem[brdof, :][:, ldof],  Mg_fem[brdof, :][:, bdof],  Mg_fem[brdof, :][:, bldof], Mg_fem[brdof, :][:, rdof],  Mg_fem[brdof, :][:, tdof],  Mg_fem[brdof, :][:, brdof], Mg_fem[brdof, :][:, tldof], Mg_fem[brdof, :][:, trdof]]),
                        hstack([Mg_fem[tldof, :][:, idof],  Mg_fem[tldof, :][:, ldof],  Mg_fem[tldof, :][:, bdof],  Mg_fem[tldof, :][:, bldof], Mg_fem[tldof, :][:, rdof],  Mg_fem[tldof, :][:, tdof],  Mg_fem[tldof, :][:, brdof], Mg_fem[tldof, :][:, tldof], Mg_fem[tldof, :][:, trdof]]),
                        hstack([Mg_fem[trdof, :][:, idof],  Mg_fem[trdof, :][:, ldof],  Mg_fem[trdof, :][:, bdof],  Mg_fem[trdof, :][:, bldof], Mg_fem[trdof, :][:, rdof],  Mg_fem[trdof, :][:, tdof],  Mg_fem[trdof, :][:, brdof], Mg_fem[trdof, :][:, tldof], Mg_fem[trdof, :][:, trdof]]),
                        ])

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

        U = zeros((fulldofs, mu.shape[0]), dtype=float64)
        freq_waves = zeros((modeEnd, mu.shape[0]), dtype=float64)
        # U0 = zeros((fulldofs), dtype=float64)
        for kw in range(mu.shape[0]):

            lambda_x = exp(mu[kw, 0])
            lambda_y = exp(mu[kw, 1])

            MAT_R = vstack([ 
                hstack([Iii, Zil, Zib, Zibl]),
                hstack([Zli, Ill, Zlb, Zlbl]),
                hstack([Zbi, Zbl, Ibb, Zbbl]),
                hstack([Zbli, Zbll, Zblb, Iblbl]),
                hstack([Zri, lambda_x*Irl, Zrb, Zrbl]),
                hstack([Zti, Ztl, lambda_y*Itb, Ztbl]),
                hstack([Zbri, Zbrl, Zbrb, lambda_x*Ibrbl]),
                hstack([Ztli, Ztll, Ztlb, lambda_y*Itlbl]),
                hstack([Ztri, Ztrl, Ztrb, lambda_x*lambda_y*Itrbl]),
            ])

            KG_cell_BF = (dot(dot(MAT_R.conjugate().transpose(), KG_cell), MAT_R)).tocsc()
            MG_cell_BF = (dot(dot(MAT_R.conjugate().transpose(), MG_cell), MAT_R)).tocsc()

            try:
                omega, __ = eigsh(
                    A=KG_cell_BF,
                    M=MG_cell_BF,
                    k=modeEnd,
                    sigma=1, 
                    which="LM",
                    maxiter=1000,
                )
            except:
                pass

            
            # U[:, kw] = dot(MAT_R.todense(), V[:, :modeEnd])
            freq_waves[:, kw] = sqrt(sort(real(omega)))       # Convertendo valores Reais de Omega
        
        solution["U"] = U
        solution["FREQ"] = freq_waves
        solution["IBZRANGE"] = tot_steps
        
        return solution
    
    def __setIBZ(cont_co):
        # Amostragem da IBZ
        img = 1j
        mu = array([img*cont_co[0, :].T])
        tot_steps = zeros(len(cont_co), dtype=int)

        step_size = 0.01*pi  # Defina o tamanho do passo desejado

        for i in range(len(cont_co) - 1):
            n_steps = int(ceil(sqrt((cont_co[i, 0] - cont_co[i + 1, 0])**2 +
                                        (cont_co[i, 1] - cont_co[i + 1, 1])**2) / step_size))
            
            ed = zeros((2, n_steps))
            
            for j in range(2):  # j = 0 -> x; j = 1 -> y
                step = (cont_co[i + 1, j] - cont_co[i, j]) / n_steps
                if step == 0:
                    ed[j, :] = cont_co[i, j]
                else:
                    ed[j, :] = linspace(cont_co[i, j], cont_co[i + 1, j], n_steps, endpoint=False)
            
            mu = concatenate((mu, img*ed.T), axis=0)
            tot_steps[i + 1] = tot_steps[i] + n_steps
        
        return mu, tot_steps