from __future__ import annotations

from numpy import (arange, array, empty, concatenate, unique, dot, float64, isin, int64,
                   sort, pi, where, zeros, sqrt, real, linspace, ceil, exp)

from scipy.sparse import coo_matrix, csc_matrix, eye, hstack, vstack
from scipy.sparse.linalg import LinearOperator, eigsh, splu, ArpackError

from myfempy.core.solver.assemblerSR import AssemblerSR
from myfempy.core.solver.assemblerMP import AssemblerMP
from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps

from concurrent.futures import ThreadPoolExecutor
from os import cpu_count

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


class PhononicCrystalPlaneBCPeriodic(Solver):
    """
    Phononic Crystal In-Plane Boundary Periodic Solver Class <ConcreteClassService>
    """
    def getMatrixAssembler(
        Model, inci = None, coord = None, tabmat = None, tabgeo = None, intgauss = None, MP = None, max_workers=None):
        matrix = dict()
        if MP:
            matrix["stiffness"] = AssemblerMP.getGlobalMatrixAssembler(
                Model, Model.element.getStifLinearMat, Model.shape.getIntNumK, inci, coord, tabmat, tabgeo, intgauss, max_workers
            )
            matrix["mass"] = AssemblerMP.getGlobalMatrixAssembler(
                Model, Model.element.getMassConsistentMat, Model.shape.getIntNumM, inci, coord, tabmat, tabgeo, intgauss, max_workers
            )
        else:
            matrix["stiffness"] = AssemblerSR.getGlobalMatrixAssembler(
                Model, Model.element.getStifLinearMat, Model.shape.getIntNumK, inci, coord, tabmat, tabgeo, intgauss,
            )
            matrix["mass"] = AssemblerSR.getGlobalMatrixAssembler(
                Model, Model.element.getMassConsistentMat, Model.shape.getIntNumM, inci, coord, tabmat, tabgeo, intgauss,
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

        testl2bl = isin(pc_left_constrain[:,0], pc_bottom_left_constrain[:,0], assume_unique=True, invert=True)
        pc_left_constrain = pc_left_constrain[testl2bl,:]
        
        testl2tl = isin(pc_left_constrain[:,0], pc_top_left_constrain[:,0], assume_unique=True, invert=True)
        pc_left_constrain = pc_left_constrain[testl2tl,:]
                        
        testb2bl = isin(pc_bottom_constrain[:,0], pc_bottom_left_constrain[:,0], assume_unique=True, invert=True)
        pc_bottom_constrain = pc_bottom_constrain[testb2bl,:]

        testb2br = isin(pc_bottom_constrain[:,0], pc_bottom_right_constrain[:,0], assume_unique=True, invert=True)
        pc_bottom_constrain = pc_bottom_constrain[testb2br,:]
        
        testt2tl = isin(pc_top_constrain[:,0], pc_top_left_constrain[:,0], assume_unique=True, invert=True)
        pc_top_constrain = pc_top_constrain[testt2tl,:]
        
        testt2tr = isin(pc_top_constrain[:,0], pc_top_right_constrain[:,0], assume_unique=True, invert=True)
        pc_top_constrain = pc_top_constrain[testt2tr,:]

        testr2tr = isin(pc_right_constrain[:,0], pc_top_right_constrain[:,0], assume_unique=True, invert=True)
        pc_right_constrain = pc_right_constrain[testr2tr,:]

        testr2br = isin(pc_right_constrain[:,0], pc_bottom_right_constrain[:,0], assume_unique=True, invert=True)
        pc_right_constrain = pc_right_constrain[testr2br,:]

        __, pc_left_dof, __ = AssemblerSR.getConstrains(pc_left_constrain, nodetot, nodedof)
        __, pc_right_dof, __ = AssemblerSR.getConstrains(pc_right_constrain, nodetot, nodedof)
        __, pc_bottom_dof, __ = AssemblerSR.getConstrains(pc_bottom_constrain, nodetot, nodedof)
        __, pc_top_dof, __ = AssemblerSR.getConstrains(pc_top_constrain, nodetot, nodedof)
        __, pc_bottom_left_dof, __ = AssemblerSR.getConstrains(pc_bottom_left_constrain, nodetot, nodedof)
        __, pc_bottom_right_dof, __ = AssemblerSR.getConstrains(pc_bottom_right_constrain, nodetot, nodedof)
        __, pc_top_left_dof, __ = AssemblerSR.getConstrains(pc_top_left_constrain, nodetot, nodedof)
        __, pc_top_right_dof, __ = AssemblerSR.getConstrains(pc_top_right_constrain, nodetot, nodedof)

        full_dofs = arange(0, nodedof * nodetot, 1, int)

        pc_dofs = concatenate((pc_left_dof, pc_bottom_dof, pc_bottom_left_dof, pc_right_dof, pc_top_dof, pc_bottom_right_dof, pc_top_left_dof, pc_top_right_dof), axis=0)
        testpc2i = isin(full_dofs, pc_dofs, assume_unique=True, invert=True)
        freedof = full_dofs[testpc2i]

        constdof = [pc_left_dof, pc_bottom_dof, pc_bottom_left_dof, pc_right_dof, pc_top_dof, pc_bottom_right_dof, pc_top_left_dof, pc_top_right_dof]
        fixedof = []

        return freedof, fixedof, constdof

    def getDirichletNH(constrains, nodetot, nodedof):
        return AssemblerSR.getDirichletNH(constrains, nodetot, nodedof)


    def runSolve(Model, Physic, assembly, constrainsdof, solverset):
        fulldofs = Model.modelinfo["fulldofs"]
        modeEnd = setSteps(solverset["STEPSET"])
        mu, tot_steps = PhononicCrystalPlaneBCPeriodic.__setIBZ(solverset["IBZ"])

        # Grupos de DOFs: internos, esquerda, base, base-esq, direita, topo, base-dir, topo-esq, topo-dir
        grupos = [constrainsdof["freedof"], *constrainsdof["constdof"][:8]]
        ordem = concatenate(grupos)

        KG_cell = assembly["stiffness"][ordem][:, ordem]
        MG_cell = assembly["mass"][ordem][:, ordem]

        R0, A, B, C = PhononicCrystalPlaneBCPeriodic.__setDenseMat([g.shape[0] for g in grupos])

        n_mu = mu.shape[0]

        # Parâmetros opcionais do solver
        sigma = 0.0
        n_workers = solverset["NCPU"]
        def resolver_ponto(kmu):
            """Frequências de um ponto da IBZ (seguro para threads: só lê dados compartilhados)."""
            lambda_x = exp(mu[kmu, 0])
            lambda_y = exp(mu[kmu, 1])

            MAT_R = R0 + lambda_x * A + lambda_y * B + lambda_x * lambda_y * C
            MAT_R_H = MAT_R.conjugate().transpose()

            KG_cell_BF = (MAT_R_H.dot(KG_cell.dot(MAT_R))).tocsc()
            MG_cell_BF = (MAT_R_H.dot(MG_cell.dot(MAT_R))).tocsc()
            del MAT_R, MAT_R_H                                   # libera memória antes da fatoração

            try:
                # Fatoração LU de (K - sigma*M) com ordenação para estrutura simétrica
                lu = splu(
                    (KG_cell_BF - sigma * MG_cell_BF).tocsc(),
                    permc_spec="MMD_AT_PLUS_A",
                    diag_pivot_thresh=0.0,                       # pivô na diagonal (K Hermitiana)
                    options={"SymmetricMode": True},
                )
                OPinv = LinearOperator(KG_cell_BF.shape, matvec=lu.solve,
                                       dtype=KG_cell_BF.dtype)

                avec = eigsh(A=KG_cell_BF, M=MG_cell_BF, k=modeEnd, sigma=sigma,
                             which="LM", OPinv=OPinv, return_eigenvectors=False)
            except (ArpackError, RuntimeError):
                return 0.0

            return sqrt(sort(real(avec)))

        freq_waves = zeros((modeEnd, n_mu), dtype=float64)

        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            for kmu, freq in enumerate(pool.map(resolver_ponto, range(n_mu))):
                freq_waves[:, kmu] = freq

        solution = {"U": zeros((fulldofs, n_mu), dtype=float64)}
        solution["FREQ"] = freq_waves
        solution["IBZRANGE"] = tot_steps
        return solution

    def __setAssembleBlocks(shape, blocks):
        """Monta uma matriz esparsa CSC posicionando blocos (linha0, coluna0, matriz)."""
        coos = [(r0, c0, coo_matrix(m)) for r0, c0, m in blocks]
        rows = concatenate([m.row + r0 for r0, _, m in coos])
        cols = concatenate([m.col + c0 for _, c0, m in coos])
        data = concatenate([m.data for _, _, m in coos])
        return coo_matrix((data, (rows, cols)), shape=shape).tocsc()


    def __setDenseMat(sizes):
        """
        Pré-calcula as matrizes esparsas tais que, para cada ponto da IBZ:
            MAT_R = R0 + lambda_x*A + lambda_y*B + lambda_x*lambda_y*C
        Ordem dos graus de liberdade: i, l, b, bl, r, t, br, tl, tr
        """
        ni, nl, nb, nbl, nr, nt, nbr, ntl, ntr = sizes
        off = [0]
        for n in sizes:
            off.append(off[-1] + n)

        shape = (off[-1], ni + nl + nb + nbl)

        R0 = PhononicCrystalPlaneBCPeriodic.__setAssembleBlocks(shape, [(off[k], off[k], eye(sizes[k])) for k in range(4)])
        A = PhononicCrystalPlaneBCPeriodic.__setAssembleBlocks(shape, [(off[4], off[1], eye(nr, nl)),     # r  <- l   (lambda_x)
                                (off[6], off[3], eye(nbr, nbl))])  # br <- bl  (lambda_x)
        B = PhononicCrystalPlaneBCPeriodic.__setAssembleBlocks(shape, [(off[5], off[2], eye(nt, nb)),     # t  <- b   (lambda_y)
                                (off[7], off[3], eye(ntl, nbl))])  # tl <- bl  (lambda_y)
        C = PhononicCrystalPlaneBCPeriodic.__setAssembleBlocks(shape, [(off[8], off[3], eye(ntr, nbl))])  # tr <- bl  (lambda_x*lambda_y)
        return R0, A, B, C
    
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