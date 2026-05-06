from __future__ import annotations

from numpy import (arange, array, zeros_like, concatenate, setdiff1d, dot, float64, in1d, int16, int64,
                   sort, pi, where, zeros, sum, real, linspace, ceil, exp, ix_, float64)

FLT64 = float64
from scipy.sparse import csc_matrix, lil_matrix, eye, hstack, vstack, bmat
from scipy.sparse.linalg import minres, spsolve

from myfempy.core.solver.assemblerfull import AssemblerFULL
from myfempy.core.solver.assemblerfull_parallel import AssemblerFULLPOOL
from myfempy.core.solver.assemblersymm import AssemblerSYMM
from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import setSteps, gauss_points


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

class HomogenizationPlaneBCPeriodic(Solver):
    """
    Homogenization Plane Boundary Periodic Solver Class <ConcreteClassService>
    """

    def getMatrixAssembler(
        Model, inci = None, coord = None, tabmat = None, tabgeo = None, intgauss = None, SYMM=None, MP=None
    ):
        matrix = dict()
        if SYMM:
            matrix["stiffness"] = AssemblerSYMM.getLinearStiffnessGlobalMatrixAssembler(
                Model, inci, coord, tabmat, tabgeo, intgauss,
            )
        else:
            if MP:
                matrix["stiffness"] = AssemblerFULLPOOL.getLinearStiffnessGlobalMatrixAssembler(
                    Model, inci, coord, tabmat, tabgeo, intgauss,
                    MP=MP,
                )
            else:
                matrix["stiffness"] = AssemblerFULL.getLinearStiffnessGlobalMatrixAssembler(
                    Model, inci, coord, tabmat, tabgeo, intgauss,
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

        dofs_bourders = concatenate((pc_left_dof, pc_bottom_dof, pc_bottom_left_dof, pc_right_dof, pc_top_dof, pc_bottom_right_dof, pc_top_left_dof, pc_top_right_dof), axis=0)
        testpc2i = in1d(full_dofs, dofs_bourders, assume_unique=True, invert=True)
        
        # nodes_constrain_XX = constrains[where(constrains[:, 3] == 1)]
        # # nodes_constrain_YY = constrains[where(constrains[:, 3] == 2)]
        # nodes_constrain_XY = constrains[where(constrains[:, 3] == 2)]

        # testfixXX2TL = in1d(nodes_constrain_XX[:,0], pc_top_left_constrain[:,0], assume_unique=True, invert=True)
        # nodes_constrain_XX = nodes_constrain_XX[testfixXX2TL,:]

        # testfixXX2BR = in1d(nodes_constrain_XX[:,0], pc_bottom_right_constrain[:,0], assume_unique=True, invert=True)
        # nodes_constrain_XX = nodes_constrain_XX[testfixXX2BR,:]

        # testfixXY2TL = in1d(nodes_constrain_XY[:,0], pc_top_left_constrain[:,0], assume_unique=True, invert=True)
        # nodes_constrain_XY = nodes_constrain_XY[testfixXY2TL,:]

        # testfixXY2BR = in1d(nodes_constrain_XY[:,0], pc_bottom_right_constrain[:,0], assume_unique=True, invert=True)
        # nodes_constrain_XY = nodes_constrain_XY[testfixXY2BR,:]
                
        freedof = full_dofs[testpc2i]
        constdof = [pc_left_dof, pc_bottom_dof, pc_bottom_left_dof, pc_right_dof, pc_top_dof, pc_bottom_right_dof, pc_top_left_dof, pc_top_right_dof]
        fixedof = []

        return freedof, fixedof, constdof

    def getDirichletNH(constrains, nodetot, nodedof):
        return AssemblerFULL.getDirichletNH(constrains, nodetot, nodedof)

    def runSolve(Model, Physic, assembly, constrainsdof, solverset):
        ndofs = Model.modelinfo["fulldofs"]
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        ntensor = len(elem_set['tensor'])
        shape_set = Model.shape.getShapeSet()
        type_shape = shape_set["key"]

        solution = dict()
        # nsteps = setSteps(solverset["STEPSET"])
        
        idof = constrainsdof["freedof"]
        ldof = constrainsdof["constdof"][0]
        bdof = constrainsdof["constdof"][1]
        bldof = constrainsdof["constdof"][2]
        rdof = constrainsdof["constdof"][3]
        tdof = constrainsdof["constdof"][4]
        brdof = constrainsdof["constdof"][5]
        tldof = constrainsdof["constdof"][6]
        trdof = constrainsdof["constdof"][7]

        full_dofs_cell = concatenate((idof, ldof, bdof, bldof, rdof, tdof, brdof, tldof, trdof), axis=0)
        fulldof_pc_red = concatenate((idof, ldof, bdof, bldof), axis=0)
        freedof_pc = where(in1d(fulldof_pc_red, bldof, assume_unique=True) == False)[0]

        Kg_fem = assembly["stiffness"]
        Fg_fem = assembly["loads"]
                
        # --- CONDENSACAO PERIODICA INF ---

        # Criando a KG_cell de forma mais eficiente
        # O uso de np.ix_ permite extrair blocos da matriz esparsa sem fatiamentos repetitivos
        # dofs = [idof, ldof, bdof, bldof, rdof, tdof, brdof, tldof, trdof]
        # Certifique-se de que todos os arrays de DOF sejam 1D
        dofs = [d.flatten() for d in [idof, ldof, bdof, bldof, rdof, tdof, brdof, tldof, trdof]]

        blocks = []
        for row_dof in dofs:
            row_blocks = []
            for col_dof in dofs:
                # Agora o np.ix_ receberá vetores 1D e funcionará corretamente
                row_blocks.append(Kg_fem[ix_(row_dof, col_dof)])
            blocks.append(row_blocks)

        KG_cell = bmat(blocks, format='csc')

        # Ajuste da FG_cell para seguir exatamente a mesma ordem de blocos
        # Criamos uma lista de sub-vetores extraídos de Fg_fem e empilhamos verticalmente
        fg_blocks = [Fg_fem[d, :] for d in dofs]
        FG_cell = vstack(fg_blocks, format='csc')

        # --- DEFINIÇÃO DA MATRIZ DE TRANSFORMAÇÃO (MAT_R) ---
        # Definindo as dimensões para facilitar a leitura
        ni, nl, nb, nbl = idof.shape[0], ldof.shape[0], bdof.shape[0], bldof.shape[0]
        nr, nt, nbr, ntl, ntr = rdof.shape[0], tdof.shape[0], brdof.shape[0], tldof.shape[0], trdof.shape[0]

        # Matrizes Identidade e Zeros (usando float64 para evitar erros de casting)
        Iii = eye(ni, ni, dtype=float64)
        Ill = eye(nl, nl, dtype=float64)
        Ibb = eye(nb, nb, dtype=float64)
        Iblbl = eye(nbl, nbl, dtype=float64)

        # Relações de Periodicidade (Identidades que ligam as faces opostas)
        Irl = eye(nr, nl, dtype=float64)    # Right -> Left
        Itb = eye(nt, nb, dtype=float64)    # Top -> Bottom
        Ibrbl = eye(nbr, nbl, dtype=float64) # Bottom-Right -> Bottom-Left
        Itlbl = eye(ntl, nbl, dtype=float64) # Top-Left -> Bottom-Left
        Itrbl = eye(ntr, nbl, dtype=float64) # Top-Right -> Bottom-Left

        # Matriz de Restrição MAT_R usando bmat (muito mais limpo que hstack/vstack manuais)
        MAT_R = bmat([
            [Iii,             None,            None,            None],  # Interior
            [None,            Ill,             None,            None],  # Left
            [None,            None,            Ibb,             None],  # Bottom
            [None,            None,            None,            Iblbl], # Bottom-Left
            [None,            Irl,             None,            None],  # Right (vinculado à Left)
            [None,            None,            Itb,             None],  # Top (vinculado à Bottom)
            [None,            None,            None,            Ibrbl], # Bottom-Right (vinculado à BL)
            [None,            None,            None,            Itlbl], # Top-Left (vinculado à BL)
            [None,            None,            None,            Itrbl]  # Top-Right (vinculado à BL)
        ], format='csc')

        # --- CONDENSAÇÃO ESTÁTICA ---
        # Usando o operador @ para multiplicação de matrizes (Python 3.5+)
        # .T é o atalho para transpose()
        KG_cell_SC = (MAT_R.T @ KG_cell @ MAT_R).tocsc()
        FG_cell_SC = (MAT_R.T @ FG_cell).tocsc()

        U0 = zeros((KG_cell_SC.shape[0], FG_cell_SC.shape[1]), dtype=float64)  # empty((fulldofs, nsteps))
        U_FULL_PC = zeros((KG_cell.shape[0], FG_cell.shape[1]), dtype=float64)

        # fixeddofs = concatenate((bldof, brdof, trdof, tldof), axis=0) #array([0, 1])
        # reduceddofs = arange(KG_cell_SC.shape[0])
        # freedof_pc = setdiff1d(reduceddofs, fixeddofs)

        # freedof_pc = arange(2, KG_cell_SC.shape[0])
        for sslv in range(ntensor):
            try:
                 U0[freedof_pc, sslv] = spsolve(A=KG_cell_SC[:, freedof_pc][freedof_pc, :], b=FG_cell_SC[freedof_pc, sslv])
            except:
                raise 'erro'
            
            U_FULL_PC[:, sslv] = dot(MAT_R.toarray(), U0[:, sslv])

        U = zeros_like(U_FULL_PC)
        U[ix_(full_dofs_cell), :] = U_FULL_PC

        inci = Model.inci
        coord = Model.coord
        tabmat = Model.tabmat
        tabgeo = Model.tabgeo
        intgauss = Model.intgauss

        CH = zeros((ntensor, ntensor))
        rhoH = 0.0
        for elm in range(inci.shape[0]):
            nodelist = Model.shape.getNodeList(inci, elm)
            elementcoord = Model.shape.getNodeCoord(coord, nodelist)
            t = tabgeo[int(inci[elm, 3] - 1)]["THICKN"]
            Ci = Model.material.getElasticTensor(tabmat, inci, elm)
            pt, wt = gauss_points(type_shape, intgauss)
            loc = Model.shape.getLocKey(nodelist, nodedof)
            ui = U[ix_(loc)]
            CHelm = zeros((ntensor, ntensor))
            rhoHelm = 0.0
            for ip in range(intgauss):
                for jp in range(intgauss):
                    detJ = Model.shape.getdetJacobi(array([pt[ip], pt[jp]]), elementcoord)
                    diffN = Model.shape.getDiffShapeFuntion(array([pt[ip], pt[jp]]), nodedof)
                    invJ = Model.shape.getinvJacobi(array([pt[ip], pt[jp]]), elementcoord, nodedof)
                    B = Model.element.getB(diffN, invJ)
                    CHelm +=  (Ci - dot(Ci, dot(B, ui))) * t * abs(detJ) * wt[ip] * wt[jp]
                    if solverset['RHOH']:
                        R = tabmat[int(inci[elm, 2]) - 1]["RHO"]
                        rhoHelm += (R) * t * abs(detJ) * wt[ip] * wt[jp]

            CH += CHelm
            rhoH += rhoHelm

        Yx = max(coord[:,1])
        Yy = max(coord[:,2])

        CH = CH/(Yx * Yy * t)
        rhoH = rhoH/(Yx * Yy * t)

        solution["U"] = U
        solution["CH"] = CH
        solution['RHOH'] = rhoH            
        return solution