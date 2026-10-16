from __future__ import annotations

from numpy import zeros, ones, array, concatenate, isin, arange, where, ravel, flatnonzero, ascontiguousarray, float64, int64

FLT64 = float64
INT64 = int64
from scipy.sparse import bmat, eye, vstack, coo_matrix, issparse
from scipy.sparse.linalg import spsolve, splu

from myfempy.core.solver.assemblerfull import AssemblerFULL
from myfempy.core.solver.solver import Solver
from myfempy.core.utilities import gauss_points, setSteps

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


class HomogenizationPlaneBCPeriodic(Solver):
    """
    Homogenization Plane Boundary Periodic Solver Class <ConcreteClassService>
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
            matrix["stiffness"] = AssemblerFULL.getGlobalMatrixAssemblerMP(
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
            matrix["stiffness"] = AssemblerFULL.getGlobalMatrixAssembler(
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

        testl2bl = isin(
            pc_left_constrain[:, 0],
            pc_bottom_left_constrain[:, 0],
            assume_unique=True,
            invert=True,
        )
        pc_left_constrain = pc_left_constrain[testl2bl, :]

        testl2tl = isin(
            pc_left_constrain[:, 0],
            pc_top_left_constrain[:, 0],
            assume_unique=True,
            invert=True,
        )
        pc_left_constrain = pc_left_constrain[testl2tl, :]

        testb2bl = isin(
            pc_bottom_constrain[:, 0],
            pc_bottom_left_constrain[:, 0],
            assume_unique=True,
            invert=True,
        )
        pc_bottom_constrain = pc_bottom_constrain[testb2bl, :]

        testb2br = isin(
            pc_bottom_constrain[:, 0],
            pc_bottom_right_constrain[:, 0],
            assume_unique=True,
            invert=True,
        )
        pc_bottom_constrain = pc_bottom_constrain[testb2br, :]

        testt2tl = isin(
            pc_top_constrain[:, 0],
            pc_top_left_constrain[:, 0],
            assume_unique=True,
            invert=True,
        )
        pc_top_constrain = pc_top_constrain[testt2tl, :]

        testt2tr = isin(
            pc_top_constrain[:, 0],
            pc_top_right_constrain[:, 0],
            assume_unique=True,
            invert=True,
        )
        pc_top_constrain = pc_top_constrain[testt2tr, :]

        testr2tr = isin(
            pc_right_constrain[:, 0],
            pc_top_right_constrain[:, 0],
            assume_unique=True,
            invert=True,
        )
        pc_right_constrain = pc_right_constrain[testr2tr, :]

        testr2br = isin(
            pc_right_constrain[:, 0],
            pc_bottom_right_constrain[:, 0],
            assume_unique=True,
            invert=True,
        )
        pc_right_constrain = pc_right_constrain[testr2br, :]

        __, pc_left_dof, __ = AssemblerFULL.getConstrains(
            pc_left_constrain, nodetot, nodedof
        )
        __, pc_right_dof, __ = AssemblerFULL.getConstrains(
            pc_right_constrain, nodetot, nodedof
        )
        __, pc_bottom_dof, __ = AssemblerFULL.getConstrains(
            pc_bottom_constrain, nodetot, nodedof
        )
        __, pc_top_dof, __ = AssemblerFULL.getConstrains(
            pc_top_constrain, nodetot, nodedof
        )
        __, pc_bottom_left_dof, __ = AssemblerFULL.getConstrains(
            pc_bottom_left_constrain, nodetot, nodedof
        )
        __, pc_bottom_right_dof, __ = AssemblerFULL.getConstrains(
            pc_bottom_right_constrain, nodetot, nodedof
        )
        __, pc_top_left_dof, __ = AssemblerFULL.getConstrains(
            pc_top_left_constrain, nodetot, nodedof
        )
        __, pc_top_right_dof, __ = AssemblerFULL.getConstrains(
            pc_top_right_constrain, nodetot, nodedof
        )

        full_dofs = arange(0, nodedof * nodetot, 1, int)

        dofs_bourders = concatenate(
            (
                pc_left_dof,
                pc_bottom_dof,
                pc_bottom_left_dof,
                pc_right_dof,
                pc_top_dof,
                pc_bottom_right_dof,
                pc_top_left_dof,
                pc_top_right_dof,
            ),
            axis=0,
        )
        testpc2i = isin(full_dofs, dofs_bourders, assume_unique=True, invert=True)

        freedof = full_dofs[testpc2i]
        constdof = [
            pc_left_dof,
            pc_bottom_dof,
            pc_bottom_left_dof,
            pc_right_dof,
            pc_top_dof,
            pc_bottom_right_dof,
            pc_top_left_dof,
            pc_top_right_dof,
        ]
        fixedof = []

        return freedof, fixedof, constdof

    def getDirichletNH(constrains, nodetot, nodedof):
        return AssemblerFULL.getDirichletNH(constrains, nodetot, nodedof)

    def runSolve(Model, Physic, assembly, constrainsdof, solverset):
        elem_set = Model.element.getElementSet()
        H = elem_set["H"]
        nodedof = Model.modelinfo["nodedof"]
        ntensor = Model.modelinfo["tensor"]
        type_shape = Model.shape.getShapeSet()["key"]

        cd = constrainsdof["constdof"]
        idof = constrainsdof["freedof"]
        ldof, bdof, bldof, rdof, tdof, brdof, tldof, trdof = cd[:8]

        # ordem completa dos dofs (i, l, b, bl, r, t, br, tl, tr)
        dofs = [ravel(d) for d in
                (idof, ldof, bdof, bldof, rdof, tdof, brdof, tldof, trdof)]
        perm = concatenate(dofs)
        full_dofs_cell = concatenate((idof, ldof, bdof, bldof, rdof, tdof, brdof, tldof, trdof), axis=0)

        fulldof_pc_red = concatenate((idof, ldof, bdof, bldof), axis=0)
        freedof_pc = flatnonzero(~isin(fulldof_pc_red, bldof, assume_unique=True))

        Kg_fem = assembly["stiffness"].tocsr()
        Fg_fem = assembly["loads"]

        # --- REORDENAÇÃO (substitui as 81 extrações ix_ + bmat) ---
        KG_cell = Kg_fem[perm][:, perm]
        FG_cell = Fg_fem[perm, :]
        if issparse(FG_cell):
            FG_cell = FG_cell.toarray()

        # --- MAT_R montada direto em COO (uma entrada 1.0 por dof) ---
        ni, nl, nb, nbl = (len(d) for d in dofs[:4])
        nr, nt, nbr, ntl, ntr = (len(d) for d in dofs[4:])
        assert nr == nl and nt == nb and nbr == ntl == ntr == nbl, \
            "dimensões periódicas inconsistentes"

        off_i, off_l = 0, ni
        off_b, off_bl = ni + nl, ni + nl + nb
        cols = concatenate((
            arange(ni) + off_i,        # interior
            arange(nl) + off_l,        # left
            arange(nb) + off_b,        # bottom
            arange(nbl) + off_bl,      # bottom-left
            arange(nr) + off_l,        # right  -> left
            arange(nt) + off_b,        # top    -> bottom
            arange(nbr) + off_bl,      # br     -> bl
            arange(ntl) + off_bl,      # tl     -> bl
            arange(ntr) + off_bl,      # tr     -> bl
        ))
        nfull, nred = perm.size, ni + nl + nb + nbl
        MAT_R = coo_matrix((ones(nfull), (arange(nfull), cols)),
                        shape=(nfull, nred)).tocsr()
        MAT_RT = MAT_R.T.tocsr()

        # --- CONDENSAÇÃO ESTÁTICA ---
        KG_SC = (MAT_RT @ KG_cell @ MAT_R).tocsr()
        FG_SC = MAT_RT @ FG_cell                      # denso (nred, ntensor)

        # --- SOLUÇÃO: uma fatoração LU, todos os ntensor RHS de uma vez ---
        Kff = KG_SC[freedof_pc][:, freedof_pc].tocsc()
        U0 = zeros((nred, FG_SC.shape[1]), dtype=float64)
        U0[freedof_pc, :] = splu(Kff).solve(ascontiguousarray(FG_SC[freedof_pc, :]))

        U_FULL_PC = MAT_R @ U0                        # esparso x denso, sem toarray()
        U = zeros_like_U = zeros(U_FULL_PC.shape, dtype=float64)
        U[full_dofs_cell, :] = U_FULL_PC

        # --- HOMOGENEIZAÇÃO ---
        inci, coord = Model.inci, Model.coord
        tabmat, tabgeo, intgauss = Model.tabmat, Model.tabgeo, Model.intgauss
        sh = Model.shape
        getNodeList, getNodeCoord = sh.getNodeList, sh.getNodeCoord
        getdetJacobi, getinvJacobi = sh.getdetJacobi, sh.getinvJacobi
        getLocKey, getDiffShapeFuntion, getB = sh.getLocKey, sh.getDiffShapeFuntion, sh.getB
        getElasticTensor = Model.material.getElasticTensor

        # invariantes do laço: pontos/pesos de Gauss e dN (não dependem do elemento)
        pt, wt = gauss_points(type_shape, intgauss)
        gp = [(array([pt[i], pt[j]]), wt[i] * wt[j])
            for i in range(intgauss) for j in range(intgauss)]
        gp = [(p, w, getDiffShapeFuntion(p, nodedof)) for p, w in gp]

        use_rho = bool(solverset["RHOH"])
        geo_id = inci[:, 3].astype(INT64) - 1
        mat_id = inci[:, 2].astype(INT64) - 1

        CH = zeros((ntensor, ntensor))
        rhoH = 0.0
        t = None
        for elm in range(inci.shape[0]):
            nodelist = getNodeList(inci, elm)
            elementcoord = getNodeCoord(coord, nodelist)
            t = tabgeo[geo_id[elm]]["THICKN"]
            Ci = getElasticTensor(tabmat, inci, elm)
            ui = U[getLocKey(nodelist, nodedof)]

            # acumula B ponderado e volume nos pontos de Gauss; um único dot por elemento
            Bsum = 0.0
            vol = 0.0
            for p, w, diffN in gp:
                f = getdetJacobi(p, elementcoord) * t * w
                B = getB(H, getinvJacobi(p, elementcoord, nodedof), diffN)
                Bsum = Bsum + B * f
                vol += f

            CH += Ci * vol - Ci @ (Bsum @ ui)
            if use_rho:
                rhoH += tabmat[mat_id[elm]]["RHO"] * vol

        Yx = coord[:, 1].max()
        Yy = coord[:, 2].max()
        den = Yx * Yy * t          # t = espessura do último elemento (comportamento original)
        solution = {"U": U, "CH": CH / den, "RHOH": rhoH / den}
        return solution