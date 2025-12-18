from __future__ import annotations

from numpy import array, float64, int32, zeros, empty
from scipy.sparse import coo_matrix

INT32 = int32
FLT64 = float64

from myfempy.core.solver.assembler import Assembler
from myfempy.core.solver.assemblerfull import AssemblerFULL
from myfempy.core.solver.assemblersymm_cython_v5 import getVectorization


class AssemblerSYMM(Assembler):
    """
    Assembler Symmetric Banded System Class <ConcreteClassService>
    """

    # @profile
    def getLinearStiffnessGlobalMatrixAssembler(Model, inci = None, coord = None, tabmat = None, tabgeo = None, intgauss = None):
        
        if inci is None:
            inci = Model.inci
        
        if coord is None:
            coord = Model.coord
        
        if tabmat is None:
            tabmat = Model.tabmat
       
        if tabgeo is None:
            tabgeo = Model.tabgeo
       
        if intgauss is None:
            intgauss = Model.intgauss
       
        else:
            pass

        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        elemdof = nodecon * nodedof
        nodetot = coord.shape[0]
        sdof = nodedof * nodetot

        dim_band = int(0.5 * (elemdof * elemdof - elemdof) * inci.shape[0])
        dim_diag = int(elemdof * inci.shape[0])

        ith_band = zeros((dim_band,), dtype=INT32)
        jth_band = zeros((dim_band,), dtype=INT32)
        val_band = zeros((dim_band,), dtype=FLT64)
        ith_diag = zeros((dim_diag,), dtype=INT32)
        val_diag = zeros((dim_diag,), dtype=FLT64)

        nb = int(0)
        nd = int(0)
        for ee in range(inci.shape[0]):
            matrix = Model.element.getStifLinearMat(
                Model, inci, coord, tabmat, tabgeo, intgauss, ee
            )
            loc = AssemblerSYMM.__getLoc(Model, inci, ee)
            ith_diag, val_diag, ith_band, jth_band, val_band, nb, nd = (
                AssemblerSYMM.__getVectorization(
                    ith_band,
                    jth_band,
                    val_band,
                    ith_diag,
                    val_diag,
                    nb,
                    nd,
                    loc,
                    matrix,
                    ee,
                    elemdof,
                )
            )

        A_sp_scipy = coo_matrix((val_band, (ith_band, jth_band)), shape=(sdof, sdof), dtype=FLT64).tocsr()
        A_sp_scipy += A_sp_scipy.transpose()
        A_sp_scipy += coo_matrix((val_diag, (ith_diag, ith_diag)), shape=(sdof, sdof), dtype=FLT64).tocsr()

        return A_sp_scipy

    def getNonLinearStiffnessGlobalMatrixAssembler():
        pass

    def getMassConsistentGlobalMatrixAssembler(Model, inci = None, coord = None, tabmat = None, tabgeo = None, intgauss = None):
        
        if inci is None:
            inci = Model.inci
        
        if coord is None:
            coord = Model.coord
        
        if tabmat is None:
            tabmat = Model.tabmat
       
        if tabgeo is None:
            tabgeo = Model.tabgeo
       
        if intgauss is None:
            intgauss = Model.intgauss
       
        else:
            pass

        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        elemdof = nodecon * nodedof
        nodetot = coord.shape[0]
        sdof = nodedof * nodetot

        dim_band = int(0.5 * (elemdof * elemdof - elemdof) * inci.shape[0])
        dim_diag = int(elemdof * inci.shape[0])

        ith_band = empty((dim_band,), dtype=INT32)
        jth_band = empty((dim_band,), dtype=INT32)
        val_band = empty((dim_band,), dtype=FLT64)
        ith_diag = empty((dim_diag,), dtype=INT32)
        val_diag = empty((dim_diag,), dtype=FLT64)

        nb = int(0)
        nd = int(0)
        for ee in range(inci.shape[0]):
            matrix = Model.element.getMassConsistentMat(
                Model, inci, coord, tabmat, tabgeo, intgauss, ee
            )
            loc = AssemblerSYMM.__getLoc(Model, inci, ee)
            ith_diag, val_diag, ith_band, jth_band, val_band, nb, nd = (
                AssemblerSYMM.__getVectorization(
                    ith_band,
                    jth_band,
                    val_band,
                    ith_diag,
                    val_diag,
                    nb,
                    nd,
                    loc,
                    matrix,
                    ee,
                    elemdof,
                )
            )

        mtKG_sp_sym = coo_matrix(
            (val_band, (ith_band, jth_band)), shape=(sdof, sdof)
        ).tocsr()
        mtKG_sp_sym += mtKG_sp_sym.transpose()
        mtKG_sp_sym += coo_matrix(
            (val_diag, (ith_diag, ith_diag)), shape=(sdof, sdof)
        ).tocsr()

        return mtKG_sp_sym

    def getMassLumpedGlobalMatrixAssembler():
        pass

    def getLoadAssembler(loadaply, nodetot, nodedof):
        return AssemblerFULL.getLoadAssembler(loadaply, nodetot, nodedof)

    def getConstrains(constrains, nodetot, nodedof, Uc):
        return AssemblerFULL.getConstrains(constrains, nodetot, nodedof, Uc)

    def getDirichletNH(constrains, nodetot, nodedof):
        return AssemblerFULL.getDirichletNH(constrains, nodetot, nodedof)

    # @profile
    def __getVectorization(
        ith_band,
        jth_band,
        val_band,
        ith_diag,
        val_diag,
        nb,
        nd,
        loc,
        matrix,
        ee,
        elemdof,
    ):
        return getVectorization(
            ith_band,
            jth_band,
            val_band,
            ith_diag,
            val_diag,
            nb,
            nd,
            loc,
            matrix,
            ee,
            elemdof,
        )

    def __getLoc(Model, inci, element_number):
        elem_set = Model.element.getElementSet()
        nodedof = int(len(elem_set["dofs"]["d"]))
        nodelist = Model.shape.getNodeList(inci, element_number)
        loc = Model.shape.getLocKey(nodelist, nodedof)
        return array(loc)
