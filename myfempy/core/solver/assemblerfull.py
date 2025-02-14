from __future__ import annotations

from numpy import array, float64, int32, zeros, empty
from scipy.sparse import coo_matrix, csc_matrix
from concurrent.futures import ThreadPoolExecutor, as_completed

INT32 = int32
FLT64 = float64

from myfempy.core.solver.assembler import Assembler
from myfempy.core.solver.assemblerfull_numpy_v1 import (getConstrains,
                                                        getDirichletNH,
                                                        getLoadAssembler)

from myfempy.core.solver.assemblerfull_cython_v5 import getVectorization

class AssemblerFULL(Assembler):
    """
    Assembler Full System Class <ConcreteClassService>
    """

    # @profile
    def getLinearStiffnessGlobalMatrixAssembler(
        Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler, MP
    ):
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        elemdof = nodecon * nodedof
        nodetot = coord.shape[0]
        sdof = nodedof * nodetot

        # Initialize lists to store sparse matrix data
        ith = zeros((inci.shape[0] * (elemdof * elemdof)), dtype=INT32)
        jth = zeros((inci.shape[0] * (elemdof * elemdof)), dtype=INT32)
        val = zeros((inci.shape[0] * (elemdof * elemdof)), dtype=FLT64)
        
        # ith_list, jth_list, val_list = [], [], []

        for ee in range(inci.shape[0]):
            matrix = Model.element.getStifLinearMat(Model, inci, coord, tabmat, tabgeo, intgauss, ee)
            loc = AssemblerFULL.__getLoc(Model, inci, ee)
            ith, jth, val = AssemblerFULL.__getVectorization(ith, jth, val, loc, matrix, ee, elemdof)

        A_sp_scipy_csc = csc_matrix((val, (ith, jth)), shape=(sdof, sdof))
        return A_sp_scipy_csc

    def getNonLinearStiffnessGlobalMatrixAssembler():
        pass

    def getMassConsistentGlobalMatrixAssembler(
        Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler, MP
    ):
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        elemdof = nodecon * nodedof
        nodetot = coord.shape[0]
        sdof = nodedof * nodetot

        # Initialize lists to store sparse matrix data
        ith = empty((inci.shape[0] * (elemdof * elemdof)), dtype=INT32)
        jth = empty((inci.shape[0] * (elemdof * elemdof)), dtype=INT32)
        val = empty((inci.shape[0] * (elemdof * elemdof)), dtype=FLT64)

        for ee in range(inci.shape[0]):
            matrix = Model.element.getMassConsistentMat(
                Model, inci, coord, tabmat, tabgeo, intgauss, ee
            )
            loc = AssemblerFULL.__getLoc(Model, inci, ee)
            ith, jth, val = AssemblerFULL.__getVectorization(
                ith, jth, val, loc, matrix, ee, elemdof
            )

        A_sp_scipy_csc = csc_matrix((val, (ith, jth)), shape=(sdof, sdof))
        return A_sp_scipy_csc

    def getMassLumpedGlobalMatrixAssembler():
        pass

    def getLoadAssembler(loadaply, nodetot, nodedof):
        return getLoadAssembler(loadaply, nodetot, nodedof)

    # Dirichlet Homogeneous https://en.wikipedia.org/wiki/Dirichlet_boundary_condition
    def getConstrains(constrains, nodetot, nodedof):
        return getConstrains(constrains, nodetot, nodedof)

    # Dirichlet Non-Homogeneous
    def getDirichletNH(constrains, nodetot, nodedof):
        return getDirichletNH(constrains, nodetot, nodedof)

    # @profile
    def __getVectorization(ith, jth, val, loc, matrix, element_number, elemdof):
        return getVectorization(ith, jth, val, loc, matrix, element_number, elemdof)

    def __getLoc(Model, inci, element_number):
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        nodelist = Model.shape.getNodeList(inci, element_number)
        loc = Model.shape.getLocKey(nodelist, nodedof)
        return array(loc)
        