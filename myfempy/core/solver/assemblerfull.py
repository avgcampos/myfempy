from __future__ import annotations

from numpy import array, float64, float32, int32, zeros, empty
from scipy.sparse import coo_matrix, csc_matrix
import scipy.sparse as sp
from concurrent.futures import ThreadPoolExecutor, as_completed

INT32 = int32
FLT64 = float64
FLT32 = float32

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

        # Initialize lists to store sparse matrix data
        ith = zeros((inci.shape[0] * (elemdof * elemdof)), dtype=INT32)
        jth = zeros((inci.shape[0] * (elemdof * elemdof)), dtype=INT32)
        val = zeros((inci.shape[0] * (elemdof * elemdof)), dtype=FLT64)
        
        # ith_list, jth_list, val_list = [], [], []

        for ee in range(inci.shape[0]):
            matrix = Model.element.getStifLinearMat(Model, inci, coord, tabmat, tabgeo, intgauss, ee)
            loc = AssemblerFULL.__getLoc(Model, inci, ee)
            ith, jth, val = AssemblerFULL.__getVectorization(ith, jth, val, loc, matrix, ee, elemdof)
        
        A_sp_scipy_csr = coo_matrix((val, (ith, jth)), shape=(sdof, sdof), dtype=FLT64).tocsr()

        # AssemblerFULL.__getSaveAssemblerFile(A_sp_scipy_csr)

        return A_sp_scipy_csr

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

        A_sp_scipy_csc = coo_matrix((val, (ith, jth)), shape=(sdof, sdof))
        A_sp_scipy_csr = A_sp_scipy_csc.tocsr()

        # AssemblerFULL.__getSaveAssemblerFile(A_sp_scipy_csr)

        return A_sp_scipy_csr

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

    def __getSaveAssemblerFile(A_sp_scipy_csr):
        ## SAVE THE ASSEMBLER IN FILE
        sp.save_npz('sparse_matrix.npz', A_sp_scipy_csr)