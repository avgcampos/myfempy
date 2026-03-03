from __future__ import annotations

from numpy import array, float64, float32, int32, zeros, empty
from scipy.sparse import coo_matrix, csc_matrix
import scipy.sparse as sp
from concurrent.futures import ThreadPoolExecutor, as_completed

INT32 = int32
FLT64 = float64
FLT32 = float32

from myfempy.core.solver.assembler import Assembler
from myfempy.core.solver.assemblerfull_numpy import (getConstrains,
                                                        getDirichletNH,
                                                        getLoadAssembler)

from myfempy.core.solver.assemblerfull_cython import getVectorization


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

        A_sp_scipy_csr = coo_matrix((val, (ith, jth)), shape=(sdof, sdof), dtype=FLT64).tocsr()
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