from __future__ import annotations

from numpy import array, float64, int32, zeros, empty
from scipy.sparse import coo_matrix

INT32 = int32
FLT64 = float64

from myfempy.core.solver.assembler import Assembler
from myfempy.core.solver.assemblerfull import AssemblerFULL
from myfempy.core.solver.assemblersymm_cython import getVectorization


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

        A_sp_scipy = coo_matrix((val_band, (ith_band, jth_band)), shape=(sdof, sdof), dtype=FLT64).tocsr()
        A_sp_scipy += A_sp_scipy.transpose()
        A_sp_scipy += coo_matrix((val_diag, (ith_diag, ith_diag)), shape=(sdof, sdof), dtype=FLT64).tocsr()
        return A_sp_scipy

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
