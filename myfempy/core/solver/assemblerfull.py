from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from scipy.sparse import coo_matrix
from numpy import zeros, float64, float32, int32

INT32 = int32
FLT64 = float64
FLT32 = float32

from myfempy.core.solver.assembler import Assembler
from myfempy.core.utilities import gauss_points
from myfempy.core.solver.assemblerfull_cython import getVectorization                                         
from myfempy.core.solver.assemblerfull_numpy import (getConstrains,
                                                    getDirichletNH,
                                                    getLoadAssembler)

__docformat__ = "google"

__doc__ = """

==========================================================================
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

    def getGlobalMatrixAssemblerMP(Model, getElemMatrix, getIntNum, inci=None, coord=None, 
                             tabmat=None, tabgeo=None, intgauss=None, max_workers=None):
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

        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        type_shape = shape_set["key"]
        nodecon = len(shape_set["nodes"])
        elemdof = nodecon * nodedof
        nodetot = coord.shape[0]
        sdof = nodedof * nodetot

        nelem = inci.shape[0]
        entries_per_elem = elemdof * elemdof
        total_entries = nelem * entries_per_elem

        # Pré-alocação dos vetores Globais em C-contiguous memory
        ith = zeros(total_entries, dtype=INT32)
        jth = zeros(total_entries, dtype=INT32)
        val = zeros(total_entries, dtype=FLT64)
        
        getNodeList = Model.shape.getNodeList
        getLocKey = Model.shape.getLocKey
        getElasticTensor = Model.material.getElasticTensor
        getNodeCoord = Model.shape.getNodeCoord

        pt, wt = gauss_points(type_shape, intgauss)
        def _process_element_chunk(start_elem: int, end_elem: int):
            for ee in range(start_elem, end_elem):
                nodelist = getNodeList(inci, ee)
                elementcoord = getNodeCoord(coord, nodelist)
                C = getElasticTensor(tabmat, inci, ee)
                matrix = getElemMatrix(inci, coord, tabmat, tabgeo, elementcoord, C, elemdof, getIntNum, intgauss, pt, wt, ee)
                loc = getLocKey(nodelist, nodedof)
                getVectorization(ith, jth, val, loc, matrix, ee, elemdof)

        # Divisão de blocos entre as threads (Python 3.14t Free-Threading)
        if max_workers is None:
            max_workers = 4
        chunk_size = max(1, nelem // (max_workers * 4))
        chunks = [(i, min(i + chunk_size, nelem)) for i in range(0, nelem, chunk_size)]

        print(f"[Paralelismo] {nelem} elementos divididos entre {max_workers} núcleos (~{nelem // max_workers} elem/núcleo em {len(chunks)} blocos de {chunk_size}).")

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(_process_element_chunk, start, end) for start, end in chunks]
            for future in futures:
                future.result()

        # Construção da matriz esparsa final
        A_sp_scipy_csr = coo_matrix((val, (ith, jth)), shape=(sdof, sdof), dtype=FLT64).tocsr()
        return A_sp_scipy_csr

    def getGlobalMatrixAssembler(Model, getElemMatrix, getIntNum, inci = None, coord = None, tabmat = None, tabgeo = None, intgauss = None):
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

        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        type_shape = shape_set["key"]
        nodecon = len(shape_set["nodes"])
        elemdof = nodecon * nodedof
        nodetot = coord.shape[0]
        sdof = nodedof * nodetot

        ith = zeros((inci.shape[0] * (elemdof * elemdof)), dtype=INT32)
        jth = zeros((inci.shape[0] * (elemdof * elemdof)), dtype=INT32)
        val = zeros((inci.shape[0] * (elemdof * elemdof)), dtype=FLT64)
        
        getNodeList = Model.shape.getNodeList
        getLocKey = Model.shape.getLocKey
        getElasticTensor = Model.material.getElasticTensor
        getNodeCoord = Model.shape.getNodeCoord

        pt, wt = gauss_points(type_shape, intgauss)

        for ee in range(inci.shape[0]):
            nodelist = getNodeList(inci, ee)
            elementcoord = getNodeCoord(coord, nodelist)
            C = getElasticTensor(tabmat, inci, ee)
            matrix = getElemMatrix(inci, coord, tabmat, tabgeo, elementcoord, C, elemdof, getIntNum, intgauss, pt, wt, ee)
            loc = getLocKey(nodelist, nodedof)
            ith, jth, val = getVectorization(ith, jth, val, loc, matrix, ee, elemdof)
        
        A_sp_scipy_csr = coo_matrix((val, (ith, jth)), shape=(sdof, sdof), dtype=FLT64).tocsr()
        return A_sp_scipy_csr


    def getLoadAssembler(loadaply, nodetot, nodedof):
        return getLoadAssembler(loadaply, nodetot, nodedof)

    # Dirichlet Homogeneous https://en.wikipedia.org/wiki/Dirichlet_boundary_condition
    def getConstrains(constrains, nodetot, nodedof):
        return getConstrains(constrains, nodetot, nodedof)

    # Dirichlet Non-Homogeneous
    def getDirichletNH(constrains, nodetot, nodedof):
        return getDirichletNH(constrains, nodetot, nodedof)