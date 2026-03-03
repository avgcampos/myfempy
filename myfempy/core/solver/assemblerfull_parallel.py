from __future__ import annotations

import os
from numpy import array, float64, int32, zeros, empty, append
from scipy.sparse import coo_matrix, csc_matrix
from concurrent.futures import ThreadPoolExecutor, as_completed
# from concurrent.futures import ProcessPoolExecutor

INT32 = int32
FLT64 = float64

from myfempy.core.solver.assembler import Assembler
from myfempy.core.solver.assemblerfull_numpy import (getConstrains,
                                                        getDirichletNH,
                                                        getLoadAssembler
                                                        )
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


class AssemblerFULLPOOL(Assembler):
    """
    Assembler Full System Class <ConcreteClassService>
    """
    # @profile
    def getLinearStiffnessGlobalMatrixAssembler(Model, inci = None, coord = None, tabmat = None, tabgeo = None, intgauss = None, MP=None):
        inci = Model.inci
        coord = Model.coord
        tabmat = Model.tabmat
        tabgeo = Model.tabgeo
        intgauss = Model.intgauss

        # Get element and shape properties
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        elemdof = nodecon * nodedof
        nodetot = coord.shape[0]
        sdof = nodedof * nodetot

        # Estimate the size of the output (adjust based on your problem)
        # Each element contributes elemdof * elemdof entries to the global matrix
        estimated_size = inci.shape[0] * elemdof * elemdof

        # Preallocate memory using NumPy arrays
        ith_array = zeros((estimated_size), dtype=INT32)
        jth_array = zeros((estimated_size), dtype=INT32)
        val_array = zeros((estimated_size), dtype=FLT64)

        # Track the current position in the arrays
        current_index = 0

        # Define the batch size for parallel processing
        batch_size = 64  # Adjust based on memory constraints and performance testing

        # Use ThreadPoolExecutor for parallel processing
        num_cores = MP if MP is not None else os.cpu_count()
        with ThreadPoolExecutor(max_workers = num_cores) as executor:
            # Submit tasks for each element in batches to reduce memory overhead
            for i in range(0, inci.shape[0], batch_size):
                batch = range(i, min(i + batch_size, inci.shape[0]))
                futures = [
                    executor.submit(
                        process_element, ee, Model, inci, coord, tabmat, tabgeo, intgauss, elemdof
                    )
                    for ee in batch
                ]

                # Collect results in batches
                batch_results = []
                for future in as_completed(futures):
                    batch_results.append(future.result())

                    # Process the batch when it reaches a certain size
                    if len(batch_results) >= batch_size:
                        ith_batch, jth_batch, val_batch = zip(*batch_results)
                        batch_length = sum(len(x) for x in ith_batch)

                        # Ensure we don't exceed the preallocated size
                        if current_index + batch_length > estimated_size:
                            ith_array.resize(current_index + batch_length)
                            jth_array.resize(current_index + batch_length)
                            val_array.resize(current_index + batch_length)

                        # Append the batch to the arrays
                        for ith, jth, val in zip(ith_batch, jth_batch, val_batch):
                            ith_array[current_index:current_index + len(ith)] = ith
                            jth_array[current_index:current_index + len(jth)] = jth
                            val_array[current_index:current_index + len(val)] = val
                            current_index += len(ith)

                        batch_results.clear()  # Clear the batch for the next set of results

                # Process any remaining results in the batch
                if batch_results:
                    ith_batch, jth_batch, val_batch = zip(*batch_results)
                    batch_length = sum(len(x) for x in ith_batch)

                    # Ensure we don't exceed the preallocated size
                    if current_index + batch_length > estimated_size:
                        ith_array.resize(current_index + batch_length)
                        jth_array.resize(current_index + batch_length)
                        val_array.resize(current_index + batch_length)

                    # Append the batch to the arrays
                    for ith, jth, val in zip(ith_batch, jth_batch, val_batch):
                        ith_array[current_index:current_index + len(ith)] = ith
                        jth_array[current_index:current_index + len(jth)] = jth
                        val_array[current_index:current_index + len(val)] = val
                        current_index += len(ith)

        # Trim the arrays to the actual size
        ith_array = ith_array[:current_index]
        jth_array = jth_array[:current_index]
        val_array = val_array[:current_index]

        # Construct the global stiffness matrix in COO format
        A_sp_scipy_csc = csc_matrix((val_array, (ith_array, jth_array)), shape=(sdof, sdof))
        return A_sp_scipy_csc
    
    # # @profile
    # def getLinearStiffnessGlobalMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler, MP):
    #     elem_set = Model.element.getElementSet()
    #     nodedof = len(elem_set["dofs"]["d"])
    #     shape_set = Model.shape.getShapeSet()
    #     nodecon = len(shape_set["nodes"])
    #     elemdof = nodecon * nodedof
    #     nodetot = coord.shape[0]
    #     sdof = nodedof * nodetot

    #     ith = []
    #     jth = []
    #     val = []

    #     @staticmethod
    #     def process_element(ee):
    #         matrix = Model.element.getStifLinearMat(
    #             Model, inci, coord, tabmat, tabgeo, intgauss, ee
    #         )
    #         loc = AssemblerFULL.__getLoc(Model, inci, ee)
            
    #         ith_local = zeros((inci.shape[0] * (elemdof * elemdof)), dtype=INT32)
    #         jth_local = zeros((inci.shape[0] * (elemdof * elemdof)), dtype=INT32)
    #         val_local = zeros((inci.shape[0] * (elemdof * elemdof)), dtype=FLT64)
            
    #         ith_local, jth_local, val_local = AssemblerFULL.__getVectorization(ith_local, jth_local, val_local, loc, matrix, ee, elemdof)
            
    #         return ith_local, jth_local, val_local

    #     with ThreadPoolExecutor() as executor:
    #         futures = [executor.submit(process_element, ee) for ee in range(inci.shape[0])]

    #     for future in futures:
    #         ith_local, jth_local, val_local = future.result()
    #         ith = append(ith, ith_local)
    #         jth = append(jth, jth_local)
    #         val = append(val, val_local)

    #     A_sp_scipy_csc = csc_matrix((val, (ith, jth)), shape=(sdof, sdof))
    #     return A_sp_scipy_csc

    def getNonLinearStiffnessGlobalMatrixAssembler():
        pass

    def getMassConsistentGlobalMatrixAssembler(Model, MP=None):
        inci = Model.inci
        coord = Model.coord
        tabmat = Model.tabmat
        tabgeo = Model.tabgeo
        intgauss = Model.intgauss

        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        elemdof = nodecon * nodedof
        nodetot = coord.shape[0]
        sdof = nodedof * nodetot

        ith = zeros((inci.shape[0] * (elemdof * elemdof)), dtype=INT32)
        jth = zeros((inci.shape[0] * (elemdof * elemdof)), dtype=INT32)
        val = zeros((inci.shape[0] * (elemdof * elemdof)), dtype=FLT64)

        for ee in range(inci.shape[0]):
            matrix = Model.element.getMassConsistentMat(
                Model, inci, coord, tabmat, tabgeo, intgauss, ee
            )
            loc = AssemblerFULLPOOL.__getLoc(Model, inci, ee)
            ith, jth, val = AssemblerFULLPOOL.__getVectorization(
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

    # https://en.wikipedia.org/wiki/Rotation_matrix
    def getRotationMatrix(node_list, coord, ndof):
        return getRotationMatrix(node_list, coord, ndof)

    
    # @staticmethod
    def getVectorization(ith, jth, val, loc, matrix, ee, elemdof):
        return getVectorization(ith, jth, val, loc, matrix, ee, elemdof)


    # @staticmethod
    def getLoc(Model, inci, element_number):
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        nodelist = Model.shape.getNodeList(inci, element_number)
        loc = Model.shape.getLocKey(nodelist, nodedof)
        return array(loc)
    
@staticmethod
def process_element(ee, Model, inci, coord, tabmat, tabgeo, intgauss, elemdof):
    matrix = Model.element.getStifLinearMat(Model, inci, coord, tabmat, tabgeo, intgauss, ee)
    loc = AssemblerFULLPOOL.getLoc(Model, inci, ee)
    ith_local = zeros((inci.shape[0] * (elemdof * elemdof)), dtype=INT32)
    jth_local = zeros((inci.shape[0] * (elemdof * elemdof)), dtype=INT32)
    val_local = zeros((inci.shape[0] * (elemdof * elemdof)), dtype=FLT64)
    ith_local, jth_local, val_local = AssemblerFULLPOOL.getVectorization(ith_local, jth_local, val_local, loc, matrix, ee, elemdof)
    return ith_local, jth_local, val_local