from abc import ABC, abstractmethod

import numpy as np
import scipy.sparse as sp

# def getSolver(set_solver):
#     if set_solver['solver'] == 'SLI':
#         from myfempy.core.staticlinear import StaticLinear
#         return StaticLinear
#     elif set_solver['solver'] == 'EIG':
#         from myfempy.core.modallinear import ModalLinear
#         return ModalLinear
#     else:
#         pass


class Solver(ABC):
    '''Solver API Class <ClassService>'''   
    
    @abstractmethod
    def Solve():
        pass


    def getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler):
                
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]['d'])

        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set['nodes'])

        elemdof = nodecon * nodedof

        elemtot = len(inci)
        nodetot = len(coord)

        ith = np.zeros((elemtot * (elemdof * elemdof)), dtype=int)
        jth = np.zeros((elemtot * (elemdof * elemdof)), dtype=int)
        val = np.zeros((elemtot * (elemdof * elemdof)), dtype=float)
                
        for ee in range(elemtot):
            mat, loc = Solver.__setassembler(Model, inci, coord, tabmat, tabgeo, intgauss, ee, type_assembler)
            ith[(elemdof*elemdof)*ee:(elemdof*elemdof)*(ee+1)] = np.tile(loc.reshape(1, elemdof).T, (1, elemdof)).flatten("F")
            jth[(elemdof*elemdof)*ee:(elemdof*elemdof)*(ee+1)] = np.transpose(np.tile(loc.reshape(1, elemdof).T, (1, elemdof))).flatten("F")
            val[(elemdof*elemdof)*ee:(elemdof*elemdof)*(ee+1)] = mat.flatten("F")
        
        A_sp_scipy_csc = sp.csc_matrix((val, (ith, jth)), shape=(nodedof * nodetot, nodedof * nodetot,),)        
        return A_sp_scipy_csc

    # assembly with PETSC
    # def assembleCSR(problem, dofs):
    #     problem.newton_update(dofs.reshape(
    #         (problem.num_total_nodes, problem.vec))).reshape(-1)
    #     A_sp_scipy = scipy.sparse.csr_array(
    #         (problem.V, (problem.I, problem.J)),
    #         shape=(problem.num_total_dofs, problem.num_total_dofs))

    #     A = PETSc.Mat().createAIJ(size=A_sp_scipy.shape,
    #                             csr=(A_sp_scipy.indptr, A_sp_scipy.indices,
    #                                 A_sp_scipy.data))
    #     for i in range(len(problem.node_inds_list)):
    #         row_inds = onp.array(problem.node_inds_list[i] * problem.vec +
    #                             problem.vec_inds_list[i],
    #                             dtype=onp.int32)
    #         A.zeroRows(row_inds)

    #     row, col, val = A.getValuesCSR()
    #     A_sp_scipy.data = val
    #     A_sp_scipy.indices = col
    #     A_sp_scipy.indptr = row

    #     return A_sp_scipy
    

    def getLoadAssembler(loadaply, nodetot, nodedof):
        forcevec = np.zeros((nodedof * nodetot,len(np.unique(loadaply[:, 3])),))
        
        for fstep in range(len(np.unique(loadaply[:, 3]))):
            forceaply = loadaply[np.where(loadaply[:, 3] == fstep + 1), :][0]
            nload = forceaply.shape[0]
            
            if nodedof == 1:
                for ii in range(nload):
                    if int(forceaply[ii, 1]) == 1:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof- 1))
                        forcevec[gdlload, fstep] += forceaply[ii, 2]
            
            elif nodedof == 2:
                for ii in range(nload):
                    if int(forceaply[ii, 1]) == 1:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
                        forcevec[gdlload, fstep] += forceaply[ii, 2]
                    
                    elif int(forceaply[ii, 1]) == 2:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
                        forcevec[gdlload, fstep] += forceaply[ii, 2]
                        
            elif nodedof == 3:
                for ii in range(nload):
                    if int(forceaply[ii, 1]) == 1:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof))
                        forcevec[gdlload, fstep] += forceaply[ii, 2]
                    
                    elif int(forceaply[ii, 1]) == 2:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 1))
                        forcevec[gdlload, fstep] += forceaply[ii, 2]
                        
                    elif int(forceaply[ii, 1]) == 3:
                        gdlload = int(nodedof * forceaply[ii, 0] - (nodedof - 2))
                        forcevec[gdlload, fstep] += forceaply[ii, 2]
        
            else:
                pass
        
        
        forcevec = sp.csc_matrix(forcevec)
        return forcevec
    

    def getConstrains(constrains, nodetot, nodedof):
        ntbc = len(constrains)
        fixedof = np.zeros((1, nodedof * nodetot))
        
        if nodedof == 1:
            for ii in range(ntbc):
                no = int(constrains[ii, 1])
                if int(constrains[ii, 0]) == 1:
                    fixedof[0, nodedof * no - 1] = (nodedof * no)
        
        elif nodedof == 2:
            for ii in range(ntbc):
                no = int(constrains[ii, 1])
                if int(constrains[ii, 0]) == 0:
                    fixedof[0, nodedof * no - 2] = (nodedof * no - 1)
                    fixedof[0, nodedof * no - 1] = (nodedof * no)
                elif int(constrains[ii, 0]) == 1:
                    fixedof[0, nodedof * no - 2] = (nodedof * no - 1)
                elif int(constrains[ii, 0]) == 2:
                    fixedof[0, nodedof * no - 1] = (nodedof * no)
            
        elif nodedof == 3:
            for ii in range(ntbc):
                no = int(constrains[ii, 1])
                if int(constrains[ii, 0]) == 0:
                    fixedof[0, nodedof * no - 3] = (nodedof * no - 2)
                    fixedof[0, nodedof * no - 2] = (nodedof * no - 1)
                    fixedof[0, nodedof * no - 1] = (nodedof * no)
                elif int(constrains[ii, 0]) == 1:
                    fixedof[0, nodedof * no - 3] = (nodedof * no - 2)
                elif int(constrains[ii, 0]) == 2:
                    fixedof[0, nodedof * no - 2] = (nodedof * no - 1)
                elif int(constrains[ii, 0]) == 3:
                    fixedof[0, nodedof * no - 1] = (nodedof * no)
            
        fixedof = fixedof[np.nonzero(fixedof)]
        fixedof = fixedof - np.ones_like(fixedof)
        alldof = np.arange(0, nodedof * nodetot, 1, int)
        freedof = np.setdiff1d(alldof, fixedof)
        return freedof, fixedof
    
    def addMatrix(A, A_add, loc):
        A[np.ix_(loc, loc)] += A_add
        return A
    
    def setSteps(steps):
        """steps setting"""
        start = steps["start"]
        end = steps["end"]
        substep = steps["step"]
        if (end - start) == 0:
            nsteps = int(end)
        else:
            nsteps = int((end - start) / substep)
        return nsteps
    

    def __setassembler(Model, inci, coord, tabmat, tabgeo, intgauss, element_number, type_assembler):

        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]['d'])
        
        nodelist = Model.shape.getNodeList(inci, element_number)
         
        if type_assembler == 'linear_stiffness':
            mat = Model.element.getStifLinearMat(Model, inci, coord, tabmat, tabgeo, intgauss, element_number)
            loc = Model.shape.getShapeKey(nodelist, nodedof)
            return mat, loc
        
        elif type_assembler == 'mass_consistent':
            mat = Model.element.getMassConsistentMat(Model, inci, coord, tabmat, tabgeo, intgauss, element_number)
            loc = Model.shape.getShapeKey(nodelist, nodedof)
            return mat, loc
        
        elif type_assembler == 'mass_lumped':
            mat = Model.element.getMassLumpedMat(Model, inci, coord, tabmat, tabgeo, intgauss, element_number)
            loc = Model.shape.getShapeKey(nodelist, nodedof)
            return mat, loc