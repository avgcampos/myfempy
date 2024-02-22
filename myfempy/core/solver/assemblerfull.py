from __future__ import annotations

import numpy as np
from scipy import empty, sparse

# from myfempy.core.alglin import linsolve_spsolve
from myfempy.core.solver.assembler import Assembler, setAssembler
from myfempy.expe.asmb_cython.import_assembler_cython2py import getMatrixAssemblerSYMM

class AssemblerFULL(Assembler):

    """
     Static Linear Solver Class <ConcreteClassService>
    """
    
    def getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler):
        """
        getMatrixAssembler Assembler module <ConcreteClassService>

        Returns:
            matrix sparse
        """    
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
            mat, loc = setAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, ee, type_assembler)
            ith[(elemdof*elemdof)*ee:(elemdof*elemdof)*(ee+1)] = np.tile(loc.reshape(1, elemdof).T, (1, elemdof)).flatten("F")
            jth[(elemdof*elemdof)*ee:(elemdof*elemdof)*(ee+1)] = np.transpose(np.tile(loc.reshape(1, elemdof).T, (1, elemdof))).flatten("F")
            val[(elemdof*elemdof)*ee:(elemdof*elemdof)*(ee+1)] = mat.flatten("F")       
        A_sp_scipy_csc = sparse.csc_matrix((val, (ith, jth)), shape=(nodedof * nodetot, nodedof * nodetot),)        
        return A_sp_scipy_csc
    
    
    def getLoadAssembler(loadaply, nodetot, nodedof):
        
        """
        getLoadAssembler Assembler module <ConcreteClassService>

        Returns:
            vector sparse
        """
        
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
        forcevec = sparse.csc_matrix(forcevec)
        return forcevec
    

    def getConstrains(constrains, nodetot, nodedof):
        
        """
        getConstrains Constrain module <ConcreteClassService>

        Returns:
            _description_
        """
              
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