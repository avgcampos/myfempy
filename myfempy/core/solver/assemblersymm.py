from __future__ import annotations

from numpy import zeros, float64
from scipy import sparse

# from myfempy.core.alglin import linsolve_spsolve
from myfempy.core.solver.assembler import Assembler, setAssembler
from myfempy.core.solver.assemblerfull import AssemblerFULL
# from myfempy.expe.asmb_cython.import_assembler_cython2py import getMatrixAssemblerSYMM

class AssemblerSYMM(Assembler):

    """
     Static Linear Solver Class <ConcreteClassService>
    """
    
    def getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler):
        """
        getMatrixAssembler Assembler Symmertric Matrix System module <ConcreteClassService>

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

        sdof = nodedof * nodetot

        ith_band = []
        jth_band = []
        val_band = []
        ith_diag = []
        jth_diag = []
        val_diag = []
        
        for ee in range(elemtot):
            mat, loc = setAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, ee, type_assembler)
            for ii in range(elemdof):
                KI = loc[ii]         
                for jj in range(elemdof):
                    KJ = loc[jj]
                    if KI==KJ:
                        ith_diag.append(KI)
                        jth_diag.append(KJ)
                        val_diag.append(mat[ii, jj])
                    elif KI<KJ:
                        ith_band.append(KI)
                        jth_band.append(KJ)
                        val_band.append(mat[ii, jj])
        A_sp_scipy = sparse.coo_matrix((val_band, (ith_band, jth_band)), shape=(sdof, sdof))
        A_sp_scipy += A_sp_scipy.transpose()
        A_sp_scipy += sparse.coo_matrix((val_diag, (ith_diag, jth_diag)), shape=(sdof, sdof))
        return A_sp_scipy
    
    def getLoadAssembler(loadaply, nodetot, nodedof):
        return AssemblerFULL.getLoadAssembler(loadaply, nodetot, nodedof)
    
    def getConstrains(constrains, nodetot, nodedof):
        return AssemblerFULL.getConstrains(constrains, nodetot, nodedof)