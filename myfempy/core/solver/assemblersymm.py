from __future__ import annotations

from numpy import zeros, float64
from scipy.sparse import coo_matrix, triu

from myfempy.core.solver.assembler import Assembler
from myfempy.core.solver.assemblerfull import AssemblerFULL
# from myfempy.core.solver.assemblersymm_numpy_v1 import getMatrixAssemblerSymm
from myfempy.core.solver.assemblersymm_cython_v4 import getMatrixAssemblerSymm

class AssemblerSYMM(Assembler):

    """
     Assembler Symmetric Banded System Class <ConcreteClassService>
    """
    # @profile
    def getMatrixAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, type_assembler, MP):
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]['d'])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set['nodes'])
        elemdof = nodecon * nodedof
        nodetot = coord.shape[0]
        sdof = nodedof * nodetot
        
        # if MP:
        #     rowsd, colsd, datad, rowsb, colsb, datab = getMatrixAssemblerSym_cy_parpool(Model, inci, coord, tabmat, tabgeo, elemdof,  intgauss, type_assembler)
        # else:
        rowsd, datad, rowsb, colsb, datab = getMatrixAssemblerSymm(Model, inci, coord, tabmat, tabgeo, elemdof,  intgauss, type_assembler)
            
        mtKG_sp_sym = coo_matrix((datab, (rowsb, colsb)), shape=(sdof, sdof)).tocsr()
        mtKG_sp_sym += mtKG_sp_sym.transpose()
        mtKG_sp_sym += coo_matrix((datad, (rowsd, rowsd)), shape=(sdof, sdof)).tocsr()
           
        return mtKG_sp_sym
    
    def getLoadAssembler(loadaply, nodetot, nodedof):
        return AssemblerFULL.getLoadAssembler(loadaply, nodetot, nodedof)
    
    def getConstrains(constrains, nodetot, nodedof):
        return AssemblerFULL.getConstrains(constrains, nodetot, nodedof)