from abc import ABC, abstractmethod
from re import I

import numpy as np
import scipy.sparse as sp

def setAssembler(Model, inci, coord, tabmat, tabgeo, intgauss, element_number, type_assembler):

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

class Assembler(ABC):
    
    """
     Assembler API Class <ClassService>
    """
    
    @abstractmethod
    def getMatrixAssembler():
        pass
    
    @abstractmethod
    def getLoadAssembler():
        pass
    
    @abstractmethod
    def getConstrains():
        pass
    
