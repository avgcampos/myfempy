from abc import ABC, abstractmethod
from re import I

import numpy as np
import scipy.sparse as sp


# @profile
def getMatrix(
    Model, inci, coord, tabmat, tabgeo, intgauss, element_number, type_assembler
):
    if type_assembler == "linear_stiffness":
        return Model.element.getStifLinearMat(
            Model, inci, coord, tabmat, tabgeo, intgauss, element_number
        )

    elif type_assembler == "mass_consistent":
        return Model.element.getMassConsistentMat(
            Model, inci, coord, tabmat, tabgeo, intgauss, element_number
        )

    elif type_assembler == "mass_lumped":
        return Model.element.getMassLumpedMat(
            Model, inci, coord, tabmat, tabgeo, intgauss, element_number
        )


def getLoc(Model, inci, element_number):
    elem_set = Model.element.getElementSet()
    nodedof = len(elem_set["dofs"]["d"])
    nodelist = Model.shape.getNodeList(inci, element_number)
    loc = Model.shape.getLocKey(nodelist, nodedof)
    return np.array(loc)


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

    @abstractmethod
    def getDirichletNH():
        pass
