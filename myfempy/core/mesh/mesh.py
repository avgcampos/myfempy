import os
from abc import ABC, abstractmethod

import numpy as np


class Mesh(ABC):
    """Element API Class <ClassService>"""

    @abstractmethod
    def getElementConection():
        pass

    @abstractmethod
    def getNodesCoord():
        pass

    @abstractmethod
    def getElementList():
        pass


class MeshADD(Mesh):
    """Mesh Add Class <ConcreteClassService>"""

    def getNodesCoord(set_mesh):
        return np.asarray(set_mesh["COORD"])

    def getElementConection(set_mesh):
        return np.asarray(set_mesh["INCI"])

    def getElementList(conec, meshset, modeldata):
        elemlist = [[None] * 3]
        for ee in range(len(conec)):
            elemlist.append(
                [
                    int(conec[ee, 0]),
                    meshset,
                    modeldata["MATERIAL"]["PROPMAT"][conec[ee, 1] - 1]["NAME"],
                    modeldata["GEOMETRY"]["PROPGEO"][conec[ee, 2] - 1]["NAME"],
                    conec[ee, 3:].astype(int).tolist(),
                ]
            )
        elemlist = elemlist[1::][::]
        return elemlist
