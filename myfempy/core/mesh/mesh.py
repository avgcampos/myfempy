import os
from abc import ABC, abstractmethod

import numpy as np

from myfempy.io.iogmsh import get_gmsh_geo, get_gmsh_msh


def setMesh(set_mesh):
    if set_mesh["TYPE"] == "add":
        return MeshADD

    elif set_mesh["TYPE"] == "legacy":
        if set_mesh["shape"] == "quad4":
            from myfempy.core.mesh.legacyquad4 import LegacyQuad4

            return LegacyQuad4
        elif set_mesh["shape"] == "tria3":
            from myfempy.core.mesh.legacytria3 import LegacyTria3

            return LegacyTria3
        else:
            pass

    elif set_mesh["TYPE"] == "gmsh":
        if "meshimport" in set_mesh.keys():
            pass
        else:
            filename = set_mesh["user_path"] + "/" + set_mesh["filename"]
            get_gmsh_geo(filename, set_mesh)
            get_gmsh_msh(filename, set_mesh)
        from myfempy.core.mesh.gmsh import MeshGmsh

        return MeshGmsh
    else:
        pass


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
