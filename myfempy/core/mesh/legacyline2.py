import numpy as np

from myfempy.core.mesh.mesh import Mesh


class LegacyLine2(Mesh):
    """Mesh Line Class <ConcreteClassService>"""

    def getElementConection(set_mesh):
        """get a linear 2 nodes mesh

        (i)------{1}------(j)

        """

        nel = set_mesh["NX"]

        conec = np.zeros((nel, 3), dtype=np.int64)
        conec[:, 0] = np.arange(1, nel + 1, 1).tolist()
        conec[:, 1] = np.arange(1, nel + 1, 1).tolist()
        conec[:, 2] = np.arange(2, nel + 2, 1).tolist()
        return conec

    def getNodesCoord(set_mesh):
        nel = set_mesh["NX"]
        lx = set_mesh["LX"]

        step = lx / nel

        coord = np.zeros((nel + 1, 4), dtype=np.float64)
        coord[:, 0] = np.arange(1, nel + 2, 1).tolist()
        coord[:, 1] = np.arange(0, lx + step, step).tolist()
        return coord

    def getElementList(conec, meshset, modeldata):
        elemlist = [[None] * 3]
        for ee in range(len(conec)):
            elemlist.append(
                [
                    int(conec[ee, 0]),
                    meshset,
                    modeldata["MATERIAL"]["PROPMAT"][0]["NAME"],
                    modeldata["GEOMETRY"]["PROPGEO"][0]["NAME"],
                    conec[ee, 1:].astype(int).tolist(),
                ]
            )
        elemlist = elemlist[1::][::]
        return elemlist
