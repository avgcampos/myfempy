import numpy as np

from myfempy.core.mesh.mesh import Mesh


class LegacyQuad4(Mesh):
    """Mesh Quad Class <ConcreteClassService>"""

    def getElementConection(set_mesh):
        """get a quadrangular 4 nodes mesh

        (l)---------------(k)
         |                 |
         |                 |
         |      {1}        |
         |                 |
         |                 |
        (i)---------------(j)

        """

        nelx = set_mesh["NX"]
        nely = set_mesh["NY"]

        nnx = nelx + 1
        nel = nelx * nely

        conec = np.zeros((nel, 5), dtype=np.int64)
        for i in range(1, nely + 1, 1):
            conec[np.arange((i - 1) * nelx, i * nelx, 1).tolist(), 0] = np.arange(
                (i - 1) * nelx + 1, i * nelx + 1, 1
            ).tolist()
            conec[np.arange((i - 1) * nelx, i * nelx, 1).tolist(), 1] = np.arange(
                (i - 1) * nnx + 1, i * nnx, 1
            ).tolist()
            conec[np.arange((i - 1) * nelx, i * nelx, 1).tolist(), 2] = np.arange(
                (i - 1) * nnx + 2, i * nnx + 1, 1
            ).tolist()
            conec[np.arange((i - 1) * nelx, i * nelx, 1).tolist(), 3] = np.arange(
                i * nnx + 2, (i + 1) * nnx + 1, 1
            ).tolist()
            conec[np.arange((i - 1) * nelx, i * nelx, 1).tolist(), 4] = np.arange(
                i * nnx + 1, (i + 1) * nnx, 1
            ).tolist()

        return conec

    def getNodesCoord(set_mesh):
        nelx = set_mesh["NX"]
        nely = set_mesh["NY"]
        lx = set_mesh["LX"]
        ly = set_mesh["LY"]

        nnx = nelx + 1
        nny = nely + 1
        nos = nnx * nny

        dx = np.arange(0, lx + lx / nelx, lx / nelx)
        dy = np.arange(0, ly + ly / nely, ly / nely)
        xv, yv = np.meshgrid(dx, dy)

        coord = np.zeros((nos, 4), dtype=np.float64)
        for i in range(1, nny + 1, 1):
            coord[np.arange((i - 1) * nnx, i * nnx, 1).tolist(), 0] = np.arange(
                (i - 1) * nnx + 1, i * nnx + 1, 1
            ).tolist()
            coord[np.arange((i - 1) * nnx, i * nnx, 1).tolist(), 1] = xv[i - 1, :]
            coord[np.arange((i - 1) * nnx, i * nnx, 1).tolist(), 2] = yv[i - 1, :]

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
