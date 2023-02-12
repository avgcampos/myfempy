#!/usr/bin/env python
from myfempy.felib.materset import get_elasticity
import numpy as np

__doc__ = """
truss21.py: Truss 2D 2-node linear Finite Element
"""


class Truss21:
    """_summary_"""

    def __init__(self, modelinfo):
        self.dofe = modelinfo["nodecon"][0] * modelinfo["nodedof"][0]
        self.fulldof = modelinfo["nodedof"][0] * len(modelinfo["coord"])
        self.nodedof = modelinfo["nodedof"][0]
        self.nelem = len(modelinfo["inci"])
        self.nnode = len(modelinfo["coord"])
        self.inci = modelinfo["inci"]
        self.coord = modelinfo["coord"]
        self.tabmat = modelinfo["tabmat"]
        self.tabgeo = modelinfo["tabgeo"]

    @staticmethod
    def elemset():
        """_summary_

        Returns:
            _description_
        """
        dofelem = {
            "key": "truss21",
            "id": 120,
            "def": "struct 1D",
            "dofs": ["ux", "uy"],
            "nnodes": ["i", "j"],
            "tensor": ["sxx"],
        }
        return dofelem

    def lockey(self, list_node):
        """_summary_

        Arguments:
            list_node -- _description_

        Returns:
            _description_
        """
        noi = list_node[0]
        noj = list_node[1]
        loc = np.array(
            [
                self.nodedof * noi - 2,
                self.nodedof * noi - 1,
                self.nodedof * noj - 2,
                self.nodedof * noj - 1,
            ]
        )
        return loc

    def stiff_linear(self, ee):
        """_summary_

        Arguments:
            ee -- _description_

        Returns:
            _description_
        """
        noi = int(self.inci[ee, 4])
        noj = int(self.inci[ee, 5])
        noix = self.coord[noi - 1, 1]
        noiy = self.coord[noi - 1, 2]
        nojx = self.coord[noj - 1, 1]
        nojy = self.coord[noj - 1, 2]
        D = get_elasticity(self.tabmat, self.inci, ee)
        E = D[0]
        A = self.tabgeo[int(self.inci[ee, 3] - 1), 0]
        L = np.sqrt((nojx - noix) ** 2 + (nojy - noiy) ** 2)
        s = (nojy - noiy) / L
        c = (nojx - noix) / L
        T = np.zeros((self.dofe, self.dofe))
        T[0, 0] = c
        T[0, 1] = s
        T[1, 0] = -s
        T[1, 1] = c
        T[2, 2] = c
        T[2, 3] = s
        T[3, 2] = -s
        T[3, 3] = c
        ket2 = np.zeros((self.dofe, self.dofe))
        ket2[0, 0] = 1.0
        ket2[0, 2] = -1.0
        ket2[2, 0] = -1.0
        ket2[2, 2] = 1.0
        ket2 = ((E * A) / L) * ket2
        ket2t = np.dot(np.dot(np.transpose(T), ket2), T)
        list_node = [noi, noj]
        loc = Truss21.lockey(self, list_node)
        return ket2t, loc

    def matrix_b(self, ee, csc):
        """_summary_

        Arguments:
            ee -- _description_
            csc -- _description_

        Returns:
            _description_
        """
        noi = int(self.inci[ee, 4])
        noj = int(self.inci[ee, 5])
        noix = self.coord[noi - 1, 1]
        noiy = self.coord[noi - 1, 2]
        nojx = self.coord[noj - 1, 1]
        nojy = self.coord[noj - 1, 2]
        D = get_elasticity(self.tabmat, self.inci, ee)
        E = D[0]
        L = np.sqrt((nojx - noix) ** 2 + (nojy - noiy) ** 2)
        s = (nojy - noiy) / L
        c = (nojx - noix) / L
        T = np.array([[c, s, 0, 0], [0, 0, c, s]])
        B = (E / L) * np.array([-1, 1])
        return B, T


if __name__ == "__main__":
    import doctest

    doctest.testmod()
