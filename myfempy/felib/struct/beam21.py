#!/usr/bin/env python
from myfempy.felib.materset import get_elasticity
import numpy as np

__doc__ = """
beam21.py: Beam 1D 2-node linear Finite Element
"""


class Beam21:
    """ """

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
        dofelem = {
            "key": "beam21",
            "id": 130,
            "def": "struct 1D",
            "dofs": ["uy", "rz"],
            "nnodes": ["i", "j"],
            "tensor": ["sxx"],
        }
        return dofelem

    def lockey(self, list_node):
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
        noi = int(self.inci[ee, 4])
        noj = int(self.inci[ee, 5])
        noix = self.coord[noi - 1, 1]
        noiy = self.coord[noi - 1, 2]
        nojx = self.coord[noj - 1, 1]
        nojy = self.coord[noj - 1, 2]
        D = get_elasticity(self.tabmat, self.inci, ee)
        E = D[0]
        I33 = self.tabgeo[int(self.inci[ee, 3] - 1), 1]
        L = np.sqrt((nojx - noix) ** 2 + (nojy - noiy) ** 2)
        keb1 = np.zeros((self.dofe, self.dofe))
        keb1[0, 0] = 12
        keb1[0, 1] = 6 * L
        keb1[0, 2] = -12
        keb1[0, 3] = 6 * L
        keb1[1, 0] = 6 * L
        keb1[1, 1] = 4 * L**2
        keb1[1, 2] = -6 * L
        keb1[1, 3] = 2 * L**2
        keb1[2, 0] = -12
        keb1[2, 1] = -6 * L
        keb1[2, 2] = 12
        keb1[2, 3] = -6 * L
        keb1[3, 0] = 6 * L
        keb1[3, 1] = 2 * L**2
        keb1[3, 2] = -6 * L
        keb1[3, 3] = 4 * L**2
        keb1 = (E * I33 / L**3) * keb1
        list_node = [noi, noj]
        loc = Beam21.lockey(self, list_node)
        return keb1, loc

    def mass(self, ee):
        noi = int(self.inci[ee, 4])
        noj = int(self.inci[ee, 5])
        noix = self.coord[noi - 1, 1]
        noiy = self.coord[noi - 1, 2]
        nojx = self.coord[noj - 1, 1]
        nojy = self.coord[noj - 1, 2]
        R = self.tabmat[int(self.inci[ee, 2] - 1), 6]
        A = self.tabgeo[int(self.inci[ee, 3] - 1), 0]
        L = np.sqrt((nojx - noix) ** 2 + (nojy - noiy) ** 2)
        meb2 = np.zeros((self.dofe, self.dofe))
        meb2[0, 0] = 156
        meb2[0, 1] = 22 * L
        meb2[0, 2] = 54
        meb2[0, 3] = -13 * L
        meb2[1, 0] = 22 * L
        meb2[1, 1] = 4 * L**2
        meb2[1, 2] = 13 * L
        meb2[1, 3] = -3 * L**2
        meb2[2, 0] = 54
        meb2[2, 1] = 13 * L
        meb2[2, 2] = 156
        meb2[2, 3] = -22 * L
        meb2[3, 0] = -13 * L
        meb2[3, 1] = -3 * L**2
        meb2[3, 2] = -22 * L
        meb2[3, 3] = 4 * L**2
        meb2 = (R * A * L / 420) * meb2
        list_node = [noi, noj]
        loc = Beam21.lockey(self, list_node)
        return meb2, loc

    def intforces(self, U, lines):
        Fint = np.zeros((self.fulldof))
        Vy = np.zeros((len(lines[0][1]), len(lines)), dtype=float)
        Mz = np.zeros((len(lines[0][1]), len(lines)), dtype=float)
        domL = np.zeros((len(lines[0][1]), len(lines)), dtype=float)
        for ee in range(self.nelem):
            keb1, loc = Beam21.stiff_linear(self, ee)
            F = np.dot(keb1, U[loc])
            Fint[loc] = [F[0], -F[1], -F[2], F[3]]
        for ed in range(len(lines)):
            nnodes = lines[ed][1]
            for nn in range(0, len(nnodes)):
                node = int(nnodes[nn])
                Vy[nn, ed] = Fint[self.nodedof * node - 2]
                Mz[nn, ed] = Fint[self.nodedof * node - 1]
                domL[nn, ed] = self.coord[node - 1, 1]
            domLIdc = np.argsort(domL[:, ed], axis=0)
            domL[:, ed] = domL[domLIdc, ed]
            Mz[:, ed] = Mz[domLIdc, ed]
            Vy[:, ed] = Vy[domLIdc, ed]
        ifb = {"le": domL, "val": [Vy, Mz]}
        title = ["VY", "MZ"]
        return ifb, title

    def matrix_B(self, ee, csc):
        y = csc[0]
        noi = int(self.inci[ee, 4])
        noj = int(self.inci[ee, 5])
        noix = self.coord[noi - 1, 1]
        noiy = self.coord[noi - 1, 2]
        nojx = self.coord[noj - 1, 1]
        nojy = self.coord[noj - 1, 2]
        D = get_elasticity(self.tabmat, self.inci, ee)
        E = D[0]
        x_mid = (nojx - noix) / 2
        L = np.sqrt((nojx - noix) ** 2 + (nojy - noiy) ** 2)
        T = np.eye(4, 4)
        B = -(y * E) * np.array(
            [
                (12 * x_mid) / (L**3) - 6 / (L**2),
                (6 * x_mid) / (L**2) - 4 / L,
                -(12 * x_mid) / (L**3) + 6 / (L**2),
                (6 * x_mid) / (L**2) - 2 / L,
            ]
        )
        return B, T


if __name__ == "__main__":
    import doctest

    doctest.testmod()
