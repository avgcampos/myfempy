#!/usr/bin/env python
from myfempy.felib.materset import get_elasticity
import numpy as np

__doc__ = """
plane31.py: Triagular Plane 3-node linear Finite Element
"""


class Plane31:
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
        self.ntensor = modelinfo["ntensor"][0]

    @staticmethod
    def elemset():
        dofelem = {
            "key": "plane31",
            "id": 210,
            "def": "struct 2D",
            "dofs": ["ux", "uy"],
            "nnodes": ["i", "j", "k"],
            "tensor": ["sxx", "syy", "sxy"],
        }
        return dofelem

    def lockey(self, nodelist):
        noi = nodelist[0]
        noj = nodelist[1]
        nok = nodelist[2]
        loc = np.array(
            [
                self.nodedof * noi - 2,
                self.nodedof * noi - 1,
                self.nodedof * noj - 2,
                self.nodedof * noj - 1,
                self.nodedof * nok - 2,
                self.nodedof * nok - 1,
            ]
        )
        return loc

    def matriz_b(self, nodelist, intpl):
        noi = nodelist[0]
        noj = nodelist[1]
        nok = nodelist[2]
        xi = self.coord[noi - 1, 1]
        yi = self.coord[noi - 1, 2]
        xj = self.coord[noj - 1, 1]
        yj = self.coord[noj - 1, 2]
        xk = self.coord[nok - 1, 1]
        yk = self.coord[nok - 1, 2]
        C = np.array([[1, xi, yi], [1, xj, yj], [1, xk, yk]])
        A = (1 / 2) * abs(np.linalg.det(C))
        beti = yj - yk
        betj = yk - yi
        betk = yi - yj
        gami = xk - xj
        gamj = xi - xk
        gamk = xj - xi
        B = np.zeros((self.ntensor, self.dofe))
        B[0, 0] = beti
        B[0, 2] = betj
        B[0, 4] = betk
        B[1, 1] = gami
        B[1, 3] = gamj
        B[1, 5] = gamk
        B[2, 0] = gami
        B[2, 1] = beti
        B[2, 2] = gamj
        B[2, 3] = betj
        B[2, 4] = gamk
        B[2, 5] = betk
        B = (1 / (2 * A)) * B
        return B, A

    def stiff_linear(self, ee):
        noi = int(self.inci[ee, 4])
        noj = int(self.inci[ee, 5])
        nok = int(self.inci[ee, 6])
        nodelist = [noi, noj, nok]
        D = get_elasticity(self.tabmat, self.inci, ee)
        L = self.tabgeo[int(self.inci[ee, 3] - 1), 4]
        intpl = 0.0
        B, A = Plane31.matriz_b(self, nodelist, intpl)
        ket3 = L * A * np.dot(np.dot(np.transpose(B), D), B)
        loc = Plane31.lockey(self, nodelist)
        return ket3, loc

    def mass(self, ee):
        noi = int(self.inci[ee, 4])
        noj = int(self.inci[ee, 5])
        nok = int(self.inci[ee, 6])
        xi = self.coord[noi - 1, 1]
        yi = self.coord[noi - 1, 2]
        xj = self.coord[noj - 1, 1]
        yj = self.coord[noj - 1, 2]
        xk = self.coord[nok - 1, 1]
        yk = self.coord[nok - 1, 2]
        C = np.array([[1, xi, yi], [1, xj, yj], [1, xk, yk]])
        A = (1 / 2) * abs(np.linalg.det(C))
        L = self.tabgeo[int(self.inci[ee, 3] - 1), 4]
        R = self.tabmat[int(self.inci[ee, 2] - 1), 6]
        mat_aux = 2 * np.eye(self.dofe)
        mat_aux[0, 2] = 1
        mat_aux[0, 4] = 1
        mat_aux[1, 3] = 1
        mat_aux[1, 5] = 1
        mat_aux[2, 0] = 1
        mat_aux[2, 4] = 1
        mat_aux[3, 1] = 1
        mat_aux[3, 5] = 1
        mat_aux[4, 0] = 1
        mat_aux[4, 2] = 1
        mat_aux[5, 1] = 1
        mat_aux[5, 3] = 1
        met3 = (R * A * L / 12) * mat_aux
        list_node = [noi, noj, nok]
        loc = Plane31.lockey(self, list_node)
        return met3, loc


if __name__ == "__main__":
    import doctest

    doctest.testmod()
