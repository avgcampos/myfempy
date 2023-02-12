#!/usr/bin/env python
from myfempy.felib.quadrature import gaussian, no_interpol
from myfempy.felib.materset import get_elasticity
import numpy as np

__doc__ = """
plate41.py: Quatrangular Isoparametric Plate Mindlin 4-node linear Finite Element
"""


class Plate41:
    """_summary_"""

    def __init__(self, modelinfo):
        self.dofe = modelinfo["nodecon"][0] * modelinfo["nodedof"][0]
        self.nodecon = modelinfo["nodecon"][0]
        self.fulldof = modelinfo["nodedof"][0] * len(modelinfo["coord"])
        self.nodedof = modelinfo["nodedof"][0]
        self.nelem = len(modelinfo["inci"])
        self.nnode = len(modelinfo["coord"])
        self.inci = modelinfo["inci"]
        self.coord = modelinfo["coord"]
        self.tabmat = modelinfo["tabmat"]
        self.tabgeo = modelinfo["tabgeo"]
        self.ntensor = modelinfo["ntensor"][0]
        if modelinfo["quadra"][0] == 1:
            self.npp = modelinfo["quadra"][1]
            self.quadra = gaussian(self.npp)
        elif modelinfo["quadra"][0] == 0:
            self.npp = modelinfo["quadra"][1]
            self.quadra = no_interpol(self.npp)

    @staticmethod
    def elemset():
        """_summary_

        Returns:
            _description_
        """
        dofelem = {
            "key": "plate41",
            "id": 230,
            "def": "struct 2D",
            "dofs": ["uz", "rx", "ry"],
            "nnodes": ["i", "j", "k", "l"],
            "tensor": ["sxx", "syy", "sxy"],
        }
        return dofelem

    def lockey(self, nodelist):
        """_summary_

        Arguments:
            nodelist -- _description_

        Returns:
            _description_
        """
        noi = nodelist[0]
        noj = nodelist[1]
        nok = nodelist[2]
        nol = nodelist[3]
        loc = np.array(
            [
                self.nodedof * noi - 2,
                self.nodedof * noi - 1,
                self.nodedof * noj - 2,
                self.nodedof * noj - 1,
                self.nodedof * nok - 2,
                self.nodedof * nok - 1,
                self.nodedof * nol - 2,
                self.nodedof * nol - 1,
            ]
        )
        return loc

    def matriz_b(self, nodelist, intpl):
        """_summary_

        Arguments:
            nodelist -- _description_
            intpl -- _description_

        Returns:
            _description_
        """
        noi = nodelist[0]
        noj = nodelist[1]
        nok = nodelist[2]
        nol = nodelist[3]
        xi = self.coord[noi - 1, 1]
        yi = self.coord[noi - 1, 2]
        xj = self.coord[noj - 1, 1]
        yj = self.coord[noj - 1, 2]
        xk = self.coord[nok - 1, 1]
        yk = self.coord[nok - 1, 2]
        xl = self.coord[nol - 1, 1]
        yl = self.coord[nol - 1, 2]
        matXY = np.array([[xi, yi], [xj, yj], [xk, yk], [xl, yl]])
        x = intpl[0]
        y = intpl[1]
        N1x = -1 + y
        N2x = 1 - y
        N3x = 1 + y
        N4x = -1 - y
        N1y = -1 + x
        N2y = -1 - x
        N3y = 1 + x
        N4y = 1 - x
        dN = (1 / 4) * np.array([[N1x, N2x, N3x, N4x], [N1y, N2y, N3y, N4y]])
        J = np.dot(dN, matXY)
        detJ = J[0, 0] * J[1, 1] - J[1, 0] * J[0, 1]
        invJ = (1 / detJ) * np.array(
            [
                [J[1, 1], -J[0, 1], 0.0, 0.0],
                [0.0, 0.0, -J[1, 0], J[0, 0]],
                [-J[1, 0], J[0, 0], J[1, 1], -J[0, 1]],
            ]
        )
        pN = np.zeros((self.nodecon, self.dofe))
        pN[0, 0] = N1x
        pN[0, 2] = N2x
        pN[0, 4] = N3x
        pN[0, 6] = N4x
        pN[1, 0] = N1y
        pN[1, 2] = N2y
        pN[1, 4] = N3y
        pN[1, 6] = N4y
        pN[2, 1] = N1x
        pN[2, 3] = N2x
        pN[2, 5] = N3x
        pN[2, 7] = N4x
        pN[3, 1] = N1y
        pN[3, 3] = N2y
        pN[3, 5] = N3y
        pN[3, 7] = N4y
        pN = (1 / 4) * pN
        B = np.dot(invJ, pN)
        return B, detJ

    def stiff_linear(self, ee):
        """_summary_

        Arguments:
            ee -- _description_

        Returns:
            _description_
        """
        noi = int(self.inci[ee, 4])
        noj = int(self.inci[ee, 5])
        nok = int(self.inci[ee, 6])
        nol = int(self.inci[ee, 7])
        nodelist = [noi, noj, nok, nol]
        D = get_elasticity(self.tabmat, self.inci, ee)
        L = self.tabgeo[int(self.inci[ee, 3] - 1), 4]
        xp, wp = self.quadra
        xpp = xp[0, :]
        ypp = xp[1, :]
        A = 0.0
        keq4 = np.zeros((self.dofe, self.dofe))
        for pp in range(0, self.npp):
            intpl = [xpp[pp], ypp[pp]]
            Bpp, detJ = Plate41.matriz_b(self, nodelist, intpl)
            A += detJ
            keq4 += (
                np.dot(np.dot(np.transpose(Bpp), D), Bpp) * L * detJ * wp[pp] * wp[pp]
            )
        loc = Plate41.lockey(self, nodelist)
        return keq4, loc

    def mass(self, ee):
        """_summary_

        Arguments:
            ee -- _description_

        Returns:
            _description_
        """
        noi = int(self.inci[ee, 4])
        noj = int(self.inci[ee, 5])
        nok = int(self.inci[ee, 6])
        nol = int(self.inci[ee, 7])
        xi = self.coord[noi - 1, 1]
        yi = self.coord[noi - 1, 2]
        xj = self.coord[noj - 1, 1]
        yj = self.coord[noj - 1, 2]
        xk = self.coord[nok - 1, 1]
        yk = self.coord[nok - 1, 2]
        xl = self.coord[nol - 1, 1]
        yl = self.coord[nol - 1, 2]
        matXY = np.array([[xi, yi], [xj, yj], [xk, yk], [xl, yl]])
        R = self.tabmat[int(self.inci[ee, 2] - 1), 6]
        L = self.tabgeo[int(self.inci[ee, 3] - 1), 4]
        xp, wp = self.quadra
        x = xp[0, :]
        y = xp[1, :]
        meq4 = np.zeros((self.dofe, self.dofe))
        for pp in range(0, self.npp):
            N1 = (1 - x[pp]) * (1 - y[pp])
            N2 = (1 + x[pp]) * (1 - y[pp])
            N3 = (1 + x[pp]) * (1 + y[pp])
            N4 = (1 - x[pp]) * (1 + y[pp])
            N = (1 / 4) * np.array(
                [[N1, 0, N2, 0, N3, 0, N4, 0], [0, N1, 0, N2, 0, N3, 0, N4]]
            )
            N1x = -(1 - y[pp])
            N2x = 1 - y[pp]
            N3x = 1 + y[pp]
            N4x = -(1 + y[pp])
            N1y = -(1 - x[pp])
            N2y = -(1 + x[pp])
            N3y = 1 + x[pp]
            N4y = 1 - x[pp]
            dN = (1 / 4) * np.array([[N1x, N2x, N3x, N4x], [N1y, N2y, N3y, N4y]])
            J = np.dot(dN, matXY)
            detJ = np.linalg.det(J)
            meq4 += np.dot(np.transpose(N), N) * R * L * detJ * wp[pp] * wp[pp]
        list_node = [noi, noj, nok, nol]
        loc = Plate41.lockey(self, list_node)
        return meq4, loc


if __name__ == "__main__":
    import doctest

    doctest.testmod()
