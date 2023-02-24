#!/usr/bin/env python
from myfempy.felib.quadrature import gaussian, no_interpol
from myfempy.felib.materset import get_elasticity
import numpy as np

__doc__ = """
solid81.py: Hexahedron Isoparametric Solid 8-node linear Finite Element
"""


class Solid81:
    """class Hexahedron Isoparametric Solid 8-node linear Finite Element"""

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
            
        """
        Arguments:
           modelinfo:dict     -- F.E. model dict with full information needed

        Parameters:
            dofe              -- element dof
            fulldof           -- total dof of model
            nodedof           -- node dof 
            nelem             -- total number of elements in mesh
            nnode             -- number of degree of freedom per node
            inci              -- elements conection and prop. list
            coord             -- nodes coordinates list in mesh
            tabmat            -- table of material prop.
            tabgeo            -- table of geometry prop.
            ntensor           -- dim. of tensor (stress-strain relat.)
            quadra            -- quadrature integration, only used in isoparametric elements
            npp               -- number of points to integrations
        
        """ 

    @staticmethod
    def elemset():
        """element setting"""
        
        dofelem = {
            "key": "solid81",
            "id": 320,
            "def": "struct 3D",
            "dofs": ["ux", "uy", "uz"],
            "nnodes": ["i", "j", "k", "l", "m", "n", "o", "p"],
            "tensor": ["sxx", "syy", "szz", "sxy", "syz", "szx"],
        }
        return dofelem

    def lockey(self, nodelist):
        """element lockey(dof)"""
        
        noi = nodelist[0]
        noj = nodelist[1]
        nok = nodelist[2]
        nol = nodelist[3]
        nom = nodelist[4]
        non = nodelist[5]
        noo = nodelist[6]
        nop = nodelist[7]
        loc = np.array(
            [
                self.nodedof * noi - 3,
                self.nodedof * noi - 2,
                self.nodedof * noi - 1,
                self.nodedof * noj - 3,
                self.nodedof * noj - 2,
                self.nodedof * noj - 1,
                self.nodedof * nok - 3,
                self.nodedof * nok - 2,
                self.nodedof * nok - 1,
                self.nodedof * nol - 3,
                self.nodedof * nol - 2,
                self.nodedof * nol - 1,
                self.nodedof * nom - 3,
                self.nodedof * nom - 2,
                self.nodedof * nom - 1,
                self.nodedof * non - 3,
                self.nodedof * non - 2,
                self.nodedof * non - 1,
                self.nodedof * noo - 3,
                self.nodedof * noo - 2,
                self.nodedof * noo - 1,
                self.nodedof * nop - 3,
                self.nodedof * nop - 2,
                self.nodedof * nop - 1,
            ]
        )
        return loc

    def matriz_b(self, nodelist, intpl):
        """shape function derivatives"""
        
        noi = nodelist[0]
        noj = nodelist[1]
        nok = nodelist[2]
        nol = nodelist[3]
        nom = nodelist[4]
        non = nodelist[5]
        noo = nodelist[6]
        nop = nodelist[7]
        xi = self.coord[noi - 1, 1]
        yi = self.coord[noi - 1, 2]
        zi = self.coord[noi - 1, 3]
        xj = self.coord[noj - 1, 1]
        yj = self.coord[noj - 1, 2]
        zj = self.coord[noj - 1, 3]
        xk = self.coord[nok - 1, 1]
        yk = self.coord[nok - 1, 2]
        zk = self.coord[nok - 1, 3]
        xl = self.coord[nol - 1, 1]
        yl = self.coord[nol - 1, 2]
        zl = self.coord[nol - 1, 3]
        xm = self.coord[nom - 1, 1]
        ym = self.coord[nom - 1, 2]
        zm = self.coord[nom - 1, 3]
        xn = self.coord[non - 1, 1]
        yn = self.coord[non - 1, 2]
        zn = self.coord[non - 1, 3]
        xo = self.coord[noo - 1, 1]
        yo = self.coord[noo - 1, 2]
        zo = self.coord[noo - 1, 3]
        xp = self.coord[nop - 1, 1]
        yp = self.coord[nop - 1, 2]
        zp = self.coord[nop - 1, 3]
        matXYZ = np.array(
            [
                [xi, yi, zi],
                [xj, yj, zj],
                [xk, yk, zk],
                [xl, yl, zl],
                [xm, ym, zm],
                [xn, yn, zn],
                [xo, yo, zo],
                [xp, yp, zp],
            ]
        )
        x = intpl[0]
        y = intpl[1]
        z = intpl[2]
        N1x = -(1 - y) * (1 - z)
        N2x = (1 - y) * (1 - z)
        N3x = (1 + y) * (1 - z)
        N4x = -(1 + y) * (1 - z)
        N5x = -(1 - y) * (1 + z)
        N6x = (1 - y) * (1 + z)
        N7x = (1 + y) * (1 + z)
        N8x = -(1 + y) * (1 + z)
        N1y = -(1 - x) * (1 - z)
        N2y = -(1 + x) * (1 - z)
        N3y = (1 + x) * (1 - z)
        N4y = (1 - x) * (1 - z)
        N5y = -(1 - x) * (1 + z)
        N6y = -(1 + x) * (1 + z)
        N7y = (1 + x) * (1 + z)
        N8y = (1 - x) * (1 + z)
        N1z = -(1 - x) * (1 - y)
        N2z = -(1 + x) * (1 - y)
        N3z = -(1 + x) * (1 + y)
        N4z = -(1 - x) * (1 + y)
        N5z = (1 - x) * (1 - y)
        N6z = (1 + x) * (1 - y)
        N7z = (1 + x) * (1 + y)
        N8z = (1 - x) * (1 + y)
        dN = (1 / 8) * np.array(
            [
                [N1x, N2x, N3x, N4x, N5x, N6x, N7x, N8x],
                [N1y, N2y, N3y, N4y, N5y, N6y, N7y, N8y],
                [N1z, N2z, N3z, N4z, N5z, N6z, N7z, N8z],
            ]
        )
        J = np.dot(dN, matXYZ)
        detJ = np.linalg.det(J)
        invJ = np.linalg.inv(J)
        dN1 = np.array([[N1x], [N1y], [N1z]])
        dN2 = np.array([[N2x], [N2y], [N2z]])
        dN3 = np.array([[N3x], [N3y], [N3z]])
        dN4 = np.array([[N4x], [N4y], [N4z]])
        dN5 = np.array([[N5x], [N5y], [N5z]])
        dN6 = np.array([[N6x], [N6y], [N6z]])
        dN7 = np.array([[N7x], [N7y], [N7z]])
        dN8 = np.array([[N8x], [N8y], [N8z]])
        pN1 = np.dot(invJ, dN1)
        pN2 = np.dot(invJ, dN2)
        pN3 = np.dot(invJ, dN3)
        pN4 = np.dot(invJ, dN4)
        pN5 = np.dot(invJ, dN5)
        pN6 = np.dot(invJ, dN6)
        pN7 = np.dot(invJ, dN7)
        pN8 = np.dot(invJ, dN8)
        B1 = np.array(
            [
                [pN1[0, 0], 0, 0],
                [0, pN1[1, 0], 0],
                [0, 0, pN1[2, 0]],
                [pN1[1, 0], pN1[0, 0], 0],
                [0, pN1[2, 0], pN1[1, 0]],
                [pN1[2, 0], 0, pN1[0, 0]],
            ]
        )
        B2 = np.array(
            [
                [pN2[0, 0], 0, 0],
                [0, pN2[1, 0], 0],
                [0, 0, pN2[2, 0]],
                [pN2[1, 0], pN2[0, 0], 0],
                [0, pN2[2, 0], pN2[1, 0]],
                [pN2[2, 0], 0, pN2[0, 0]],
            ]
        )
        B3 = np.array(
            [
                [pN3[0, 0], 0, 0],
                [0, pN3[1, 0], 0],
                [0, 0, pN3[2, 0]],
                [pN3[1, 0], pN3[0, 0], 0],
                [0, pN3[2, 0], pN3[1, 0]],
                [pN3[2, 0], 0, pN3[0, 0]],
            ]
        )
        B4 = np.array(
            [
                [pN4[0, 0], 0, 0],
                [0, pN4[1, 0], 0],
                [0, 0, pN4[2, 0]],
                [pN4[1, 0], pN4[0, 0], 0],
                [0, pN4[2, 0], pN4[1, 0]],
                [pN4[2, 0], 0, pN4[0, 0]],
            ]
        )
        B5 = np.array(
            [
                [pN5[0, 0], 0, 0],
                [0, pN5[1, 0], 0],
                [0, 0, pN5[2, 0]],
                [pN5[1, 0], pN5[0, 0], 0],
                [0, pN5[2, 0], pN5[1, 0]],
                [pN5[2, 0], 0, pN5[0, 0]],
            ]
        )
        B6 = np.array(
            [
                [pN6[0, 0], 0, 0],
                [0, pN6[1, 0], 0],
                [0, 0, pN6[2, 0]],
                [pN6[1, 0], pN6[0, 0], 0],
                [0, pN6[2, 0], pN6[1, 0]],
                [pN6[2, 0], 0, pN6[0, 0]],
            ]
        )
        B7 = np.array(
            [
                [pN7[0, 0], 0, 0],
                [0, pN7[1, 0], 0],
                [0, 0, pN7[2, 0]],
                [pN7[1, 0], pN7[0, 0], 0],
                [0, pN7[2, 0], pN7[1, 0]],
                [pN7[2, 0], 0, pN7[0, 0]],
            ]
        )
        B8 = np.array(
            [
                [pN8[0, 0], 0, 0],
                [0, pN8[1, 0], 0],
                [0, 0, pN8[2, 0]],
                [pN8[1, 0], pN8[0, 0], 0],
                [0, pN8[2, 0], pN8[1, 0]],
                [pN8[2, 0], 0, pN8[0, 0]],
            ]
        )
        B = (1 / 8) * np.concatenate((B1, B2, B3, B4, B5, B6, B7, B8), axis=1)
        return B, detJ

    def stiff_linear(self, ee):
        """stiffness linear matrix"""
        
        noi = int(self.inci[ee, 4])
        noj = int(self.inci[ee, 5])
        nok = int(self.inci[ee, 6])
        nol = int(self.inci[ee, 7])
        nom = int(self.inci[ee, 8])
        non = int(self.inci[ee, 9])
        noo = int(self.inci[ee, 10])
        nop = int(self.inci[ee, 11])
        nodelist = [noi, noj, nok, nol, nom, non, noo, nop]
        D = get_elasticity(self.tabmat, self.inci, ee)
        xp, wp = self.quadra
        xpp = xp[0, :]
        ypp = xp[1, :]
        zpp = xp[2, :]
        keh8 = np.zeros((self.dofe, self.dofe))
        V = 0
        for pp in range(0, self.npp):
            intpl = [xpp[pp], ypp[pp], zpp[pp]]
            B, detJ = Solid81.matriz_b(self, nodelist, intpl)
            keh8 += (
                np.dot(np.dot(np.transpose(B), D), B) * detJ * wp[pp] * wp[pp] * wp[pp]
            )
            V += detJ
        loc = Solid81.lockey(self, nodelist)
        return keh8, loc

    def mass(self, ee):
        """consistent mass matrix"""
        
        noi = int(self.inci[ee, 4])
        noj = int(self.inci[ee, 5])
        nok = int(self.inci[ee, 6])
        nol = int(self.inci[ee, 7])
        nom = int(self.inci[ee, 8])
        non = int(self.inci[ee, 9])
        noo = int(self.inci[ee, 10])
        nop = int(self.inci[ee, 11])
        xi = self.coord[noi - 1, 1]
        yi = self.coord[noi - 1, 2]
        zi = self.coord[noi - 1, 3]
        xj = self.coord[noj - 1, 1]
        yj = self.coord[noj - 1, 2]
        zj = self.coord[noj - 1, 3]
        xk = self.coord[nok - 1, 1]
        yk = self.coord[nok - 1, 2]
        zk = self.coord[nok - 1, 3]
        xl = self.coord[nol - 1, 1]
        yl = self.coord[nol - 1, 2]
        zl = self.coord[nol - 1, 3]
        xm = self.coord[nom - 1, 1]
        ym = self.coord[nom - 1, 2]
        zm = self.coord[nom - 1, 3]
        xn = self.coord[non - 1, 1]
        yn = self.coord[non - 1, 2]
        zn = self.coord[non - 1, 3]
        xo = self.coord[noo - 1, 1]
        yo = self.coord[noo - 1, 2]
        zo = self.coord[noo - 1, 3]
        xp = self.coord[nop - 1, 1]
        yp = self.coord[nop - 1, 2]
        zp = self.coord[nop - 1, 3]
        matXY = np.array(
            [
                [xi, yi, zi],
                [xj, yj, zj],
                [xk, yk, zk],
                [xl, yl, zl],
                [xm, ym, zm],
                [xn, yn, zn],
                [xo, yo, zo],
                [xp, yp, zp],
            ]
        )
        R = self.tabmat[int(self.inci[ee, 2] - 1), 6]
        xp, wp = self.quadra
        x = xp[0, :]
        y = xp[1, :]
        z = xp[2, :]
        meh8 = np.zeros((self.dofe, self.dofe))
        for pp in range(0, self.npp):
            N1 = (1 - x[pp]) * (1 - y[pp]) * (1 - z[pp])
            N2 = (1 + x[pp]) * (1 - y[pp]) * (1 - z[pp])
            N3 = (1 + x[pp]) * (1 + y[pp]) * (1 - z[pp])
            N4 = (1 - x[pp]) * (1 + y[pp]) * (1 - z[pp])
            N5 = (1 - x[pp]) * (1 - y[pp]) * (1 + z[pp])
            N6 = (1 + x[pp]) * (1 - y[pp]) * (1 + z[pp])
            N7 = (1 + x[pp]) * (1 + y[pp]) * (1 + z[pp])
            N8 = (1 - x[pp]) * (1 + y[pp]) * (1 + z[pp])
            N = (1 / 8) * np.array(
                [
                    [
                        N1,
                        0,
                        0,
                        N2,
                        0,
                        0,
                        N3,
                        0,
                        0,
                        N4,
                        0,
                        0,
                        N5,
                        0,
                        0,
                        N6,
                        0,
                        0,
                        N7,
                        0,
                        0,
                        N8,
                        0,
                        0,
                    ],
                    [
                        0,
                        N1,
                        0,
                        0,
                        N2,
                        0,
                        0,
                        N3,
                        0,
                        0,
                        N4,
                        0,
                        0,
                        N5,
                        0,
                        0,
                        N6,
                        0,
                        0,
                        N7,
                        0,
                        0,
                        N8,
                        0,
                    ],
                    [
                        0,
                        0,
                        N1,
                        0,
                        0,
                        N2,
                        0,
                        0,
                        N3,
                        0,
                        0,
                        N4,
                        0,
                        0,
                        N5,
                        0,
                        0,
                        N6,
                        0,
                        0,
                        N7,
                        0,
                        0,
                        N8,
                    ],
                ]
            )
            N1x = -(1 - y[pp]) * (1 - z[pp])
            N2x = (1 - y[pp]) * (1 - z[pp])
            N3x = (1 + y[pp]) * (1 - z[pp])
            N4x = -(1 + y[pp]) * (1 - z[pp])
            N5x = -(1 - y[pp]) * (1 + z[pp])
            N6x = (1 - y[pp]) * (1 + z[pp])
            N7x = (1 + y[pp]) * (1 + z[pp])
            N8x = -(1 + y[pp]) * (1 + z[pp])
            N1y = -(1 - x[pp]) * (1 - z[pp])
            N2y = -(1 + x[pp]) * (1 - z[pp])
            N3y = (1 + x[pp]) * (1 - z[pp])
            N4y = (1 - x[pp]) * (1 - z[pp])
            N5y = -(1 - x[pp]) * (1 + z[pp])
            N6y = -(1 + x[pp]) * (1 + z[pp])
            N7y = (1 + x[pp]) * (1 + z[pp])
            N8y = (1 - x[pp]) * (1 + z[pp])
            N1z = -(1 - x[pp]) * (1 - y[pp])
            N2z = -(1 + x[pp]) * (1 - y[pp])
            N3z = -(1 + x[pp]) * (1 + y[pp])
            N4z = -(1 - x[pp]) * (1 + y[pp])
            N5z = (1 - x[pp]) * (1 - y[pp])
            N6z = (1 + x[pp]) * (1 - y[pp])
            N7z = (1 + x[pp]) * (1 + y[pp])
            N8z = (1 - x[pp]) * (1 + y[pp])
            dN = (1 / 8) * np.array(
                [
                    [N1x, N2x, N3x, N4x, N5x, N6x, N7x, N8x],
                    [N1y, N2y, N3y, N4y, N5y, N6y, N7y, N8y],
                    [N1z, N2z, N3z, N4z, N5z, N6z, N7z, N8z],
                ]
            )
            J = np.dot(dN, matXY)
            detJ = np.linalg.det(J)
            meh8 += np.dot(np.transpose(N), N) * R * detJ * wp[pp] * wp[pp]
        list_node = [noi, noj, nok, nol, nom, non, noo, nop]
        loc = Solid81.lockey(self, list_node)
        return meh8, loc


if __name__ == "__main__":
    import doctest

    doctest.testmod()
