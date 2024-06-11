from __future__ import annotations

import numdifftools as nd
import numpy as np

from myfempy.core.shapes.shape import Shape


class Hexa8(Shape):
    """Hexaedro 8-Node Shape Class <ConcreteClassService>"""

    def getShapeSet():
        shapeset = {
            "def": "8-nodes_conec 1-first_order",
            "key": "hexa8",
            "id": 81,
            "nodes": ["i", "j", "k", "l", "m", "n", "o", "p"],
        }
        return shapeset

    def N(r_coord):
        N = np.zeros((1, 8))
        N[0, 0] = 0.125 * (1 - r_coord[0]) * (1 - r_coord[1]) * (1 - r_coord[2])
        N[0, 1] = 0.125 * (1 + r_coord[0]) * (1 - r_coord[1]) * (1 - r_coord[2])
        N[0, 2] = 0.125 * (1 + r_coord[0]) * (1 + r_coord[1]) * (1 - r_coord[2])
        N[0, 3] = 0.125 * (1 - r_coord[0]) * (1 + r_coord[1]) * (1 - r_coord[2])
        N[0, 4] = 0.125 * (1 - r_coord[0]) * (1 - r_coord[1]) * (1 + r_coord[2])
        N[0, 5] = 0.125 * (1 + r_coord[0]) * (1 - r_coord[1]) * (1 + r_coord[2])
        N[0, 6] = 0.125 * (1 + r_coord[0]) * (1 + r_coord[1]) * (1 + r_coord[2])
        N[0, 7] = 0.125 * (1 - r_coord[0]) * (1 + r_coord[1]) * (1 + r_coord[2])
        return N

    def diffN(shape_function, r_coord):
        dN = nd.Gradient(shape_function, n=1)
        return dN(r_coord)

    def getShapeFunctions(r_coord, nodedof):
        shape_function = Hexa8.N(r_coord)
        mat_N = np.zeros((nodedof, 8 * nodedof))
        for block in range(8):
            for dof in range(nodedof):
                mat_N[dof, block * nodedof + dof] = shape_function[0, block]
        return mat_N

    def getDiffShapeFuntion(shape_function, r_coord, nodedof):
        diff_shape_function = Hexa8.diffN(shape_function, r_coord)
        mat_diff_N = np.zeros((3 * nodedof, 8 * nodedof))
        for block in range(8):
            for dof in range(nodedof):
                mat_diff_N[
                    nodedof * dof - dof * (nodedof - 3), block * nodedof + dof
                ] = diff_shape_function[0, block]
                mat_diff_N[
                    nodedof * dof - dof * (nodedof - 3) + 1, block * nodedof + dof
                ] = diff_shape_function[1, block]
                mat_diff_N[
                    nodedof * dof - dof * (nodedof - 3) + 2, block * nodedof + dof
                ] = diff_shape_function[2, block]
        return mat_diff_N

    def getJacobian(shape_function, r_coord, element_coord):
        return np.dot(Hexa8.diffN(shape_function, r_coord), element_coord)

    def invJacobi(shape_function, r_coord, element_coord, nodedof):
        J = Hexa8.getJacobian(shape_function, r_coord, element_coord)
        invJ = np.linalg.inv(J)
        mat_invJ = np.kron(np.eye(nodedof, dtype=int), invJ)
        return mat_invJ

    def detJacobi(shape_function, r_coord, element_coord):
        J = Hexa8.getJacobian(shape_function, r_coord, element_coord)
        return np.linalg.det(J)

    def getNodeList(inci, element_number):
        noi = int(inci[element_number, 4])
        noj = int(inci[element_number, 5])
        nok = int(inci[element_number, 6])
        nol = int(inci[element_number, 7])
        nom = int(inci[element_number, 8])
        non = int(inci[element_number, 9])
        noo = int(inci[element_number, 10])
        nop = int(inci[element_number, 11])

        node_list = [noi, noj, nok, nol, nom, non, noo, nop]

        return node_list

    def getNodeCoord(coord, nodelist):
        noi = nodelist[0]
        noj = nodelist[1]
        nok = nodelist[2]
        nol = nodelist[3]
        nom = nodelist[4]
        non = nodelist[5]
        noo = nodelist[6]
        nop = nodelist[7]

        xi = coord[noi - 1, 1]
        yi = coord[noi - 1, 2]
        zi = coord[noi - 1, 3]
        xj = coord[noj - 1, 1]
        yj = coord[noj - 1, 2]
        zj = coord[noj - 1, 3]
        xk = coord[nok - 1, 1]
        yk = coord[nok - 1, 2]
        zk = coord[nok - 1, 3]
        xl = coord[nol - 1, 1]
        yl = coord[nol - 1, 2]
        zl = coord[nol - 1, 3]
        xm = coord[nom - 1, 1]
        ym = coord[nom - 1, 2]
        zm = coord[nom - 1, 3]
        xn = coord[non - 1, 1]
        yn = coord[non - 1, 2]
        zn = coord[non - 1, 3]
        xo = coord[noo - 1, 1]
        yo = coord[noo - 1, 2]
        zo = coord[noo - 1, 3]
        xp = coord[nop - 1, 1]
        yp = coord[nop - 1, 2]
        zp = coord[nop - 1, 3]

        element_coord = np.array(
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

        return element_coord

    def getShapeKey(nodelist, nodedof):
        """element lockey(dof)"""

        shape_key = np.zeros(8 * nodedof, dtype=int)

        for node in range(len(nodelist)):
            for dof in range(nodedof):
                shape_key[nodedof * node + dof] = nodedof * nodelist[node] - (
                    nodedof - dof
                )

        return shape_key
