from __future__ import annotations

import numdifftools as nd
import numpy as np

from myfempy.core.shapes.shape import Shape


class Tetra4(Shape):
    """Tetraedro 4-Node Shape Class <ConcreteClassService>"""

    def getShapeSet():
        shapeset = {
            "def": "4-nodes_conec 1-first_order",
            "key": "tetr4",
            "id": 41,
            "nodes": ["i", "j", "k", "l"],
        }
        return shapeset

    def N(r_coord):
        N = np.zeros((1, 4))
        N[0, 0] = 1 - r_coord[0] - r_coord[1] - r_coord[2]
        N[0, 1] = r_coord[0]
        N[0, 2] = r_coord[1]
        N[0, 3] = r_coord[2]
        return N

    def diffN(shape_function, r_coord):
        dN = nd.Gradient(shape_function, n=1)
        return dN(r_coord)

    def getShapeFunctions(r_coord, nodedof):
        shape_function = Tetra4.N(r_coord)
        mat_N = np.zeros((nodedof, 4 * nodedof))
        for block in range(4):
            for dof in range(nodedof):
                mat_N[dof, block * nodedof + dof] = shape_function[0, block]
        return mat_N

    def getDiffShapeFuntion(shape_function, r_coord, nodedof):
        diff_shape_function = Tetra4.diffN(shape_function, r_coord)
        mat_diff_N = np.zeros((3 * nodedof, 4 * nodedof))
        for block in range(4):
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
        return np.dot(Tetra4.diffN(shape_function, r_coord), element_coord)

    def invJacobi(shape_function, r_coord, element_coord, nodedof):
        J = Tetra4.getJacobian(shape_function, r_coord, element_coord)
        invJ = np.linalg.inv(J)
        mat_invJ = np.kron(np.eye(nodedof, dtype=int), invJ)
        return mat_invJ

    def detJacobi(shape_function, r_coord, element_coord):
        J = Tetra4.getJacobian(shape_function, r_coord, element_coord)
        return (1 / 6) * np.linalg.det(J)

    def getNodeList(inci, element_number):
        noi = int(inci[element_number, 4])
        noj = int(inci[element_number, 5])
        nok = int(inci[element_number, 6])
        nol = int(inci[element_number, 7])

        node_list = [noi, noj, nok, nol]

        return node_list

    def getNodeCoord(coord, nodelist):
        noi = nodelist[0]
        noj = nodelist[1]
        nok = nodelist[2]
        nol = nodelist[3]

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

        element_coord = np.array(
            [
                [xi, yi, zi],
                [xj, yj, zj],
                [xk, yk, zk],
                [xl, yl, zl],
            ]
        )

        return element_coord

    def getShapeKey(nodelist, nodedof):
        """element lockey(dof)"""

        shape_key = np.zeros(4 * nodedof, dtype=int)

        for node in range(len(nodelist)):
            for dof in range(nodedof):
                shape_key[nodedof * node + dof] = nodedof * nodelist[node] - (
                    nodedof - dof
                )

        return shape_key
