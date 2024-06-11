from __future__ import annotations

import numdifftools as nd
import numpy as np

from myfempy.core.shapes.shape import Shape


class Line2(Shape):
    """Linear 2-Node Shape Class <ConcreteClassService>"""

    def getShapeSet():
        shapeset = {
            "def": "2-nodes_conec 1-first_order",
            "key": "line2",
            "id": 21,
            "nodes": ["i", "j"],
        }
        return shapeset

    def N(r_coord):
        N = np.zeros((1, 2))
        N[0, 0] = 0.5 * (1.0 - r_coord)
        N[0, 1] = 0.5 * (1.0 + r_coord)
        return N

    def diffN(shape_function, r_coord):
        dN = nd.Gradient(shape_function, n=1)
        return dN(r_coord)

    def getShapeFunctions(r_coord, nodedof):
        shape_function = Line2.N(r_coord)

        mat_N = np.zeros((nodedof, 1 * nodedof))
        for block in range(1):
            for dof in range(nodedof):
                mat_N[dof, block * nodedof + dof] = shape_function[0, block]

        return mat_N

    def getDiffShapeFuntion(shape_function, r_coord, nodedof):
        diff_shape_function = Line2.diffN(shape_function, r_coord)

        mat_diff_N = np.zeros((1 * nodedof, 1 * nodedof))

        for block in range(1):
            for dof in range(nodedof):
                mat_diff_N[
                    nodedof * dof - dof * (nodedof - 1), block * nodedof + dof
                ] = diff_shape_function[0, block]

        return mat_diff_N

    def getJacobian(shape_function, r_coord, element_coord):
        return np.dot(Line2.diffN(shape_function, r_coord), element_coord)

    def invJacobi(shape_function, r_coord, element_coord, nodedof):
        J = Line2.getJacobian(shape_function, r_coord, element_coord)

        invJ = np.linalg.inv(J)

        mat_invJ = np.kron(np.eye(nodedof, dtype=int), invJ)

        return mat_invJ

    def detJacobi(shape_function, r_coord, element_coord):
        J = Line2.getJacobian(shape_function, r_coord, element_coord)
        return np.linalg.det(J)

    def getNodeList(inci, element_number):
        noi = int(inci[element_number, 4])
        noj = int(inci[element_number, 5])

        node_list = [noi, noj]

        return node_list

    def getNodeCoord(coord, node_list):
        noi = node_list[0]
        noj = node_list[1]

        xi = coord[noi - 1, 1]
        yi = coord[noi - 1, 2]
        xj = coord[noj - 1, 1]
        yj = coord[noj - 1, 2]

        element_coord = np.array([[xi, yi], [xj, yj]])
        return element_coord

    def getShapeKey(node_list, nodedof):
        """element lockey(dof)"""
        shape_key = np.zeros(2 * nodedof, dtype=int)

        for node in range(len(node_list)):
            for dof in range(nodedof):
                shape_key[nodedof * node + dof] = nodedof * node_list[node] - (
                    nodedof - dof
                )
        return shape_key
