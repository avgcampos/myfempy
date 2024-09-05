from __future__ import annotations

from numpy import sqrt

from myfempy.core.shapes.quad4_tasks import (DiffShapeFuntion, Jacobian,
                                             LocKey, NodeCoord, NodeList,
                                             ShapeFunctions, detJacobi,
                                             invJacobi)
from myfempy.core.shapes.shape import Shape


class Quad4(Shape):
    """Quadrilateral 4-Node Shape Class <ConcreteClassService>"""

    def getShapeSet():
        shapeset = {
            "def": "4-nodes_conec 1-first_order",
            "key": "quad4",
            "id": 41,
            "nodes": ["i", "j", "k", "l"],
            "sidenorm": {"0": [0, -1], "1": [1, 0], "2": [0, 1], "3": [-1, 0]},
        }
        return shapeset

    # quad4 sides
    def getIsoParaSide(side, r):
        # [r_valor, r_axis]
        # r = 0/ s = 1/ t = 2
        isops = {
            "0": [r, -1.0],  # [-1, s]
            "1": [1.0, r],  # [+1, r]
            "2": [r, 1.0],  # [+1, s]
            "3": [-1.0, r],  # [-1, r]
        }

        return isops[side]

    def getShapeFunctions(r_coord, nodedof):
        return ShapeFunctions(r_coord, nodedof)

    def getDiffShapeFuntion(r_coord, nodedof):
        return DiffShapeFuntion(r_coord, nodedof)

    def getJacobian(r_coord, element_coord):
        return Jacobian(r_coord, element_coord)

    def getinvJacobi(r_coord, element_coord, nodedof):
        return invJacobi(r_coord, element_coord, nodedof)

    def getdetJacobi(r_coord, element_coord):
        return detJacobi(r_coord, element_coord)

    def getNodeList(inci, element_number):
        return NodeList(inci, element_number)

    def getNodeCoord(coord, node_list):
        return NodeCoord(coord, node_list)

    def getLocKey(node_list, nodedof):
        return LocKey(node_list, nodedof)

    def getEdgeLength(J, side):
        #   J = [dx/dr       dy/dr]
        #       [dx/ds       dy/ds]
        #   detJ_r = sqrt[dx/ds^2 + dy/ds^2]
        #   detJ_s = sqrt[dx/dr^2 + dy/dr^2]
        if side == "0" or side == "2":
            return sqrt(J[0, 0] ** 2 + J[0, 1] ** 2)
        elif side == "1" or side == "3":
            return sqrt(J[1, 0] ** 2 + J[1, 1] ** 2)
        else:
            return 0.0
