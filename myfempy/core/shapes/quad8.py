from __future__ import annotations

from numpy import sqrt, array, zeros
from numpy.linalg import norm

from myfempy.core.shapes.quad8_tasks import (DiffShapeFuntion, Jacobian,
                                             LocKey, NodeCoord, NodeList,
                                             ShapeFunctions, detJacobi,
                                             invJacobi)
from myfempy.core.shapes.shape import Shape


class Quad8(Shape):
    """Quadrilateral 8-Node Shape Class <ConcreteClassService>"""

    def getShapeSet():
        shapeset = {
            "def": "8-nodes_conec 2-interpol_order",
            "key": "quad8",
            "id": 82,
            "nodes": ["i", "j", "k", "l", "m", "n", "o", "p"],
            "sidenorm": {"0": [0, -1], "1": [1, 0], "2": [0, 1], "3": [-1, 0]},
            "nodesconecedge": 3,
        }
        return shapeset

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

    def getSideAxis(set_side):
        side = {
            "0 1 4": "0",
            # "1 0 4": "0",
            "1 2 5": "1",
            # "2 1 5": "1",
            "2 3 6": "2",
            # "3 2 6": "2",
            # "3 0 7": "3",
            "0 3 7": "3",
        }
        return side[set_side]
    
    def getNormalEdge(elementcoord, side):
        
        nodes_conec_dic = {
            '0': [0, 1],
            '1': [1, 2],
            '2': [2, 3],
            '3': [0, 3],
        }
        
        nodes_conec = nodes_conec_dic[side]
        
        noi = nodes_conec[0]
        noj = nodes_conec[1]
        
        normal = zeros((2))
        dx = elementcoord[noj, 0] - elementcoord[noi, 0]
        dy = elementcoord[noj, 1] - elementcoord[noi, 1]
        L = sqrt(dx**2 + dy**2)

        normal[0] = -dy / L
        normal[1] = dx / L
        
        return normal

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