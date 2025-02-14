from __future__ import annotations

from numpy import sqrt, array, cross
from numpy.linalg import norm

from myfempy.core.shapes.shape import Shape
from myfempy.core.shapes.tetr4_tasks import (DiffShapeFuntion, Jacobian,
                                             LocKey, NodeCoord, NodeList,
                                             ShapeFunctions, detJacobi,
                                             invJacobi)
from myfempy.core.utilities import poly_area

class Tetra4(Shape):
    """Tetrahedron 4-Node Shape Class <ConcreteClassService>"""

    def getShapeSet():
        shapeset = {
            "def": "4-nodes_conec 1-interpol_order",
            "key": "tetr4",
            "id": 41,
            "nodes": ["i", "j", "k", "l"],
            "sidenorm": {
                "0": [0, 0, -1],
                "1": [0, -1, 0],
                "2": [-1, 0, 0],
                "3": [1, 1, 1],
            },
            "nodesconecface": 3,
        }
        return shapeset

    # tetra4 sides
    def getIsoParaSide(side, r):
        # r = 0/ s = 1/ t = 2
        isops = {
            "0": [r[0], r[1], 0.0],  # [r, s, 0]
            "1": [r[0], 0.0, r[1]],  # [r, 0, t]
            "2": [0.0, r[0], r[1]],  # [0, s, t]
            "3": [r[0], 1 - r[0] - r[1], r[1]],  # [r, 1 - r, 1 - r]
        }
        return isops[side]
    
    def getAreaLength(side, elementcoord):
        
        nodes_conec_dic = {
            '0': [0, 1, 2],
            '1': [0, 1, 3],
            '2': [0, 2, 3],
            '3': [1, 2, 3],
        }
        
        nodes_conec = nodes_conec_dic[side]
                        
        no1 = nodes_conec[0]
        no2 = nodes_conec[1]
        no3 = nodes_conec[2]
        
        coord_x_no1 = elementcoord[no1, 0]
        coord_y_no1 = elementcoord[no1, 1]
        coord_z_no1 = elementcoord[no1, 2]
        coord_x_no2 = elementcoord[no2, 0]
        coord_y_no2 = elementcoord[no2, 1]
        coord_z_no2 = elementcoord[no2, 2]
        coord_x_no3 = elementcoord[no3, 0]
        coord_y_no3 = elementcoord[no3, 1]
        coord_z_no3 = elementcoord[no3, 2]
        
        poly = array(
            [
                [coord_x_no1, coord_y_no1, coord_z_no1],
                [coord_x_no2, coord_y_no2, coord_z_no2],
                [coord_x_no3, coord_y_no3, coord_z_no3],
            ]
        )
        area_surf = poly_area(poly)
        return area_surf

    def getSideAxis(set_side):
        side = {
            '0 1 2': '0',
            # '0 2 1': '0',
            # '1 0 2': '0',
            # '1 2 0': '0',
            # '2 0 1': '0',
            # '2 1 0': '0',
            '0 1 3': '1',
            # '1 0 3': '1',
            # '3 1 0': '1',
            # '3 0 1': '1',
            # '0 1 3': '1',
            # '0 3 1': '1',
            '0 2 3': '2',
            # '0 3 2': '2',
            # '2 0 3': '2',
            # '2 3 0': '2',
            # '3 0 2': '2',
            # '3 2 0': '2',
            '1 2 3': '3',
            # '1 3 2': '3',
            # '2 1 3': '3',
            # '2 3 1': '3',
            # '3 1 2': '3',
            # '3 2 1': '3',
        }
        return side[set_side]

    def getNormalFace(elementcoord, get_side):
        faces = array([
                [elementcoord[0], elementcoord[1], elementcoord[2]],  # Face 0
                [elementcoord[0], elementcoord[1], elementcoord[3]],  # Face 1
                [elementcoord[0], elementcoord[2], elementcoord[3]],  # Face 2
                [elementcoord[1], elementcoord[2], elementcoord[3]]   # Face 3
            ])
        face = faces[int(get_side)]
        # Compute the outward normal vector of the face
        v1 = face[1] - face[0]
        v2 = face[2] - face[0]
        normal = cross(v1, v2)
        normal *= 1.0/norm(normal)  # Normalize
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
