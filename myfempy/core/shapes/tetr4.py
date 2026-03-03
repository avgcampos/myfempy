from __future__ import annotations

from numpy import sqrt, array, cross
from numpy.linalg import norm

from myfempy.core.shapes.shape import Shape
from myfempy.core.shapes.tetr4_tasks import (DiffShapeFuntion, Jacobian,
                                             LocKey, NodeCoord, NodeList,
                                             ShapeFunctions, detJacobi,
                                             invJacobi)
from myfempy.core.utilities import poly_area, unit_normal


__docformat__ = "google"

__doc__ = """

==========================================================================
                            __                                
         _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
        | '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
        | | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
        |_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
                    |___/                       |_|     |___/ 
        myfempy -- MultiphYsics Finite Element Module to PYthon    
                    COMPUTATIONAL ANALYSIS PROGRAM                   
        Copyright (C) 2022-2026 Antonio Vinicius Garcia Campos        
==========================================================================
This Python file is part of myfempy project.

myfempy is a python package based on finite element method to multiphysics
analysis. The code is open source and *intended for educational and scientific
purposes only, not recommended to commercial use. The name myfempy is an acronym
for MultiphYsics Finite Elements Module to PYthon. You can help us by contributing
with the main project, send us a mensage on https://github.com/avgcampos/myfempy/discussions/10
If you use myfempy in your research, the  developers would be grateful if you 
could cite in your work.
																		
The code is written by Antonio Vinicius Garcia Campos.                                  
																		
A github repository, with the most up to date version of the code,      
can be found here: https://github.com/avgcampos/myfempy.                 
																		
The code is open source and intended for educational and scientific     
purposes only. If you use myfempy in your research, the developers      
would be grateful if you could cite this. The myfempy project is published
under the GPLv3, see the myfempy LICENSE on
https://github.com/avgcampos/myfempy/blob/main/LICENSE.
																		
Disclaimer:                                                             
The authors reserve all rights but do not guarantee that the code is    
free from errors. Furthermore, the authors shall not be liable in any   
event caused by the use of the program.

"""


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
    
    def getIsoParaSide(side, r):
        isops = {
            "0": [r[0], r[1], 0.0],  
            "1": [r[0], 0.0, r[1]],  
            "2": [0.0, r[0], r[1]],  
            "3": [r[0], 1 - r[0] - r[1], r[1]], 
        }
        return isops[side]
    
    def getAreaLength(side, elementcoord):
        
        nodes_conec_dic = {
            '0': [0, 1, 2],
            '1': [0, 3, 1],
            '2': [0, 2, 3],
            '3': [3, 2, 1],
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
            '0 1 3': '1',
            '0 2 3': '2',
            '1 2 3': '3',
        }
        return side[set_side]

    def getNormalFace(elementcoord, side):
        
        nodes_conec_dic = {
            '0': [0, 1, 2],
            '1': [0, 3, 1],
            '2': [0, 2, 3],
            '3': [3, 2, 1],
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

        normal = unit_normal(poly[0], poly[1], poly[2])
        normal = array(normal)
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
