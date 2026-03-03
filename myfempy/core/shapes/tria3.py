from __future__ import annotations

from numpy import sqrt, array, zeros
from numpy.linalg import norm

from myfempy.core.shapes.shape import Shape
from myfempy.core.shapes.tria3_tasks import (DiffShapeFuntion, Jacobian,
                                             LocKey, NodeCoord, NodeList,
                                             ShapeFunctions, detJacobi,
                                             invJacobi)



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


class Tria3(Shape):
    """Triangular 3-Node Shape Class <ConcreteClassService>"""

    def getShapeSet():
        shapeset = {
            "def": "3-nodes_conec 1-interpol_order",
            "key": "tria3",
            "id": 31,
            "nodes": ["i", "j", "k"],
            "sidenorm": {"0": [0, -1], "1": [1, 1], "2": [-1, 0]},
            "nodesconecedge": 2,
        }
        return shapeset

    # tria3 sides
    def getIsoParaSide(side, r):
        # [r_valor, r_axis]
        # r = 0/ s = 1/ t = 2
        isops = {
            "0": [r, 0.0],  # [0, s]
            "1": [r, 1 - r],  # [1-r, s]
            "2": [0.0, r],  # [+1, r]
        }

        return isops[side]

    def getSideAxis(set_side):
        side = {
            "0 2": "2",
            "0 1": "0",
            "1 2": "1",
        }
        return side[set_side]
    
    def getEdgeLength(J, side):
        #   J = [dx/dr       dy/dr]
        #       [dx/ds       dy/ds]
        #   detJ_r = sqrt[dx/ds^2 + dy/ds^2]
        #   detJ_s = sqrt[dx/dr^2 + dy/dr^2]
        if side == "0":
            return sqrt(J[0, 0] ** 2 + J[0, 1] ** 2)
        elif side == "1":
            return sqrt((J[0, 0] - J[1, 0]) ** 2 + (J[0, 1] - J[1, 1]) ** 2)
        elif side == "2":
            return sqrt(J[1, 0] ** 2 + J[1, 1] ** 2)
        else:
            return 0.0
    
    
    def getNormalEdge(elementcoord, side):
        
        nodes_conec_dic = {
            '0': [0, 1],
            '1': [1, 2],
            '2': [2, 0],
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




