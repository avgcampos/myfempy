from __future__ import annotations

from numpy import sqrt

from myfempy.core.shapes.line2_tasks import (DiffDiffShapeFuntion,
                                             DiffShapeFuntion, Jacobian,
                                             LocKey, NodeCoord, NodeList,
                                             ShapeFunctions, detJacobi,
                                             invJacobi)
from myfempy.core.shapes.shape import Shape

# from myfempy.core.utilities import getRotational_3dVector


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


class Line2(Shape):
    """Line 2-Node Shape Class <ConcreteClassService>"""

    def getShapeSet():
        shapeset = {
            "def": "2-nodes_conec 1-interpol_order",
            "key": "line2",
            "id": 21,
            "nodes": ["i", "j"],
            "sidenorm": {"0": [0, 1]},
        }
        return shapeset

    # quad4 sides
    def getIsoParaSide(side, r):
        # [r_valor, r_axis]
        # r = 0/ s = 1/ t = 2
        isops = {
            "0": [r, 0.0],  # [r, 0]
        }

        return isops[side]

    def getShapeFunctions(r_coord, nodedof, detJ):
        shapeN = ShapeFunctions(r_coord, nodedof)
        shapeN[1,5] = 2 * detJ * shapeN[1,5] #* wt[ip]
        shapeN[1,11] = 2 * detJ * shapeN[1,11] #* wt[ip]
        shapeN[2,4] = 2 * detJ * shapeN[2,4] #* wt[ip]
        shapeN[2,10] = 2 * detJ * shapeN[2,10] #* wt[ip]
        return shapeN

    def getDiffShapeFuntion(r_coord, nodedof):
        return DiffShapeFuntion(r_coord, nodedof)

    def getDiffDiffShapeFuntion(r_coord, nodedof, detJ):
        diffN = DiffDiffShapeFuntion(r_coord, nodedof)
        diffN[1,5] = 2 * detJ * diffN[1,5] #* wt[ip]
        diffN[1,11] = 2 * detJ * diffN[1,11] #* wt[ip]
        diffN[2,4] = 2 * detJ * diffN[2,4] #* wt[ip]
        diffN[2,10] = 2 * detJ * diffN[2,10] #* wt[ip]
        return diffN

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

        if side == "0":
            return J[0, 0]
        else:
            return 0.0

    def getSideAxis(set_side):
        side = {
            "0 1": "0",
            # "1 0": "0",
        }
        return side[set_side]
