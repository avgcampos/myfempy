from __future__ import annotations

from numpy import (
    abs,
    array,
    concatenate,
    dot,
    matmul,
    float64,
    float32,
    int32,
    ix_,
    sqrt,
    zeros,
)

INT32 = int32
FLT64 = float64
FLT32 = float32

from myfempy.core.elements.element import Element
from myfempy.core.utilities import gauss_points

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

_ELEMENT_SET = {
"def": "2D-space 2-node_dofs",
"key": "plane",
"id": 22,
"dofs": {
    "d": {"ux": 1, "uy": 2},
    "f": {
        "fx": 1,
        "fy": 2,
        "masspoint": 15,
        "spring2ground": 16,
        "damper2ground": 17,
    },
},
"tensor": ["sxx", "syy", "sxy"],
"H": array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 1, 0]], dtype=INT32),
}

class StructuralPlane(Element):
    """Plane Structural Element Class <ConcreteClassService>"""

    def getElementSet():
        return _ELEMENT_SET

    #@profile
    def getStifLinearMat(inci, coord, tabmat, tabgeo, elementcoord, C, elemdof, getIntNumK, intgauss, pt, wt, element_number):
        elem_set = StructuralPlane.getElementSet()
        H = elem_set['H']
        nodedof = len(elem_set["dofs"]["d"])
        t = tabgeo[int(inci[element_number, 3] - 1)]["THICKN"]        
        K_elem_mat = zeros((elemdof, elemdof), dtype=FLT64)
        K_elem_mat = getIntNumK(pt, wt, intgauss, elementcoord, elemdof, nodedof, H, C, t)
        return K_elem_mat


    def getMassConsistentMat(inci, coord, tabmat, tabgeo, elementcoord, C, elemdof, getIntNumM, intgauss, pt, wt, element_number):
        elem_set = StructuralPlane.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        R = tabmat[int(inci[element_number, 2]) - 1]["RHO"]
        t = tabgeo[int(inci[element_number, 3] - 1)]["THICKN"]
        M_elem_mat = zeros((elemdof, elemdof), dtype=FLT64)
        M_elem_mat = getIntNumM(pt, wt, intgauss, elementcoord, elemdof, nodedof, R, t)
        return M_elem_mat

    def getUpdateMatrix(Model, matrix, addval):
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        if int(addval[0, 1]) == 16:
            matrix_update = array([[1.0, -1.0], [-1.0, 1.0]])
        elif int(addval[0, 1]) == 15:
            matrix_update = 0.5*array([[1.0, 0.0], [0.0, 1.0]])
        for ii in range(len(addval)):
            A_add = addval[ii, 2] * matrix_update
            loc = Model.shape.getLocKey(array([addval[ii, 0]], dtype=int32), nodedof)[0:nodedof]
            matrix[ix_(loc, loc)] += A_add
        return matrix

    def getElementDeformation(U, modelinfo):
        nodetot = modelinfo["nnode"]
        nodedof = modelinfo["nodedof"]
        Udef = zeros((nodetot, 3), dtype=FLT64)
        for nn in range(1, nodetot + 1):
            Udef[nn - 1, 0] = U[nodedof * nn - 2]
            Udef[nn - 1, 1] = U[nodedof * nn - 1]
        return Udef

    def setTitleDeformation():
        return "DISPLACEMENT"

    def getElementVolume(Model, inci, coord, tabgeo, element_number):
        t = tabgeo[int(inci[element_number, 3] - 1)]["THICKN"]
        shape_set = Model.shape.getShapeSet()
        type_shape = shape_set["key"]
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        pt, wt = gauss_points(type_shape, 1)
        detJ = 0.0
        for ip in range(1):
            for jp in range(1):
                detJ += (
                    abs(Model.shape.getdetJacobi(array([pt[ip], pt[jp]]), elementcoord))
                    * wt[ip]
                    * wt[jp]
                )
        return detJ * t
