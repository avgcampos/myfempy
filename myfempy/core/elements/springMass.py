from __future__ import annotations

from numpy import (ndarray, array, ascontiguousarray, dot, eye, float64, int32, ix_,
                   sqrt, transpose, zeros)
from scipy.linalg import block_diag

INT32 = int32
FLT64 = float64

from myfempy.core.elements.element import Element
from myfempy.core.utilities import (gauss_points, get3D_LocalVector,
                                    getRotational_Matrix)

__docformat__ = "google"

__doc__ = """
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
    "def": "3D-space 4-node_dofs",
    "key": "springmass",
    "id": 34,
    "dofs": {
        "d": {"ux": 1, "uy": 2, "uz": 3, "rx": 4,
               },
        "f": {
            "fx": 1,
            "fy": 2,
            "fz": 3,
            "tx": 4,
            "masspoint": 15,
            "spring2ground": 16,
            "damper2ground": 17,
            "torsion2ground": 18,
        },
    },
    "tensor": ["null"],
    "H": array([[1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]], dtype=INT32),
}

def _element_coord(elementcoord):
    """(6,1) float64 C-contiguous [x1 y1 z1 x2 y2 z2] e comprimento L."""
    ec = ascontiguousarray(elementcoord, dtype=FLT64).reshape(6, 1)
    d = ec[3:6, 0] - ec[0:3, 0]
    L = float(sqrt(dot(d, d)))
    if L == 0.0:
        raise ValueError("SpringMass: elemento com comprimento nulo (nós coincidentes)")
    return L

class SpringMass(Element):
    """
    Spring Mass Linear Element Class <ConcreteClassService>
    """

    def getElementSet():
        return _ELEMENT_SET

    def getStifLinearMat(
            inci,
            coord,
            tabmat,
            tabgeo,
            elementcoord,
            D,
            elemdof,
            getIntNumK,
            intgauss,
            pt,
            wt,
            element_number,
        ):
        elem_set = SpringMass.getElementSet()
        H = elem_set["H"]
        nodedof = len(elem_set["dofs"]["d"])
        elementcoord_local = get3D_LocalVector(array(elementcoord), 2)
        L = _element_coord(elementcoord)
        STIFK = tabmat[int(inci[element_number, 2] - 1)]["STIFK"]
        STIFT = tabmat[int(inci[element_number, 2] - 1)]["STIFT"]
        STIFR = tabmat[int(inci[element_number, 2] - 1)]["STIFR"]
        C = array([[STIFK, STIFT / L, STIFT / L, STIFR]]) * eye(4)
        K_elem_mat = zeros((elemdof, elemdof), dtype=FLT64)
        K_elem_mat = getIntNumK(
            pt, wt, intgauss, elementcoord_local, elemdof, nodedof, H, C
        )
        lamb = getRotational_Matrix(elementcoord, 1)
        lambNode = block_diag(lamb, 1.0)
        R = block_diag(lambNode, lambNode)
        K_elem_mat = R.transpose().dot(K_elem_mat).dot(R)
        return K_elem_mat

    def getMassConsistentMat(
        inci,
        coord,
        tabmat,
        tabgeo,
        elementcoord,
        rho,
        elemdof,
        getIntNumM,
        intgauss,
        pt,
        wt,
        element_number,
    ):
        elem_set = SpringMass.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        elementcoord_local = get3D_LocalVector(elementcoord, 2)
        L = _element_coord(elementcoord)
        MASS = tabmat[int(inci[element_number, 2] - 1)]["MASS"]
        RHO = tabmat[int(inci[element_number, 2] - 1)]["RHO"]
        IXX = tabgeo[int(inci[element_number, 3] - 1)]["INERXX"]
        R = array([[MASS, RHO * L, RHO * L, IXX]]) * eye(4)
        M_elem_mat = zeros((elemdof, elemdof), dtype=FLT64)
        M_elem_mat = getIntNumM(
            pt, wt, intgauss, elementcoord_local, elemdof, nodedof, R
        )
        return M_elem_mat

    def getUpdateMatrix(Model, matrix, addval):
        nodedof = Model.modelinfo["nodedof"]

        if int(addval[0, 1]) == 16:
            matrix_update = array(
                [
                    [1.0, 0.0, -1.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0],
                    [-1.0, 0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0],
                ]
            )

        elif int(addval[0, 1]) == 17:
            matrix_update = array(
                [
                    [1.0, 0.0, -1.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0],
                    [-1.0, 0.0, 1.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0],
                ]
            )

        elif int(addval[0, 1]) == 18:
            matrix_update = array(
                [
                    [0.0, 0.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0, -1.0],
                    [0.0, 0.0, 0.0, 0.0],
                    [0.0, -1.0, 0.0, 1.0],
                ]
            )

        elif int(addval[0, 1]) == 15:
            matrix_update = 0.5 * array(
                [
                    [1.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0],
                    [0.0, 0.0, 0.0, 0.0],
                ]
            )

        for ii in range(len(addval)):

            A_add = addval[ii, 2] * matrix_update

            loc = Model.shape.getLocKey(array([addval[ii, 0]], dtype=int32), nodedof)[
                0:nodedof
            ]

            matrix[ix_(loc, loc)] += A_add
        return matrix


    def getElementDeformation(U, modelinfo):
        nodetot = modelinfo["nnode"]
        nodedof = modelinfo["nodedof"]
        Udef = zeros((nodetot, 3), dtype=FLT64)
        for nn in range(1, nodetot + 1):
            Udef[nn - 1, 0] = U[nodedof * nn - 6]
            Udef[nn - 1, 1] = U[nodedof * nn - 5]
            Udef[nn - 1, 2] = U[nodedof * nn - 4]
        return Udef

    def setTitleDeformation():
        return "DISPLACEMENT"

    def getElementVolume(
        inci, tabgeo, getVOL, type_shape, element_coord, element_number
    ):
        return 0.0