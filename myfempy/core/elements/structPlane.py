from __future__ import annotations

from numpy import (abs, array, concatenate, dot, float64, int32, ix_, sqrt,
                   zeros)

# from os import environ
# environ["OMP_NUM_THREADS"] = "8"


INT32 = int32
FLT64 = float64

from myfempy.core.elements.element import Element
from myfempy.core.utilities import gauss_points


class StructuralPlane(Element):
    """Plane Structural Element Class <ConcreteClassService>"""

    def getElementSet():
        elemset = {
            "def": "2D-space 2-node_dofs",
            "key": "plane",
            "id": 22,
            "dofs": {
                "d": {"ux": 1, "uy": 2},
                "f": {
                    "fx": 1,
                    "fy": 2,
                    "massaadd": 15,
                    "spring2ground": 16,
                    "damper2ground": 17,
                },
            },
            "tensor": ["sxx", "syy", "sxy"],
        }
        return elemset

    # def getH():
    #     return array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 1, 0]], dtype=INT32)

    # @profile
    def getB(Model, elementcoord, pt, nodedof):
        diffN = Model.shape.getDiffShapeFuntion(pt, nodedof)
        invJ = Model.shape.getinvJacobi(pt, elementcoord, nodedof)
        H = array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 1, 0]], dtype=INT32)
        B = dot(H, dot(invJ, diffN))
        return B

    # @profile
    def getStifLinearMat(Model, inci, coord, tabmat, tabgeo, intgauss, element_number):
        elem_set = StructuralPlane.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        type_shape = shape_set["key"]
        edof = nodecon * nodedof
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        E = tabmat[int(inci[element_number, 2]) - 1][
            "EXX"
        ]  # tabmat[int(inci[element_number, 2]) - 1, 0]  # material elasticity
        v = tabmat[int(inci[element_number, 2]) - 1][
            "VXX"
        ]  # tabmat[int(inci[element_number, 2]) - 1, 1]  # material poisson ratio
        C = Model.material.getElasticTensor(E, v)
        t = tabgeo[int(inci[element_number, 3] - 1)][
            "THICKN"
        ]  # tabgeo[int(inci[element_number, 3] - 1), 4]
        pt, wt = gauss_points(type_shape, intgauss)
        K_elem_mat = zeros((edof, edof), dtype=FLT64)
        for ip in range(intgauss):
            for jp in range(intgauss):
                detJ = Model.shape.getdetJacobi(array([pt[ip], pt[jp]]), elementcoord)
                B = StructuralPlane.getB(
                    Model, elementcoord, array([pt[ip], pt[jp]]), nodedof
                )
                BCB = dot(dot(B.transpose(), C), B)
                K_elem_mat += BCB * t * abs(detJ) * wt[ip] * wt[jp]
        return K_elem_mat

    def getMassConsistentMat(
        Model, inci, coord, tabmat, tabgeo, intgauss, element_number
    ):
        elem_set = StructuralPlane.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        type_shape = shape_set["key"]
        edof = nodecon * nodedof
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        R = tabmat[int(inci[element_number, 2]) - 1][
            "RHO"
        ]  # tabmat[int(inci[element_number, 2]) - 1, 6]  # material density
        t = tabgeo[int(inci[element_number, 3] - 1)][
            "THICKN"
        ]  # tabgeo[int(inci[element_number, 3] - 1), 4]
        pt, wt = gauss_points(type_shape, intgauss)
        M_elem_mat = zeros((edof, edof), dtype=FLT64)
        for ip in range(intgauss):
            for jp in range(intgauss):
                detJ = Model.shape.getdetJacobi(array([pt[ip], pt[jp]]), elementcoord)
                N = Model.shape.getShapeFunctions(array([pt[ip], pt[jp]]), nodedof)
                NRN = dot(dot(N.transpose(), R), N)
                M_elem_mat += NRN * t * abs(detJ) * wt[ip] * wt[jp]
        return M_elem_mat

    def getUpdateMatrix(Model, matrix, addval):
        elem_set = Model.element.getElementSet()
        shape_set = Model.shape.getShapeSet()
        dofe = len(shape_set["nodes"]) * len(elem_set["dofs"]["d"])
        for ii in range(len(addval)):

            A_add = addval[ii, 2] * array([[1, -1], [-1, 1]])

            loc = array(
                [
                    int(dofe * addval[ii, 0] - (dofe)),
                    int(dofe * addval[ii, 0] - (dofe - 1)),
                ]
            )

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
