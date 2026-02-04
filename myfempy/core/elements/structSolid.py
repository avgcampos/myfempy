from __future__ import annotations

from os import environ

environ["OMP_NUM_THREADS"] = "1"

from numpy import (abs, array, concatenate, dot, float64, int32, ix_, sqrt,
                   zeros)

INT32 = int32
FLT64 = float64

from myfempy.core.elements.element import Element
from myfempy.core.utilities import gauss_points

class StructuralSolid(Element):
    """Solid Structural Element Class <ConcreteClassService>"""
    
    def getElementSet():
        elemset = {
            "def": "3D-space 3-node_dofs",
            "key": "solid",
            "id": 33,
            "dofs": {
                "d": {"ux": 1, "uy": 2, "uz": 3},
                "f": {
                    "fx": 1,
                    "fy": 2,
                    "fz": 3,
                    "massaadd": 15,
                    "spring2ground": 16,
                    "damper2ground": 17,
                },
            },
            "tensor": ["sxx", "syy", "szz", "sxy", "syz", "szx"],
        }
        return elemset


    def getB(diffN, invJ):
        H = array(
            [
                [1, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 1],
                [0, 1, 0, 1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1, 0, 1, 0],
                [0, 0, 1, 0, 0, 0, 1, 0, 0],
            ],
            dtype=INT32,
        )
        # B = dot(H, dot(invJ, diffN))
        B = H.dot(invJ).dot(diffN)
        # B = fastDOT(H, fastDOT(invJ, diffN))
        return B

    # @profile
    def getStifLinearMat(Model, inci, coord, tabmat, tabgeo, intgauss, element_number):
        elem_set = StructuralSolid.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        type_shape = shape_set["key"]
        edof = nodecon * nodedof
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        C = Model.material.getElasticTensor(tabmat, inci, element_number)
        pt, wt = gauss_points(type_shape, intgauss)
        K_elem_mat = zeros((edof, edof), dtype=FLT64)
        for ip in range(intgauss):
            for jp in range(intgauss):
                for kp in range(intgauss):
                    detJ = Model.shape.getdetJacobi(array([pt[ip], pt[jp], pt[kp]]), elementcoord)
                    diffN = Model.shape.getDiffShapeFuntion(array([pt[ip], pt[jp], pt[kp]]), nodedof)
                    invJ = Model.shape.getinvJacobi(array([pt[ip], pt[jp], pt[kp]]), elementcoord, nodedof)
                    B = StructuralSolid.getB(diffN, invJ)
                    BCB = B.transpose().dot(C).dot(B)
                    K_elem_mat += BCB * abs(detJ) * wt[ip] * wt[jp] * wt[kp]
        return K_elem_mat

    def getMassConsistentMat(
        Model, inci, coord, tabmat, tabgeo, intgauss, element_number
    ):
        elem_set = StructuralSolid.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        type_shape = shape_set["key"]
        edof = nodecon * nodedof
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        R = tabmat[int(inci[element_number, 2]) - 1]["RHO"]
        pt, wt = gauss_points(type_shape, intgauss)
        M_elem_mat = zeros((edof, edof), dtype=FLT64)
        for ip in range(intgauss):
            for jp in range(intgauss):
                for kp in range(intgauss):
                    detJ = Model.shape.getdetJacobi(array([pt[ip], pt[jp], pt[kp]]), elementcoord)
                    N = Model.shape.getShapeFunctions(array([pt[ip], pt[jp], pt[kp]]), nodedof)
                    NRN = N.transpose().dot(R).dot(N)
                    M_elem_mat += NRN * abs(detJ) * wt[ip] * wt[jp] * wt[kp]
        return M_elem_mat

    # def getUpdateMatrix(Model, matrix, addval):
    #     elem_set = Model.element.getElementSet()
    #     shape_set = Model.shape.getShapeSet()
    #     dofe = len(shape_set["nodes"]) * len(elem_set["dofs"]["d"])
    #     for ii in range(len(addval)):

    #         A_add = addval[ii, 2] * array([[1, -1],
    #                                        [-1, 1]])

    #         loc = array([int(dofe * addval[ii, 0] - (dofe)),
    #                     int(dofe * addval[ii, 0]  - (dofe - 1)),])

    #         matrix[ix_(loc, loc)] += A_add
    #     return matrix

    def getElementDeformation(U, modelinfo):
        nodetot = modelinfo["nnode"]
        nodedof = modelinfo["nodedof"]
        Udef = zeros((nodetot, 3), dtype=FLT64)
        for nn in range(1, nodetot + 1):
            Udef[nn - 1, 0] = U[nodedof * nn - 3]
            Udef[nn - 1, 1] = U[nodedof * nn - 2]
            Udef[nn - 1, 2] = U[nodedof * nn - 1]
        return Udef

    def setTitleDeformation():
        return "DISPLACEMENT"

    def getElementVolume(Model, inci, coord, tabgeo, element_number):
        shape_set = Model.shape.getShapeSet()
        type_shape = shape_set["key"]
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        pt, wt = gauss_points(type_shape, 1)
        Vol = 0.0
        for ip in range(1):
            for jp in range(1):
                for kp in range(1):
                    Vol += (
                        abs(
                            Model.shape.getdetJacobi(
                                array([pt[ip], pt[jp], pt[kp]]), elementcoord
                            )
                        )
                        * wt[ip]
                        * wt[jp]
                        * wt[kp]
                    )
        return Vol