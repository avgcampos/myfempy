from __future__ import annotations

from os import environ

environ["OMP_NUM_THREADS"] = "1"

from numpy import abs, array, concatenate, dot, ix_, float64, int32, sqrt, zeros

INT32 = int32
FLT64 = float64

from myfempy.core.elements.element import Element
from myfempy.core.utilities import gauss_points


def HDIFFNINVJ(H, diffN, invJ):
    invJdiffN = dot(invJ, diffN)
    B = dot(H, invJdiffN)
    return B


def BTCB(diffN, H, invJ, C):
    B = HDIFFNINVJ(H, diffN, invJ)
    BT = B.transpose()
    BTC = dot(BT, C)
    BCB = dot(BTC, B)
    return BCB


def NTRN(N, R):
    NT = N.transpose()
    NTR = dot(NT, R)
    NRN = dot(NTR, N)
    return NRN


class StructuralPlane(Element):
    """Plane Structural Element Class <ConcreteClassService>"""

    def getElementSet():
        elemset = {
            "def": "2D-space 2-node_dofs",
            "key": "plane",
            "id": 22,
            "dofs": {
                "d": {"ux": 1, "uy": 2},
                "f": {"fx": 1, "fy": 2,
                      "massaadd": 15, "spring2ground": 16, "damper2ground": 17},
            },
            "tensor": ["sxx", "syy", "sxy"],
        }
        return elemset

    def getH():
        return array([[1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 1, 0]], dtype=INT32)

    def getB(Model, elementcoord, ptg, nodedof):
        diffN = Model.shape.getDiffShapeFuntion(ptg, nodedof)
        invJ = Model.shape.getinvJacobi(ptg, elementcoord, nodedof)
        H = StructuralPlane.getH()
        B = HDIFFNINVJ(H, diffN, invJ)
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
        E = tabmat[int(inci[element_number, 2]) - 1]["EXX"] #tabmat[int(inci[element_number, 2]) - 1, 0]  # material elasticity
        v = tabmat[int(inci[element_number, 2]) - 1]["VXX"] #tabmat[int(inci[element_number, 2]) - 1, 1]  # material poisson ratio
        C = Model.material.getElasticTensor(E, v)
        t = tabgeo[int(inci[element_number, 3] - 1)]["THICKN"] #tabgeo[int(inci[element_number, 3] - 1), 4]
        H = StructuralPlane.getH()
        pt, wt = gauss_points(type_shape, intgauss)
        K_elem_mat = zeros((edof, edof), dtype=FLT64)
        for pp in range(intgauss):
            detJ = Model.shape.getdetJacobi(pt[pp], elementcoord)
            diffN = Model.shape.getDiffShapeFuntion(pt[pp], nodedof)
            invJ = Model.shape.getinvJacobi(pt[pp], elementcoord, nodedof)
            BCB = BTCB(diffN, H, invJ, C)
            K_elem_mat += BCB * t * abs(detJ) * wt[pp]
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
        R = tabmat[int(inci[element_number, 2]) - 1]["RHO"] #tabmat[int(inci[element_number, 2]) - 1, 6]  # material density
        t = tabgeo[int(inci[element_number, 3] - 1)]["THICKN"] #tabgeo[int(inci[element_number, 3] - 1), 4]
        pt, wt = gauss_points(type_shape, intgauss)
        M_elem_mat = zeros((edof, edof), dtype=FLT64)
        for pp in range(intgauss):
            detJ = Model.shape.getdetJacobi(pt[pp], elementcoord)
            N = Model.shape.getShapeFunctions(pt[pp], nodedof)
            NRN = NTRN(N, R)
            M_elem_mat += NRN * t * abs(detJ) * wt[pp]
        return M_elem_mat

    def getUpdateMatrix(Model, matrix, addval):
        elem_set = Model.element.getElementSet()
        shape_set = Model.shape.getShapeSet()
        dofe = len(shape_set["nodes"]) * len(elem_set["dofs"]["d"])
        for ii in range(len(addval)):
            
            A_add = addval[ii, 2] * array([[1, -1],
                                           [-1, 1]])
            
            loc = array([int(dofe * addval[ii, 0] - (dofe)),
                        int(dofe * addval[ii, 0]  - (dofe - 1)),])
           
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

    def getElementVolume(Model, inci, coord, tabgeo, intgauss, element_number):
        t = tabgeo[int(inci[element_number, 3] - 1)]["THICKN"] #tabgeo[int(inci[element_number, 3] - 1), 4]
        shape_set = Model.shape.getShapeSet()
        type_shape = shape_set["key"]
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        pt, wt = gauss_points(type_shape, intgauss)
        detJ = 0.0
        for pp in range(intgauss):
            detJ += abs(Model.shape.getdetJacobi(pt[pp], elementcoord))
        return detJ * t
