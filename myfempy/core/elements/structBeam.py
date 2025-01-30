from __future__ import annotations

from os import environ

environ["OMP_NUM_THREADS"] = "1"

from numpy import (abs, array, concatenate, dot, eye, float64, int32, ix_,
                   sqrt, transpose, zeros)

INT32 = int32
FLT64 = float64

from myfempy.core.elements.element import Element
from myfempy.core.utilities import (gauss_points, get3D_LocalVector,
                                    getRotational_Matrix)


class StructuralBeam(Element):
    """Beam Structural Element Class <ConcreteClassService>"""

    def getElementSet():
        elemset = {
            "def": "1D-space 6-node_dofs",
            "key": "beam",
            "id": 16,
            "dofs": {
                "d": {"ux": 1, "uy": 2, "uz": 3, "rx": 4, "ry": 5, "rz": 6},
                "f": {
                    "fx": 1,
                    "fy": 2,
                    "fz": 3,
                    "tx": 4,
                    "ty": 5,
                    "tz": 6,
                    "massaadd": 15,
                    "spring2ground": 16,
                    "damper2ground": 17,
                },
            },
            "tensor": [
                "sntxx",
                "snbxymax",
                "snbxymin",
                "snbxzmax",
                "snbxzmin",
                "sstxy",
            ],
        }
        return elemset

    # def getH():
    #     return array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1],], dtype=INT32)

    def getB(Model, elementcoord, pt, nodedof):
        diffN = Model.shape.getDiffDiffShapeFuntion(pt, nodedof)
        invJ = Model.shape.getinvJacobi(pt, elementcoord, nodedof)
        H = array(
            [
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1],
            ],
            dtype=INT32,
        )  # StructuralBeam.getH()
        B = dot(H, dot(invJ, diffN))
        return B

    # @profile
    def getStifLinearMat(Model, inci, coord, tabmat, tabgeo, intgauss, element_number):
        elem_set = StructuralBeam.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        type_shape = shape_set["key"]
        edof = nodecon * nodedof
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        if type_shape == "line3":
            elementcoord_local = get3D_LocalVector(elementcoord, 3)
        else:
            elementcoord_local = get3D_LocalVector(array(elementcoord), 2)
        E = tabmat[int(inci[element_number, 2]) - 1]["EXX"]
        G = tabmat[int(inci[element_number, 2]) - 1]["GXX"]
        D = Model.material.getElasticTensor(E, G)
        AREA = tabgeo[int(inci[element_number, 3] - 1)][
            "AREACS"
        ]  # tabgeo[int(inci[element_number, 3] - 1), 4]
        IZZ = tabgeo[int(inci[element_number, 3] - 1)]["INERZZ"]
        IYY = tabgeo[int(inci[element_number, 3] - 1)]["INERYY"]
        IXX = tabgeo[int(inci[element_number, 3] - 1)]["INERXX"]
        C = array([[AREA, IZZ, IYY, IXX]]) * eye(4) * D
        pt, wt = gauss_points(type_shape, intgauss)
        K_elem_mat = zeros((edof, edof), dtype=FLT64)
        for ip in range(intgauss):
            detJ = Model.shape.getdetJacobi(array([pt[ip]]), elementcoord_local)
            B = StructuralBeam.getB(Model, elementcoord, array([pt[ip]]), nodedof)
            BCB = dot(dot(B.transpose(), C), B)
            K_elem_mat += BCB * abs(detJ) * wt[ip]
        if type_shape == "line3":
            R = getRotational_Matrix(elementcoord, 6)
        else:  # line2
            R = getRotational_Matrix(elementcoord, 4)
        K_elem_mat = dot(dot(transpose(R), K_elem_mat), R)
        return K_elem_mat

    def getMassConsistentMat(
        Model, inci, coord, tabmat, tabgeo, intgauss, element_number
    ):
        elem_set = StructuralBeam.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        type_shape = shape_set["key"]
        edof = nodecon * nodedof
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        if type_shape == "line3":
            elementcoord_local = get3D_LocalVector(elementcoord, 3)
        else:
            elementcoord_local = get3D_LocalVector(elementcoord, 2)
        rho = tabmat[int(inci[element_number, 2]) - 1][
            "RHO"
        ]  # tabmat[int(inci[element_number, 2]) - 1, 6]  # material density
        AREA = tabgeo[int(inci[element_number, 3] - 1)][
            "AREACS"
        ]  # tabgeo[int(inci[element_number, 3] - 1), 4]
        IXX = tabgeo[int(inci[element_number, 3] - 1)]["INERXX"]
        R = array([[rho * AREA, rho * AREA, rho * AREA, rho * IXX]]) * eye(4)
        pt, wt = gauss_points(type_shape, intgauss)
        M_elem_mat = zeros((edof, edof), dtype=FLT64)
        for ip in range(intgauss):
            detJ = Model.shape.getdetJacobi(array([pt[ip]]), elementcoord_local)
            N = Model.shape.getShapeFunctions(array([pt[ip]]), nodedof)
            NRN = dot(dot(N.transpose(), R), N)
            M_elem_mat += NRN * abs(detJ) * wt[ip]
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
            Udef[nn - 1, 0] = U[nodedof * nn - 6]
            Udef[nn - 1, 1] = U[nodedof * nn - 5]
            Udef[nn - 1, 2] = U[nodedof * nn - 4]
        return Udef

    def setTitleDeformation():
        return "DISPLACEMENT"

    def getElementVolume(Model, inci, coord, tabgeo, element_number):
        return 0.0
