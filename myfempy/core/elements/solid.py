from __future__ import annotations

import numdifftools as nd
import numpy as np

from myfempy.core.elements.element import Element
# from myfempy.felib.material.setpropmat import getElasticity
from myfempy.core.utilities import gauss_points

# from myfempy.felib.quadrature import gaussian, no_interpol
# from myfempy.felib.materset import get_elasticity


class Solid(Element):
    """Solid Structural Element Class <ConcreteClassService>"""

    def getElementSet():
        elemset = {
            "def": "3D-space 3-node_dofs",
            "key": "solid",
            "id": 33,
            "dofs": {
                "d": {"ux": 1, "uy": 2, "uz": 3},
                "f": {"fx": 1, "fy": 2, "fz": 3},
            },
            "tensor": ["sxx", "syy", "szz", "sxy", "syz", "szx"],
        }
        return elemset

    def getL():
        L = np.array(
            [
                [1, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 1, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 1],
                [0, 1, 0, 1, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 1, 0, 1, 0],
                [0, 0, 1, 0, 0, 0, 1, 0, 0],
            ]
        )
        return L

    def getB(Model, elementcoord, ptg, nodedof):
        diffN = Model.shape.getDiffShapeFuntion(Model.shape.N, ptg, nodedof)
        invJ = Model.shape.invJacobi(Model.shape.N, ptg, elementcoord, nodedof)

        return np.dot(Solid.getL(), np.dot(invJ, diffN))

    def getStifLinearMat(Model, inci, coord, tabmat, tabgeo, intgauss, element_number):
        elem_set = Solid.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])

        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        type_shape = shape_set["key"]

        # nelem = len(inci)
        # nnode = len(coord)

        # fulldof = nodedof * nnode
        edof = nodecon * nodedof

        nodelist = Model.shape.getNodeList(inci, element_number)

        elementcoord = Model.shape.getNodeCoord(coord, nodelist)

        E = tabmat[int(inci[element_number, 2]) - 1, 0]  # material elasticity
        v = tabmat[int(inci[element_number, 2]) - 1, 1]  # material poisson ratio
        C = Model.material.getElasticTensor(E, v)

        pt, wt = gauss_points(type_shape, intgauss)

        K_elem_mat = np.zeros((edof, edof))
        for pp in range(intgauss):
            detJ = Model.shape.detJacobi(Model.shape.N, pt[pp], elementcoord)

            B = Solid.getB(Model, elementcoord, pt[pp], nodedof)

            K_elem_mat += (
                np.dot(np.dot(np.transpose(B), C), B) * detJ * wt[pp] * wt[pp] * wt[pp]
            )

        return K_elem_mat

    def getMassConsistentMat(
        Model, inci, coord, tabmat, tabgeo, intgauss, element_number
    ):
        elem_set = Solid.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])

        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        type_shape = shape_set["key"]

        edof = nodecon * nodedof

        nodelist = Model.shape.getNodeList(inci, element_number)

        elementcoord = Model.shape.getNodeCoord(coord, nodelist)

        R = tabmat[int(inci[element_number, 2]) - 1, 6]  # material density

        pt, wt = gauss_points(type_shape, intgauss)

        M_elem_mat = np.zeros((edof, edof))
        for pp in range(intgauss):
            detJ = Model.shape.detJacobi(Model.shape.N, pt[pp], elementcoord)
            matN = Model.shape.getShapeFunctions(pt[pp], nodedof)
            M_elem_mat += (
                np.dot(np.dot(np.transpose(matN), R), matN)
                * detJ
                * wt[pp]
                * wt[pp]
                * wt[pp]
            )

        return M_elem_mat

    def getElementDeformation(U, modelinfo):
        nodetot = modelinfo["nnode"]
        nodedof = modelinfo["nodedof"]

        Udef = np.zeros((nodetot, 3), dtype=float)
        Umag = np.zeros((nodetot, 1), dtype=float)
        for nn in range(1, nodetot + 1):
            Udef[nn - 1, 0] = U[nodedof * nn - 3]
            Udef[nn - 1, 1] = U[nodedof * nn - 2]
            Udef[nn - 1, 2] = U[nodedof * nn - 1]
            Umag[nn - 1, 0] = np.sqrt(
                U[nodedof * nn - 3] ** 2
                + U[nodedof * nn - 2] ** 2
                + U[nodedof * nn - 1] ** 2
            )
        result = np.concatenate((Umag, Udef), axis=1)

        return result

    def getTitleDeformation():
        title = ["DISPL_X", "DISPL_Y", "DISPL_Z"]
        return title

    def getElementVolume(Model, inci, coord, tabgeo, intgauss, element_number):
        shape_set = Model.shape.getShapeSet()
        type_shape = shape_set["key"]
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        pt, wt = gauss_points(type_shape, intgauss)
        detJ = 0.0
        for pp in range(intgauss):
            detJ += Model.shape.detJacobi(Model.shape.N, pt[pp], elementcoord)
        return detJ
