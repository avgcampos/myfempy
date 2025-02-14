from __future__ import annotations

from numpy import (abs, array, array2string, concatenate, dot, float64, in1d,
                   int32, ix_, sqrt, unique, where, zeros)

from myfempy.core.utilities import (gauss_points, get_elemen_from_nodelist,
                                    get_nodes_from_list)

INT32 = int32
FLT64 = float64

from myfempy.core.elements.element import Element
from myfempy.core.utilities import gauss_points


class HeatPlane(Element):
    """Plane Heat Element Class <ConcreteClassService>"""

    def getElementSet():
        elemset = {
            "def": "2D-space 1-node_dofs",
            "key": "plane",
            "id": 21,
            "dofs": {
                "d": {"t": 1},
                "f": {"heatflux": 1, "convection": 15},
            },
            "tensor": ["qxx", "qyy"],
        }
        return elemset

    def getB(diffN, invJ):
        H = array([[1, 0], [0, 1]], dtype=FLT64)
        B = H.dot(invJ).dot(diffN)
        return B

    # @profile
    def getStifLinearMat(Model, inci, coord, tabmat, tabgeo, intgauss, element_number):
        elem_set = HeatPlane.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        type_shape = shape_set["key"]
        edof = nodecon * nodedof
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        C = Model.material.getElasticTensor(tabmat, inci, element_number)
        t = tabgeo[int(inci[element_number, 3] - 1)]["THICKN"]
        pt, wt = gauss_points(type_shape, intgauss)
        K_elem_mat = zeros((edof, edof), dtype=FLT64)
        for ip in range(intgauss):
            for jp in range(intgauss):
                detJ = Model.shape.getdetJacobi(array([pt[ip], pt[jp]]), elementcoord)
                diffN = Model.shape.getDiffShapeFuntion(array([pt[ip], pt[jp]]), nodedof)
                invJ = Model.shape.getinvJacobi(array([pt[ip], pt[jp]]), elementcoord, nodedof)
                B = HeatPlane.getB(diffN, invJ)
                BCB = B.transpose().dot(C).dot(B)
                K_elem_mat += BCB * t * abs(detJ) * wt[ip] * wt[jp]
        return K_elem_mat

    # def getMassConsistentMat(
    #     Model, inci, coord, tabmat, tabgeo, intgauss, element_number
    # ):
    #     elem_set = HeatPlane.getElementSet()
    #     nodedof = len(elem_set["dofs"]["d"])
    #     shape_set = Model.shape.getShapeSet()
    #     nodecon = len(shape_set["nodes"])
    #     type_shape = shape_set["key"]
    #     edof = nodecon * nodedof
    #     nodelist = Model.shape.getNodeList(inci, element_number)
    #     elementcoord = Model.shape.getNodeCoord(coord, nodelist)
    #     R = tabmat[int(inci[element_number, 2]) - 1, 6]  # material density
    #     t = tabgeo[int(inci[element_number, 3] - 1), 4]
    #     pt, wt = gauss_points(type_shape, intgauss)
    #     M_elem_mat = zeros((edof, edof), dtype=FLT64)
    #     for pp in range(intgauss):
    #         detJ = Model.shape.getdetJacobi(pt[pp], elementcoord)
    #         N = Model.shape.getShapeFunctions(pt[pp], nodedof)
    #         NRN = NTRN(N, R)
    #         M_elem_mat += NRN * t * abs(detJ) * wt[pp]
    #     return M_elem_mat

    def getUpdateMatrix(Model, matrix, addval):
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])
        shape_set = Model.shape.getShapeSet()
        nodecon = len(shape_set["nodes"])
        type_shape = shape_set["key"]
        edof = nodecon * nodedof

        nodelistconv = unique(addval[:, 0])
        elmlist = get_elemen_from_nodelist(Model.inci, nodelistconv)
        for ee in range(len(elmlist)):

            nodelist = Model.shape.getNodeList(Model.inci, elmlist[ee] - 1)
            elementcoord = Model.shape.getNodeCoord(Model.coord, nodelist)
            t = Model.tabgeo[int(Model.inci[elmlist[ee] - 1, 3] - 1)]["THICKN"] 
            test = in1d(nodelist, nodelistconv, assume_unique=True)

            # nodes = array(nodelist)[test]
            nodes_conec = where(test == True)[0]
            if len(nodes_conec) < 2:
                pass
            else:
                idx_conec = array2string(nodes_conec)
                get_side = Model.shape.getSideAxis(nodes_conec[1:-1])
                pt, wt = gauss_points(type_shape, Model.intgauss)
                loc = Model.shape.getLocKey(nodelist, nodedof)
                h = addval[0, 2]
                Kh = zeros((edof, edof))
                for ip in range(2):
                    for jp in range(2):
                        points = Model.shape.getIsoParaSide(get_side, pt[ip])
                        N = Model.shape.getShapeFunctions(array(points), nodedof)
                        J = Model.shape.getJacobian(array(points), elementcoord)
                        detJ_e = Model.shape.getEdgeLength(J, get_side)
                        Kh += dot(N.transpose(), N) * h * t * abs(detJ_e) * wt[ip] * wt[jp]

                matrix[ix_(loc, loc)] += Kh
        return matrix

    def getElementDeformation(U, modelinfo):
        nodetot = modelinfo["nnode"]
        nodedof = modelinfo["nodedof"]
        Udef = zeros((nodetot, 1), dtype=FLT64)
        for nn in range(1, nodetot + 1):
            Udef[nn - 1, 0] = U[nodedof * nn - 1]
        return Udef

    def setTitleDeformation():
        return "TEMPERATURE"

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
