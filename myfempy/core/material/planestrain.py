import numpy as np

INT32 = np.uint32
FLT64 = np.float64

from myfempy.core.material.material import Material


class PlaneStrainIsotropic(Material):
    """Plane Strain Material Class <ConcreteClassService>"""

    def getMaterialSet():
        matset = {
            "mat": "planestrain",
            "idmat": 4,
            "type": "isotropic",
            "idtyp": 20,
        }
        return matset

    def getElasticTensor(E, v):
        D = np.zeros((3, 3), dtype=FLT64)
        D[0, 0] = E * (1.0 - v) / ((1 + v) * (1.0 - 2.0 * v))
        D[0, 1] = D[0, 0] * v / (1.0 - v)
        D[1, 0] = D[0, 1]
        D[1, 1] = D[0, 0]
        D[2, 2] = D[0, 0] * 0.5 * (1.0 - 2.0 * v) / (1.0 - v)
        return D

    def getElementStrain(Model, U, ptg, element_number):
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"])

        nodelist = Model.shape.getNodeList(Model.inci, element_number)

        loc = Model.shape.getShapeKey(nodelist, nodedof)

        elementcoord = Model.shape.getNodeCoord(Model.coord, nodelist)

        B = Model.element.getB(Model, elementcoord, ptg)

        epsilon = np.dot(B, U[loc])  # B @ (U[loc])

        strn_elm_xx = epsilon[0]
        strn_elm_yy = epsilon[1]
        strn_elm_xy = epsilon[2]
        strn_elm_vm = np.sqrt(
            epsilon[0] ** 2
            - epsilon[0] * epsilon[1]
            + epsilon[1] ** 2
            + 3 * epsilon[2] ** 2
        )

        strain = [strn_elm_vm, strn_elm_xx, strn_elm_yy, strn_elm_xy]

        return epsilon, strain

    def getTitleStrain():
        title = ["STRAIN_VM", "STRAIN_XX", "STRAIN_YY", "STRAIN_XY"]
        return title

    def getElementStress(Model, epsilon, element_number):
        E = Model.tabmat[
            int(Model.inci[element_number, 2]) - 1, 0
        ]  # material elasticity
        v = Model.tabmat[
            int(Model.inci[element_number, 2]) - 1, 1
        ]  # material poisson ratio

        C = PlaneStrainIsotropic.getElasticTensor(E, v)

        sigma = np.dot(C, epsilon)

        strs_elm_xx = sigma[0]
        strs_elm_yy = sigma[1]
        strs_elm_xy = sigma[2]

        strs_elm_vm = np.sqrt(
            sigma[0] ** 2 - sigma[0] * sigma[1] + sigma[1] ** 2 + 3 * sigma[2] ** 2
        )

        stress = [strs_elm_vm, strs_elm_xx, strs_elm_yy, strs_elm_xy]

        return stress

    def getTitleStress():
        title = ["STRESS_VM", "STRESS_XX", "STRESS_YY", "STRESS_XY"]
        return title
