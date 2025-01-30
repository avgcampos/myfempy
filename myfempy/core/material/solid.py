import numpy as np

INT32 = np.uint32
FLT64 = np.float64

from myfempy.core.material.material import Material


class SolidIsotropic(Material):
    """Solid Stress Isotropic Material Class <ConcreteClassService>"""

    def getMaterialSet():
        matset = {
            "mat": "solidstress",
            "type": "isotropic",
        }
        return matset

    # def getProMaterial(modeldata):
    #     tabmat = io.samethings()
    #     pass

    def getElasticTensor(E, v):
        D = np.zeros((6, 6))
        fac = 1.0 / (2.0 * v * v + v - 1.0)
        D[0, 0] = fac * E * (v - 1.0)
        D[0, 1] = -1.0 * fac * E * v
        D[0, 2] = D[0, 1]
        D[1, 0] = D[0, 1]
        D[1, 1] = D[0, 0]
        D[1, 2] = D[0, 1]
        D[2, 0] = D[0, 1]
        D[2, 1] = D[0, 1]
        D[2, 2] = D[0, 0]
        D[3, 3] = E / (2.0 + 2.0 * v)
        D[4, 4] = D[3, 3]
        D[5, 5] = D[3, 3]
        return D

    def getElementStrain(Model, U, ptg, element_number):
        elem_set = Model.element.getElementSet()
        nodedof = len(elem_set["dofs"]["d"])

        nodelist = Model.shape.getNodeList(Model.inci, element_number)

        loc = Model.shape.getShapeKey(nodelist, nodedof)

        elementcoord = Model.shape.getNodeCoord(Model.coord, nodelist)

        B = Model.element.getB(Model, elementcoord, ptg, nodedof)

        epsilon = np.dot(B, U[loc])  # B @ (U[loc])

        strn_elm_xx = epsilon[0]
        strn_elm_yy = epsilon[1]
        strn_elm_zz = epsilon[2]
        strn_elm_xy = epsilon[3]
        strn_elm_yz = epsilon[4]
        strn_elm_zx = epsilon[5]

        strn_elm_eqv = np.sqrt(
            0.5
            * (
                (epsilon[0] - epsilon[1]) ** 2
                + (epsilon[1] - epsilon[2]) ** 2
                + (epsilon[2] - epsilon[0]) ** 2
                + 6 * (epsilon[3] ** 2 + epsilon[4] ** 2 + epsilon[5] ** 2)
            )
        )

        strain = [
            strn_elm_eqv,
            strn_elm_xx,
            strn_elm_yy,
            strn_elm_zz,
            strn_elm_xy,
            strn_elm_yz,
            strn_elm_zx,
        ]

        return epsilon, strain

    def getTitleStrain():
        title = [
            "STRAIN_VM",
            "STRAIN_XX",
            "STRAIN_YY",
            "STRAIN_ZZ",
            "STRAIN_XY",
            "STRAIN_YZ",
            "STRAIN_ZX",
        ]
        return title

    def getElementStress(Model, epsilon, element_number):
        E = Model.tabmat[
            int(Model.inci[element_number, 2]) - 1, 0
        ]  # material elasticity
        v = Model.tabmat[
            int(Model.inci[element_number, 2]) - 1, 1
        ]  # material poisson ratio

        C = SolidIsotropic.getElasticTensor(E, v)

        sigma = np.dot(C, epsilon)

        strs_elm_xx = sigma[0]
        strs_elm_yy = sigma[1]
        strs_elm_zz = sigma[2]
        strs_elm_xy = sigma[3]
        strs_elm_yz = sigma[4]
        strs_elm_zx = sigma[5]
        strs_elm_eqv = np.sqrt(
            0.5
            * (
                (sigma[0] - sigma[1]) ** 2
                + (sigma[1] - sigma[2]) ** 2
                + (sigma[2] - sigma[0]) ** 2
                + 6 * (sigma[3] ** 2 + sigma[4] ** 2 + sigma[5] ** 2)
            )
        )

        stress = [
            strs_elm_eqv,
            strs_elm_xx,
            strs_elm_yy,
            strs_elm_zz,
            strs_elm_xy,
            strs_elm_yz,
            strs_elm_zx,
        ]

        return sigma, stress

    def getTitleStress():
        title = [
            "STRESS_VM",
            "STRESS_XX",
            "STRESS_YY",
            "STRESS_ZZ",
            "STRESS_XY",
            "STRESS_YZ",
            "STRESS_ZX",
        ]
        return title

    def getStrainEnergyDensity(sigma, epsilon, elemvol):
        strain_energy = 0.5 * np.dot(sigma.transpose(), epsilon) / elemvol
        return strain_energy

    def getTitleCompliance():
        title = ["STRAIN_ENERGY_DENSITY"]
        return title

    def getFailureCriteria(sigma):
        return 0.0

    def getTitleFoS():
        title = ["FoS_YIELD_VON_MISES"]
        return title
