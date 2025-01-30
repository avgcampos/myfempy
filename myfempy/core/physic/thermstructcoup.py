from __future__ import annotations

import numpy as np
from scipy.special import roots_legendre

from myfempy.core.physic.loadstruct import LoadStructural
from myfempy.core.physic.structural import Structural
from myfempy.core.utilities import (gauss_points, get_elemen_from_nodelist,
                                    get_nodes_from_list, poly_area)


class ThermalStructuralCoupling(Structural):
    """Thermal Structural Coupled field analysis <ConcreteClassService>"""

    def getLoadApply(Model, modelinfo, coupling):
        forcenodeaply = np.zeros((1, 4))
        if coupling["TYPE"] == "thermalstress":  # thermo stress mechanical
            fapp = ThermalStructuralCoupling.__ForceThermalStress(
                Model, modelinfo, coupling
            )
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        else:
            pass
        forcenodeaply = forcenodeaply[1::][::]
        return forcenodeaply

    def __ForceThermalStress(Model, modelinfo, coupling):
        forcenodedof = np.zeros((1, 4))
        inci = modelinfo["inci"]
        coord = modelinfo["coord"]
        tabmat = modelinfo["tabmat"]
        tabgeo = modelinfo["tabgeo"]
        intgauss = modelinfo["intgauss"]
        # fc_type_dof = modelinfo["dofs"]["f"][forcelist['DOF']]

        strain_thermal = np.zeros((inci.shape[0], 3))

        dT = coupling["GRADTEMP"]  # np.mean(coupling['TEMPERATURE'])
        strain_thermal[:, 0] = dT
        strain_thermal[:, 1] = dT

        for ee in range(inci.shape[0]):
            force_value_vector, nodelist = (
                ThermalStructuralCoupling.__body_thermal_stress(
                    Model,
                    inci,
                    coord,
                    tabmat,
                    tabgeo,
                    intgauss,
                    strain_thermal[ee, :],
                    ee,
                )
            )

            fc_type_dof = np.tile(
                [modelinfo["dofs"]["f"]["fx"], modelinfo["dofs"]["f"]["fy"]],
                len(nodelist),
            )
            nodelist = np.repeat(nodelist, 2)
            for j in range(len(nodelist)):
                fcdof = np.array(
                    [
                        [
                            int(nodelist[j]),
                            fc_type_dof[j],
                            force_value_vector[j],
                            int(coupling["STEP"]),
                        ]
                    ]
                )
                forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        # forcenodedof[np.nonzero(forcenodedof)]
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    def getUpdateMatrix(Model, matrix, modelinfo, loadaply):
        return LoadStructural.getUpdateMatrix(Model, matrix, modelinfo, loadaply)

    def getUpdateLoad(self):
        return LoadStructural.getUpdateLoad(self)

    def __body_thermal_stress(
        Model,
        inci,
        coord,
        tabmat,
        tabgeo,
        intgauss,
        strain_thermal,
        element_number,
    ):
        # body force
        elem_set = Model.element.getElementSet()
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
        a = tabmat[int(inci[element_number, 2]) - 1]["CTE"]
        t = tabgeo[int(inci[element_number, 3] - 1)][
            "THICKN"
        ]  # tabgeo[int(inci[element_number, 3] - 1), 4]
        C = Model.material.getElasticTensor(E, v)
        pt, wt = gauss_points(type_shape, intgauss)
        W = np.zeros((nodedof, 1))
        force_value_vector = np.zeros((edof, 1))
        for ip in range(intgauss):
            for jp in range(intgauss):
                detJ = Model.shape.getdetJacobi(
                    np.array([pt[ip], pt[jp]]), elementcoord
                )
                B = Model.element.getB(
                    Model, elementcoord, np.array([pt[ip], pt[jp]]), nodedof
                )
                force_value_vector += np.dot(
                    np.dot(B.transpose(), C), np.reshape(strain_thermal, (3, 1))
                ) * (a * t * abs(detJ) * wt[ip] * wt[jp])
        force_value_vector = np.reshape(force_value_vector, (edof))
        return force_value_vector, nodelist
