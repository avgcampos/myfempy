from __future__ import annotations

import numpy as np
from scipy.special import roots_legendre

from myfempy.core.physic.physics import Physics
from myfempy.core.utilities import (gauss_points, get3D_LocalVector,
                                    get_elemen_from_nodelist,
                                    get_nodes_from_list, getRotational_Matrix)

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


class LoadStructural(Physics):
    """Structural Load Class <ConcreteClassService>"""

    def getLoadApply(Model, forcelist):
        forcenodeaply = np.zeros((1, 4))
        if forcelist["TYPE"] == "forcenode":
            fapp = LoadStructural.ForceNodeLoadApply(Model, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        elif forcelist["TYPE"] == "forceedge":
            fapp = LoadStructural.ForceEdgeLoadApply(Model, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        elif forcelist["TYPE"] == "forcesurf":
            fapp = LoadStructural.ForceSurfLoadApply(Model, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        elif forcelist["TYPE"] == "forcebeam":
            fapp = LoadStructural.ForceBeamLoadApply(Model, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        elif forcelist["TYPE"] == "bodyforce":
            fapp = LoadStructural.BodyForce(Model, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        elif forcelist["TYPE"] == "strainzero":
            fapp = LoadStructural.ForceStrainZero(Model, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        else:
            pass
        forcenodeaply = forcenodeaply[1::][::]
        return forcenodeaply

    def getUpdateMatrix(Model, matrix, loadaply):

        addSpring = np.where(loadaply[:, 1] == 16)
        addMass = np.where(loadaply[:, 1] == 15)

        if addSpring[0].size:
            addLoad = loadaply[addSpring, :][0]
            matrix["stiffness"] = Model.element.getUpdateMatrix(
                Model, matrix["stiffness"], addLoad
            )

        if addMass[0].size:
            addLoad = loadaply[addMass, :][0]
            matrix["mass"] = Model.element.getUpdateMatrix(
                Model, matrix["mass"], addLoad
            )

        return matrix

    def getUpdateLoad(self):
        return None

    def ForceStrainZero(Model, forcelist):
        forcenodedof = np.zeros((1, 4))
        strain_zero = np.array(forcelist["VAL"])
        inci = Model.inci
        coord = Model.coord
        tabmat = Model.tabmat
        tabgeo = Model.tabgeo
        intgauss = Model.intgauss
        for ee in range(inci.shape[0]):
            if (
                Model.modelinfo["shape"] == "hexa8"
                or Model.modelinfo["shape"] == "tetr4"
            ):
                force_value_vector, nodelist = LoadStructural.__solid_strain_zero(
                    Model, inci, coord, tabmat, tabgeo, intgauss, ee, strain_zero
                )
                fc_type_dof = np.tile(
                    [
                        Model.modelinfo["dofs"]["f"]["fx"],
                        Model.modelinfo["dofs"]["f"]["fy"],
                        Model.modelinfo["dofs"]["f"]["fz"],
                    ],
                    len(nodelist),
                )
                nodelist = np.repeat(nodelist, 3)

            elif (
                Model.modelinfo["shape"] == "tria3"
                or Model.modelinfo["shape"] == "tria6"
                or Model.modelinfo["shape"] == "quad4"
                or Model.modelinfo["shape"] == "quad8"
            ):
                force_value_vector, nodelist = LoadStructural.__plane_strain_zero(
                    Model, inci, coord, tabmat, tabgeo, intgauss, ee, strain_zero
                )
                fc_type_dof = np.tile(
                    [
                        Model.modelinfo["dofs"]["f"]["fx"],
                        Model.modelinfo["dofs"]["f"]["fy"],
                    ],
                    len(nodelist),
                )
                nodelist = np.repeat(nodelist, 2)

            else:
                force_value_vector, nodelist = [0], [0]

            for j in range(len(nodelist)):
                fcdof = np.array(
                    [
                        [
                            int(nodelist[j]),
                            fc_type_dof[j],
                            force_value_vector[j],
                            int(forcelist["STEP"]),
                        ]
                    ]
                )
                forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    def ForceNodeLoadApply(Model, forcelist):
        forcenodedof = np.zeros((1, 4))
        nodelist = [
            forcelist["DIR"],
            forcelist["LOCX"],
            forcelist["LOCY"],
            forcelist["LOCZ"],
            forcelist["TAG"],
            forcelist["MESHNODE"],
        ]
        node_list_fc, dir_fc = get_nodes_from_list(nodelist, Model.coord, Model.regions)
        force_value_vector = np.ones_like(node_list_fc, dtype=float) * float(forcelist["VAL"])
        fc_type_dof = Model.modelinfo["dofs"]["f"][forcelist["DOF"]] * np.ones_like(
            node_list_fc
        )
        for j in range(len(node_list_fc)):
            fcdof = np.array(
                [
                    [
                        int(node_list_fc[j]),
                        fc_type_dof[j],
                        force_value_vector[j],
                        int(forcelist["STEP"]),
                    ]
                ]
            )
            forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    def BodyForce(Model, forcelist):
        forcenodedof = np.zeros((1, 4))
        gravity_value = float(forcelist["VAL"])
        inci = Model.inci
        coord = Model.coord
        tabmat = Model.tabmat
        tabgeo = Model.tabgeo
        intgauss = Model.intgauss
        fc_type_dof = Model.modelinfo["dofs"]["f"][forcelist["DOF"]]
        for ee in range(inci.shape[0]):

            if (
                Model.modelinfo["shape"] == "hexa8"
                or Model.modelinfo["shape"] == "tetr4"
            ):
                force_value_vector, nodelist = (
                    LoadStructural.__solid_body_force_volumetric(
                        Model,
                        inci,
                        coord,
                        tabmat,
                        tabgeo,
                        intgauss,
                        ee,
                        gravity_value,
                        fc_type_dof,
                    )
                )
            elif (
                Model.modelinfo["shape"] == "tria3"
                or Model.modelinfo["shape"] == "tria6"
                or Model.modelinfo["shape"] == "quad4"
                or Model.modelinfo["shape"] == "quad8"
            ):
                force_value_vector, nodelist = (
                    LoadStructural.__plane_body_force_volumetric(
                        Model,
                        inci,
                        coord,
                        tabmat,
                        tabgeo,
                        intgauss,
                        ee,
                        gravity_value,
                        fc_type_dof,
                    )
                )
            else:
                force_value_vector, nodelist = [0.0], [0]

            for j in range(len(nodelist)):
                fcdof = np.array(
                    [
                        [
                            int(nodelist[j]),
                            fc_type_dof,
                            force_value_vector[j],
                            int(forcelist["STEP"]),
                        ]
                    ]
                )
                forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    def ForceEdgeLoadApply(Model, forcelist):
        nodelist = [
            forcelist["DIR"],
            forcelist["LOCX"],
            forcelist["LOCY"],
            forcelist["LOCZ"],
            forcelist["TAG"],
            forcelist["MESHNODE"],
        ]  # forcelist[3:]
        node_list_fc, dir_fc = get_nodes_from_list(nodelist, Model.coord, Model.regions)
        force_value = float(forcelist["VAL"])
        force_dirc = forcelist["DOF"]
        inci = Model.inci
        coord = Model.coord
        tabmat = Model.tabmat
        tabgeo = Model.tabgeo
        intgauss = Model.intgauss
        fc_type = forcelist["DOF"]
        elmlist = get_elemen_from_nodelist(inci, node_list_fc)
        forcenodedof = np.zeros((1, 4))
        for ee in range(len(elmlist)):
            force_value_vector, nodelist, norm = (
                LoadStructural.__line_force_distribuition(
                    Model,
                    inci,
                    coord,
                    tabgeo,
                    intgauss,
                    node_list_fc,
                    elmlist[ee],
                    force_value,
                    fc_type,
                )
            )

            elem_set = Model.element.getElementSet()
            nodedof = len(elem_set["dofs"]["d"])

            if len(force_value_vector) > len(nodelist):
                nodelist = np.repeat(nodelist, nodedof)

            if forcelist["DOF"] == "pressure":
                shape_set = Model.shape.getShapeSet()
                nodesce = shape_set["nodesconecedge"]
                if int(norm[0]) == 1 and int(norm[1]) == 0:
                    fc_type_dof = Model.modelinfo["dofs"]["f"]["fx"] * np.ones_like(
                        nodelist
                    )
                elif int(norm[0]) == 0 and int(norm[1]) == 1:
                    fc_type_dof = Model.modelinfo["dofs"]["f"]["fy"] * np.ones_like(
                        nodelist
                    )
                else:
                    fc_type_dof = np.tile(
                        [
                            Model.modelinfo["dofs"]["f"]["fx"],
                            Model.modelinfo["dofs"]["f"]["fy"],
                        ],
                        nodesce,
                    )
            else:
                fc_type_dof = Model.modelinfo["dofs"]["f"][
                    forcelist["DOF"]
                ] * np.ones_like(nodelist)

            for j in range(len(nodelist)):
                fcdof = np.array(
                    [
                        [
                            int(nodelist[j]),
                            fc_type_dof[j],
                            force_value_vector[j],
                            int(forcelist["STEP"]),
                        ]
                    ]
                )
                forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    def ForceSurfLoadApply(Model, forcelist):
        nodelist = [
            forcelist["DIR"],
            forcelist["LOCX"],
            forcelist["LOCY"],
            forcelist["LOCZ"],
            forcelist["TAG"],
            forcelist["MESHNODE"],
        ]  # forcelist[3:]
        node_list_fc, dir_fc = get_nodes_from_list(nodelist, Model.coord, Model.regions)
        force_value = float(forcelist["VAL"])
        force_dirc = forcelist["DOF"]
        inci = Model.inci
        coord = Model.coord
        tabmat = Model.tabmat
        tabgeo = Model.tabgeo
        intgauss = Model.intgauss
        fc_type = forcelist["DOF"]
        elmlist = get_elemen_from_nodelist(inci, node_list_fc)
        forcenodedof = np.zeros((1, 4))
        for ee in range(len(elmlist)):
            force_value_vector, nodelist, norm = (
                LoadStructural.__surf_force_distribuition(
                    Model,
                    inci,
                    coord,
                    tabgeo,
                    intgauss,
                    node_list_fc,
                    elmlist[ee],
                    force_value,
                    fc_type,
                )
            )

            elem_set = Model.element.getElementSet()
            nodedof = len(elem_set["dofs"]["d"])

            if len(force_value_vector) > len(nodelist):
                nodelist = np.repeat(nodelist, int(len(force_value_vector) / nodedof))

            if forcelist["DOF"] == "pressure":
                shape_set = Model.shape.getShapeSet()
                nodescf = shape_set["nodesconecface"]
                if int(norm[0]) == 1 and int(norm[1]) == 0 and int(norm[2]) == 0:
                    fc_type_dof = Model.modelinfo["dofs"]["f"]["fx"] * np.ones_like(
                        nodelist
                    )
                elif int(norm[0]) == 0 and int(norm[1]) == 1 and int(norm[2]) == 0:
                    fc_type_dof = Model.modelinfo["dofs"]["f"]["fy"] * np.ones_like(
                        nodelist
                    )
                elif int(norm[0]) == 0 and int(norm[1]) == 0 and int(norm[2]) == 1:
                    fc_type_dof = Model.modelinfo["dofs"]["f"]["fz"] * np.ones_like(
                        nodelist
                    )
                elif int(norm[0]) == 1 and int(norm[1]) == 1 and int(norm[2]) == 0:
                    fc_type_dof = np.tile(
                        [
                            Model.modelinfo["dofs"]["f"]["fx"],
                            Model.modelinfo["dofs"]["f"]["fy"],
                        ],
                        nodescf,
                    )
                elif int(norm[0]) == 1 and int(norm[1]) == 0 and int(norm[2]) == 1:
                    fc_type_dof = np.tile(
                        [
                            Model.modelinfo["dofs"]["f"]["fx"],
                            Model.modelinfo["dofs"]["f"]["fz"],
                        ],
                        nodescf,
                    )
                elif int(norm[0]) == 0 and int(norm[1]) == 1 and int(norm[2]) == 1:
                    fc_type_dof = np.tile(
                        [
                            Model.modelinfo["dofs"]["f"]["fy"],
                            Model.modelinfo["dofs"]["f"]["fz"],
                        ],
                        nodescf,
                    )
                else:
                    fc_type_dof = np.tile(
                        [
                            Model.modelinfo["dofs"]["f"]["fx"],
                            Model.modelinfo["dofs"]["f"]["fy"],
                            Model.modelinfo["dofs"]["f"]["fz"],
                        ],
                        nodescf,
                    )
            else:
                fc_type_dof = Model.modelinfo["dofs"]["f"][
                    forcelist["DOF"]
                ] * np.ones_like(nodelist)

            for j in range(len(nodelist)):
                fcdof = np.array(
                    [
                        [
                            int(nodelist[j]),
                            fc_type_dof[j],
                            force_value_vector[j],
                            int(forcelist["STEP"]),
                        ]
                    ]
                )
                forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    def ForceBeamLoadApply(Model, forcelist):
        nodelist = [
            forcelist["DIR"],
            forcelist["LOCX"],
            forcelist["LOCY"],
            forcelist["LOCZ"],
            forcelist["TAG"],
            forcelist["MESHNODE"],
        ]  # forcelist[3:]
        node_list_fc, dir_fc = get_nodes_from_list(nodelist, Model.coord, Model.regions)
        force_value = float(forcelist["VAL"])
        force_dirc = forcelist["DOF"]
        inci = Model.inci
        coord = Model.coord
        tabmat = Model.tabmat
        tabgeo = Model.tabgeo
        intgauss = Model.intgauss
        fc_type = forcelist["DOF"]
        elmlist = get_elemen_from_nodelist(inci, node_list_fc)
        forcenodedof = np.zeros((1, 4))
        for ee in range(len(elmlist)):
            force_value_vector, nodeslist, norm = (
                LoadStructural.__line_beam_distribuition(
                    Model,
                    inci,
                    coord,
                    tabgeo,
                    intgauss,
                    node_list_fc,
                    elmlist[ee],
                    force_value,
                    fc_type,
                )
            )

            nodedof = Model.modelinfo["nodedof"]

            nodes = nodeslist
            if len(force_value_vector) > len(nodeslist):
                nodeslist = np.repeat(nodeslist, 2)

            if forcelist["DOF"] == "pressure":
                if int(norm[0]) == 1 and int(norm[1]) == 0:
                    fc_type_dof = Model.modelinfo["dofs"]["f"]["fx"] * np.ones_like(
                        nodeslist
                    )
                elif int(norm[0]) == 0 and int(norm[1]) == 1:
                    fc_type_dof = Model.modelinfo["dofs"]["f"]["fy"] * np.ones_like(
                        nodeslist
                    )
                else:
                    fc_type_dof = np.tile(
                        [
                            Model.modelinfo["dofs"]["f"]["fx"],
                            Model.modelinfo["dofs"]["f"]["fy"],
                        ],
                        int(len(nodes) / nodedof),
                    )
            else:
                if force_dirc == "fy":
                    fc_type_dof = np.tile([2, 6], len(nodes))
                elif force_dirc == "fz":
                    fc_type_dof = np.tile([3, 5], len(nodes))
                elif force_dirc == "fx":
                    fc_type_dof = np.tile([1], len(nodes))
                elif force_dirc == "tx":
                    fc_type_dof = np.tile([4], len(nodes))
                else:
                    pass

            for j in range(len(nodeslist)):
                fcdof = np.array(
                    [
                        [
                            int(nodeslist[j]),
                            fc_type_dof[j],
                            force_value_vector[j],
                            int(forcelist["STEP"]),
                        ]
                    ]
                )
                forcenodedof = np.append(forcenodedof, fcdof, axis=0)
        forcenodedof = forcenodedof[1::][::]
        return forcenodedof

    def __plane_strain_zero(
        Model,
        inci,
        coord,
        tabmat,
        tabgeo,
        intgauss,
        element_number,
        strain_zero,
    ):
        # internal balance forces
        elem_set = Model.element.getElementSet()
        H = elem_set["H"]
        nodedof = Model.modelinfo["nodedof"]
        shape = Model.modelinfo["shape"]
        edof = Model.modelinfo["elemdof"]
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        C = Model.material.getElasticTensor(tabmat, inci, element_number)
        pt, wt = gauss_points(shape, intgauss)
        len_sigma = len(elem_set["tensor"])
        force_value_vector = np.zeros((edof, 1))
        for ip in range(intgauss):
            for jp in range(intgauss):
                detJ = Model.shape.getdetJacobi(
                    np.array([pt[ip], pt[jp]]), elementcoord
                )
                diffN = Model.shape.getDiffShapeFuntion(
                    np.array([pt[ip], pt[jp]]), nodedof
                )
                invJ = Model.shape.getinvJacobi(
                    np.array([pt[ip], pt[jp]]), elementcoord, nodedof
                )
                B = Model.shape.getB(H, invJ, diffN)
                force_value_vector += (
                    np.dot(np.dot(B.transpose(), C), strain_zero.reshape(-1, 1))
                    * abs(detJ)
                    * wt[ip]
                    * wt[jp]
                )
        force_value_vector = np.reshape(force_value_vector, (edof))
        return force_value_vector, nodelist

    def __solid_strain_zero(
        Model,
        inci,
        coord,
        tabmat,
        tabgeo,
        intgauss,
        element_number,
        strain_zero,
    ):
        # internal balance forces
        elem_set = Model.element.getElementSet()
        nodedof = Model.modelinfo["nodedof"]
        shape = Model.modelinfo["shape"]
        edof = Model.modelinfo["elemdof"]
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        C = Model.material.getElasticTensor(tabmat, inci, element_number)
        pt, wt = gauss_points(shape, intgauss)
        force_value_vector = np.zeros((edof, 1))
        for ip in range(intgauss):
            for jp in range(intgauss):
                for kp in range(intgauss):
                    detJ = Model.shape.getdetJacobi(
                        np.array([pt[ip], pt[jp], pt[kp]]), elementcoord
                    )
                    diffN = Model.shape.getDiffShapeFuntion(
                        np.array([pt[ip], pt[jp], pt[kp]]), nodedof
                    )
                    invJ = Model.shape.getinvJacobi(
                        np.array([pt[ip], pt[jp], pt[kp]]), elementcoord, nodedof
                    )
                    B = Model.element.getB(diffN, invJ)
                    force_value_vector += (
                        np.dot(np.dot(B.transpose(), C), strain_zero.reshape(-1, 1))
                        * abs(detJ)
                        * wt[ip]
                        * wt[jp]
                        * wt[kp]
                    )
        force_value_vector = np.reshape(force_value_vector, (edof))
        return force_value_vector, nodelist

    def __plane_body_force_volumetric(
        Model,
        inci,
        coord,
        tabmat,
        tabgeo,
        intgauss,
        element_number,
        gravity_value,
        fc_type_dof,
    ):
        # body force plane
        nodedof = Model.modelinfo["nodedof"]
        shape = Model.modelinfo["shape"]
        edof = Model.modelinfo["elemdof"]
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        R = tabmat[int(inci[element_number, 2]) - 1]["RHO"]
        t = tabgeo[int(inci[element_number, 3] - 1)]["THICKN"]
        pt, wt = gauss_points(shape, intgauss)
        G = gravity_value
        W = np.zeros((nodedof, 1))
        W[fc_type_dof - 1, 0] = R * G
        force_value_vector = np.zeros((edof, 1))
        for ip in range(intgauss):
            for jp in range(intgauss):
                detJ = Model.shape.getdetJacobi(
                    np.array([pt[ip], pt[jp]]), elementcoord
                )
                N = Model.shape.getShapeFunctions(np.array([pt[ip], pt[jp]]), nodedof)
                force_value_vector += (
                    np.dot(N.transpose(), W) * t * abs(detJ) * wt[ip] * wt[jp]
                )
        # Mantém os valores (inclusive 0.0) apenas nos DOFs correspondentes a fc_type_dof
        force_value_vector = force_value_vector[fc_type_dof - 1 :: nodedof, 0]
        return force_value_vector, nodelist

    def __solid_body_force_volumetric(
        Model,
        inci,
        coord,
        tabmat,
        tabgeo,
        intgauss,
        element_number,
        gravity_value,
        fc_type_dof,
    ):
        # body force solid
        nodedof = Model.modelinfo["nodedof"]
        shape = Model.modelinfo["shape"]
        edof = Model.modelinfo["elemdof"]
        nodelist = Model.shape.getNodeList(inci, element_number)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        R = tabmat[int(inci[element_number, 2]) - 1]["RHO"]
        pt, wt = gauss_points(shape, intgauss)
        G = gravity_value
        W = np.zeros((nodedof, 1))
        W[fc_type_dof - 1, 0] = R * G
        force_value_vector = np.zeros((edof, 1))
        for ip in range(intgauss):
            for jp in range(intgauss):
                for kp in range(intgauss):
                    detJ = Model.shape.getdetJacobi(
                        np.array([pt[ip], pt[jp], pt[kp]]), elementcoord
                    )
                    N = Model.shape.getShapeFunctions(
                        np.array([pt[ip], pt[jp], pt[kp]]), nodedof
                    )
                    force_value_vector += (
                        np.dot(N.transpose(), W) * abs(detJ) * wt[ip] * wt[jp] * wt[kp]
                    )
        # Mantém os valores (inclusive 0.0) apenas nos DOFs correspondentes a fc_type_dof
        force_value_vector = force_value_vector[fc_type_dof - 1 :: nodedof, 0]
        return force_value_vector, nodelist

    def __line_force_distribuition(
        Model,
        inci,
        coord,
        tabgeo,
        intgauss,
        node_list_fc,
        element_number,
        force_value,
        fc_type,
    ):
        nodedof = Model.modelinfo["nodedof"]
        shape = Model.modelinfo["shape"]
        edof = Model.modelinfo["elemdof"]
        nodelist = Model.shape.getNodeList(inci, element_number - 1)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        t = tabgeo[int(inci[element_number - 1, 3] - 1)]["THICKN"]
        test = np.isin(nodelist, node_list_fc, assume_unique=True)
        nodes = np.array(nodelist)[test]
        nodes_conec = np.where(test == True)[0]
        norm = np.zeros((2))
        if len(nodes_conec) < 2:
            nodes = np.repeat(nodes, 2)
            force_value_vector = np.zeros((len(nodes)))
            pass
        else:
            idx_conec = np.array2string(nodes_conec)
            get_side = Model.shape.getSideAxis(idx_conec[1:-1])
            if fc_type == "fx":
                T = np.array([[force_value], [0.0]])  # force_value
                norm[0] = 1
            elif fc_type == "fy":
                T = np.array([[0.0], [force_value]])  # force_value
                norm[1] = 1
            elif fc_type == "pressure":  #  -->[+]<--
                normal = Model.shape.getNormalEdge(elementcoord, get_side)
                PN = force_value * normal
                T = np.array(PN).reshape(-1, 1)
                norm = np.abs(np.array(normal))
                norm[np.abs(norm) < 1e-6] = 0
                norm = (norm >= 1e-6).astype(int)  # norm/norm
            else:
                T = np.array([[0.0], [0.0]])
            pt, wt = gauss_points(shape, intgauss)
            force_value_vector = np.zeros((edof, 1))
            for ip in range(intgauss):
                points = Model.shape.getIsoParaSide(get_side, pt[ip])
                N = Model.shape.getShapeFunctions(np.array(points), nodedof)
                J = Model.shape.getJacobian(np.array(points), elementcoord)
                detJ_e = Model.shape.getEdgeLength(J, get_side)
                force_value_vector += (
                    np.dot(np.array(N).transpose(), T) * t * abs(detJ_e) * wt[ip]
                )

            # Seleciona os DOFs ativos nos nós da borda sem eliminar valores nulos (0.0)
            active_dofs = np.where(norm == 1)[0]
            if active_dofs.size == 0:
                active_dofs = np.arange(nodedof)
            dof_indices = (
                nodes_conec[:, None] * nodedof + active_dofs[None, :]
            ).ravel()
            force_value_vector = force_value_vector[dof_indices, 0]
        return force_value_vector, nodes, norm

    def __surf_force_distribuition(
        Model,
        inci,
        coord,
        tabgeo,
        intgauss,
        node_list_fc,
        element_number,
        force_value,
        fc_type,
    ):
        shape_obj = Model.shape
        nodedof = Model.modelinfo["nodedof"]

        nodelist = shape_obj.getNodeList(inci, element_number - 1)
        elementcoord = shape_obj.getNodeCoord(coord, nodelist)

        nodelist_arr = np.asarray(nodelist)
        test = np.isin(nodelist_arr, node_list_fc, assume_unique=True)
        nodes = nodelist_arr[test]
        nodes_conec = np.flatnonzero(test)
        norm = np.zeros(3)

        if len(nodes_conec) < 3:
            nodes = np.repeat(nodes, 3)
            force_value_vector = np.zeros(len(nodes))
            return force_value_vector, nodes, norm

        idx_conec = np.array2string(nodes_conec)
        get_side = shape_obj.getSideAxis(idx_conec[1:-1])

        if fc_type == "fx":
            T = np.array([[force_value], [0.0], [0.0]])
            norm[0] = 1
        elif fc_type == "fy":
            T = np.array([[0.0], [force_value], [0.0]])
            norm[1] = 1
        elif fc_type == "fz":
            T = np.array([[0.0], [0.0], [force_value]])
            norm[2] = 1
        elif fc_type == "pressure":  #  -->[+]<--
            normal = shape_obj.getNormalFace(elementcoord, get_side)
            T = (force_value * normal).reshape(-1, 1)
            norm = (np.abs(np.asarray(normal)) >= 1e-6).astype(int)
        else:
            T = np.array([[0.0], [0.0], [0.0]])

        pt, wt = gauss_points(Model.modelinfo["shape"], intgauss)

        # A área da face não depende do ponto de Gauss: calculada uma vez
        detJ_a = abs(shape_obj.getAreaLength(get_side, elementcoord))

        get_iso = shape_obj.getIsoParaSide
        get_N = shape_obj.getShapeFunctions

        # T é constante: acumulo Σ N·w e faço Nᵀ·T uma única vez no final
        accN = 0.0
        for ip in range(intgauss):
            wi = wt[ip]
            for jp in range(intgauss):
                points = np.array(get_iso(get_side, [pt[ip], pt[jp]]))
                N = np.asarray(get_N(points, nodedof))
                accN = accN + N * (wi * wt[jp])

        force_value_vector = (accN.T @ T).ravel() * detJ_a

        # Seleciona os DOFs ativos nos nós da face sem eliminar valores nulos (0.0)
        active_dofs = np.where(norm == 1)[0]
        if active_dofs.size == 0:
            active_dofs = np.arange(nodedof)
        dof_indices = (
            nodes_conec[:, None] * nodedof + active_dofs[None, :]
        ).ravel()
        force_value_vector = force_value_vector[dof_indices]
        return force_value_vector, nodes, norm

    def __line_beam_distribuition(
        Model,
        inci,
        coord,
        tabgeo,
        intgauss,
        node_list_fc,
        element_number,
        force_value,
        fc_type,
    ):
        nodedof = Model.modelinfo["nodedof"]
        shape = Model.modelinfo["shape"]
        edof = Model.modelinfo["elemdof"]
        nodelist = Model.shape.getNodeList(inci, element_number - 1)
        elementcoord = Model.shape.getNodeCoord(coord, nodelist)
        test = np.isin(nodelist, node_list_fc, assume_unique=True)
        nodes = np.array(nodelist)[test]
        idx_conec = np.where(test == True)[0]
        norm = np.zeros((3))
        if len(idx_conec) < 2:
            nodes = np.repeat(nodes, 2)
            force_value_vector = np.zeros((len(nodes)))
            pass
        else:
            if fc_type == "fx":
                T = np.array([[force_value], [0.0], [0.0], [0.0]])  # force_value
                norm[0] = 1
            elif fc_type == "fy":
                T = np.array([[0.0], [force_value], [0.0], [0.0]])  # force_value
                norm[1] = 1
            elif fc_type == "fz":
                T = np.array([[0.0], [0.0], [force_value], [0.0]])  # force_value
                norm[2] = 1
            elif fc_type == "tx":
                T = np.array([[0.0], [0.0], [0.0], [force_value]])  # force_value
                norm[0] = 1
            elif fc_type == "pressure":  #  -->[+]<--
                noi = idx_conec[0]
                noj = idx_conec[1]
                dx = abs(elementcoord[noj, 0] - elementcoord[noi, 0])
                dy = abs(elementcoord[noj, 1] - elementcoord[noi, 1])
                L = np.sqrt(dx**2 + dy**2)
                norm[0] = dy / L
                norm[1] = dx / L
                T = -1 * (np.array([norm]).T) * force_value
            else:
                T = np.array([[0.0], [0.0], [0.0], [0.0]])
            idx_conec_str = np.array2string(idx_conec)
            get_side = Model.shape.getSideAxis(idx_conec_str[1:-1])
            pt, wt = gauss_points(shape, intgauss)
            force_value_vector = np.zeros((edof, 1))
            for ip in range(intgauss):
                points = Model.shape.getIsoParaSide(get_side, pt[ip])
                N = Model.shape.getShapeFunctions(np.array(points), nodedof)
                J = Model.shape.getJacobian(np.array(points), elementcoord)
                detJ_e = Model.shape.getEdgeLength(J, get_side)
                force_value_vector += (
                    np.dot(np.array(N).transpose(), T) * abs(detJ_e) * wt[ip]
                )

            # Seleciona os DOFs ativos na viga sem eliminar valores nulos (0.0)
            if fc_type == "pressure":
                active_dofs = np.where(norm[:2] == 1)[0]
                if active_dofs.size == 0:
                    active_dofs = np.array([0, 1])
            elif fc_type == "fy":
                active_dofs = np.array([1, 5])  # DOFs 2 e 6 (base 0)
            elif fc_type == "fz":
                active_dofs = np.array([2, 4])  # DOFs 3 e 5 (base 0)
            elif fc_type == "fx":
                active_dofs = np.array([0])     # DOF 1 (base 0)
            elif fc_type == "tx":
                active_dofs = np.array([3])     # DOF 4 (base 0)
            else:
                active_dofs = np.arange(nodedof)

            dof_indices = (
                idx_conec[:, None] * nodedof + active_dofs[None, :]
            ).ravel()
            force_value_vector = force_value_vector[dof_indices, 0]
        return force_value_vector, nodes, norm