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
        shape = Model.modelinfo["shape"]
        fdofs = Model.modelinfo["dofs"]["f"]

        # Decisão de tipo de elemento feita UMA vez, fora do loop
        if shape in ("hexa8", "tetr4"):
            elem_func = LoadStructural.__solid_strain_zero
            dof_types = np.array([fdofs["fx"], fdofs["fy"], fdofs["fz"]])
        elif shape in ("tria3", "tria6", "quad4", "quad8"):
            elem_func = LoadStructural.__plane_strain_zero
            dof_types = np.array([fdofs["fx"], fdofs["fy"]])
        else:
            raise ValueError(f"Elemento '{shape}' não suportado em ForceStrainZero")

        ndim = dof_types.size
        strain_zero = np.array(forcelist["VAL"])
        step = int(forcelist["STEP"])
        inci = Model.inci
        coord = Model.coord
        tabmat = Model.tabmat
        tabgeo = Model.tabgeo
        intgauss = Model.intgauss

        nodes_list = []
        vals_list = []
        for ee in range(inci.shape[0]):
            force_value_vector, nodelist = elem_func(
                Model, inci, coord, tabmat, tabgeo, intgauss, ee, strain_zero
            )
            nodes_list.append(np.repeat(nodelist, ndim))
            vals_list.append(np.asarray(force_value_vector, dtype=float).ravel())

        if not nodes_list:
            return np.empty((0, 4))

        # Montagem única, sem np.append dentro do loop
        nodes = np.concatenate(nodes_list).astype(np.int64)
        vals = np.concatenate(vals_list)
        n = nodes.size

        forcenodedof = np.empty((n, 4))
        forcenodedof[:, 0] = nodes
        forcenodedof[:, 1] = np.tile(dof_types, n // ndim)
        forcenodedof[:, 2] = vals
        forcenodedof[:, 3] = step
        return forcenodedof

    # def ForceStrainZero(Model, forcelist):
    #     forcenodedof = np.zeros((1, 4))
    #     strain_zero = np.array(forcelist["VAL"])
    #     inci = Model.inci
    #     coord = Model.coord
    #     tabmat = Model.tabmat
    #     tabgeo = Model.tabgeo
    #     intgauss = Model.intgauss
    #     for ee in range(inci.shape[0]):
    #         if (
    #             Model.modelinfo["shape"] == "hexa8"
    #             or Model.modelinfo["shape"] == "tetr4"
    #         ):
    #             force_value_vector, nodelist = LoadStructural.__solid_strain_zero(
    #                 Model, inci, coord, tabmat, tabgeo, intgauss, ee, strain_zero
    #             )
    #             fc_type_dof = np.tile(
    #                 [
    #                     Model.modelinfo["dofs"]["f"]["fx"],
    #                     Model.modelinfo["dofs"]["f"]["fy"],
    #                     Model.modelinfo["dofs"]["f"]["fz"],
    #                 ],
    #                 len(nodelist),
    #             )
    #             nodelist = np.repeat(nodelist, 3)

    #         elif (
    #             Model.modelinfo["shape"] == "tria3"
    #             or Model.modelinfo["shape"] == "tria6"
    #             or Model.modelinfo["shape"] == "quad4"
    #             or Model.modelinfo["shape"] == "quad8"
    #         ):
    #             force_value_vector, nodelist = LoadStructural.__plane_strain_zero(
    #                 Model, inci, coord, tabmat, tabgeo, intgauss, ee, strain_zero
    #             )
    #             fc_type_dof = np.tile(
    #                 [
    #                     Model.modelinfo["dofs"]["f"]["fx"],
    #                     Model.modelinfo["dofs"]["f"]["fy"],
    #                 ],
    #                 len(nodelist),
    #             )
    #             nodelist = np.repeat(nodelist, 2)

    #         else:
    #             force_value_vector, nodelist = [0], [0]

    #         for j in range(len(nodelist)):
    #             fcdof = np.array(
    #                 [
    #                     [
    #                         int(nodelist[j]),
    #                         fc_type_dof[j],
    #                         force_value_vector[j],
    #                         int(forcelist["STEP"]),
    #                     ]
    #                 ]
    #             )
    #             forcenodedof = np.append(forcenodedof, fcdof, axis=0)
    #     forcenodedof = forcenodedof[1::][::]
    #     return forcenodedof

    def ForceNodeLoadApply(Model, forcelist):
        nodelist = [
            forcelist["DIR"],
            forcelist["LOCX"],
            forcelist["LOCY"],
            forcelist["LOCZ"],
            forcelist["TAG"],
            forcelist["MESHNODE"],
        ]
        node_list_fc, dir_fc = get_nodes_from_list(nodelist, Model.coord, Model.regions)

        nodes = np.asarray(node_list_fc).ravel()
        n = nodes.size

        # Montagem direta, sem loop e sem np.append
        forcenodedof = np.empty((n, 4))
        forcenodedof[:, 0] = nodes.astype(np.int64)
        forcenodedof[:, 1] = Model.modelinfo["dofs"]["f"][forcelist["DOF"]]
        forcenodedof[:, 2] = float(forcelist["VAL"])
        forcenodedof[:, 3] = int(forcelist["STEP"])
        return forcenodedof

    # def ForceNodeLoadApply(Model, forcelist):
    #     forcenodedof = np.zeros((1, 4))
    #     nodelist = [
    #         forcelist["DIR"],
    #         forcelist["LOCX"],
    #         forcelist["LOCY"],
    #         forcelist["LOCZ"],
    #         forcelist["TAG"],
    #         forcelist["MESHNODE"],
    #     ]
    #     node_list_fc, dir_fc = get_nodes_from_list(nodelist, Model.coord, Model.regions)
    #     force_value_vector = np.ones_like(node_list_fc) * float(forcelist["VAL"])
    #     fc_type_dof = Model.modelinfo["dofs"]["f"][forcelist["DOF"]] * np.ones_like(
    #         node_list_fc
    #     )
    #     for j in range(len(node_list_fc)):
    #         fcdof = np.array(
    #             [
    #                 [
    #                     int(node_list_fc[j]),
    #                     fc_type_dof[j],
    #                     force_value_vector[j],
    #                     int(forcelist["STEP"]),
    #                 ]
    #             ]
    #         )
    #         forcenodedof = np.append(forcenodedof, fcdof, axis=0)
    #     forcenodedof = forcenodedof[1::][::]
    #     return forcenodedof

    def BodyForce(Model, forcelist):
        shape = Model.modelinfo["shape"]

        # Decisão do tipo de elemento feita UMA vez, fora do loop
        if shape in ("hexa8", "tetr4"):
            elem_func = LoadStructural.__solid_body_force_volumetric
        elif shape in ("tria3", "tria6", "quad4", "quad8"):
            elem_func = LoadStructural.__plane_body_force_volumetric
        else:
            elem_func = None

        gravity_value = float(forcelist["VAL"])
        step = int(forcelist["STEP"])
        fc_type_dof = Model.modelinfo["dofs"]["f"][forcelist["DOF"]]
        inci = Model.inci
        coord = Model.coord
        tabmat = Model.tabmat
        tabgeo = Model.tabgeo
        intgauss = Model.intgauss
        nelem = inci.shape[0]

        if nelem == 0:
            return np.empty((0, 4))

        # Forma não suportada: o original gerava uma linha [0, dof, 0, step] por elemento
        if elem_func is None:
            forcenodedof = np.zeros((nelem, 4))
            forcenodedof[:, 1] = fc_type_dof
            forcenodedof[:, 3] = step
            return forcenodedof

        nodes_list = []
        vals_list = []
        for ee in range(nelem):
            force_value_vector, nodelist = elem_func(
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
            nodes_list.append(np.asarray(nodelist))
            vals_list.append(np.asarray(force_value_vector, dtype=float).ravel())

        # Montagem única, sem np.append dentro do loop
        nodes = np.concatenate(nodes_list).astype(np.int64)
        vals = np.concatenate(vals_list)

        forcenodedof = np.empty((nodes.size, 4))
        forcenodedof[:, 0] = nodes
        forcenodedof[:, 1] = fc_type_dof
        forcenodedof[:, 2] = vals
        forcenodedof[:, 3] = step
        return forcenodedof

    # def BodyForce(Model, forcelist):
    #     forcenodedof = np.zeros((1, 4))
    #     gravity_value = float(forcelist["VAL"])
    #     inci = Model.inci
    #     coord = Model.coord
    #     tabmat = Model.tabmat
    #     tabgeo = Model.tabgeo
    #     intgauss = Model.intgauss
    #     fc_type_dof = Model.modelinfo["dofs"]["f"][forcelist["DOF"]]
    #     for ee in range(inci.shape[0]):

    #         if (
    #             Model.modelinfo["shape"] == "hexa8"
    #             or Model.modelinfo["shape"] == "tetr4"
    #         ):
    #             force_value_vector, nodelist = (
    #                 LoadStructural.__solid_body_force_volumetric(
    #                     Model,
    #                     inci,
    #                     coord,
    #                     tabmat,
    #                     tabgeo,
    #                     intgauss,
    #                     ee,
    #                     gravity_value,
    #                     fc_type_dof,
    #                 )
    #             )
    #         elif (
    #             Model.modelinfo["shape"] == "tria3"
    #             or Model.modelinfo["shape"] == "tria6"
    #             or Model.modelinfo["shape"] == "quad4"
    #             or Model.modelinfo["shape"] == "quad8"
    #         ):
    #             force_value_vector, nodelist = (
    #                 LoadStructural.__plane_body_force_volumetric(
    #                     Model,
    #                     inci,
    #                     coord,
    #                     tabmat,
    #                     tabgeo,
    #                     intgauss,
    #                     ee,
    #                     gravity_value,
    #                     fc_type_dof,
    #                 )
    #             )
    #         else:
    #             force_value_vector, nodelist = [0], [0]

    #         for j in range(len(nodelist)):
    #             fcdof = np.array(
    #                 [
    #                     [
    #                         int(nodelist[j]),
    #                         fc_type_dof,
    #                         force_value_vector[j],
    #                         int(forcelist["STEP"]),
    #                     ]
    #                 ]
    #             )
    #             forcenodedof = np.append(forcenodedof, fcdof, axis=0)
    #     forcenodedof = forcenodedof[1::][::]
    #     return forcenodedof


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
        inci = Model.inci
        coord = Model.coord
        tabgeo = Model.tabgeo
        intgauss = Model.intgauss
        fc_type = forcelist["DOF"]
        step = int(forcelist["STEP"])
        elmlist = get_elemen_from_nodelist(inci, node_list_fc)

        if len(elmlist) == 0:
            return np.empty((0, 4))

        # Tudo que é constante entre elementos é calculado UMA vez, fora do loop
        dofs_f = Model.modelinfo["dofs"]["f"]
        nodedof = len(Model.element.getElementSet()["dofs"]["d"])
        is_pressure = fc_type == "pressure"
        if is_pressure:
            nodesce = Model.shape.getShapeSet()["nodesconecedge"]
            fx = dofs_f["fx"]
            fy = dofs_f["fy"]
            fxfy_tile = np.tile([fx, fy], nodesce)
        else:
            dof_const = dofs_f[fc_type]

        line_force = LoadStructural.__line_force_distribuition

        nodes_list = []
        dofs_list = []
        vals_list = []
        for elm in elmlist:
            force_value_vector, nodes, norm = line_force(
                Model,
                inci,
                coord,
                tabgeo,
                intgauss,
                node_list_fc,
                elm,
                force_value,
                fc_type,
            )

            if len(force_value_vector) > len(nodes):
                nodes = np.repeat(nodes, nodedof)
            else:
                nodes = np.asarray(nodes)
            n = nodes.size

            if is_pressure:
                n0 = int(norm[0])
                n1 = int(norm[1])
                if n0 == 1 and n1 == 0:
                    dofs = np.full(n, fx)
                elif n0 == 0 and n1 == 1:
                    dofs = np.full(n, fy)
                else:
                    dofs = fxfy_tile
            else:
                dofs = np.full(n, dof_const)

            # n linhas por elemento, como no loop original (range(len(nodelist)))
            nodes_list.append(nodes)
            dofs_list.append(dofs[:n])
            vals_list.append(np.asarray(force_value_vector, dtype=float).ravel()[:n])

        # Montagem única, sem np.append dentro do loop
        nodes_all = np.concatenate(nodes_list).astype(np.int64)

        forcenodedof = np.empty((nodes_all.size, 4))
        forcenodedof[:, 0] = nodes_all
        forcenodedof[:, 1] = np.concatenate(dofs_list)
        forcenodedof[:, 2] = np.concatenate(vals_list)
        forcenodedof[:, 3] = step
        return forcenodedof

    # def ForceEdgeLoadApply(Model, forcelist):
    #     nodelist = [
    #         forcelist["DIR"],
    #         forcelist["LOCX"],
    #         forcelist["LOCY"],
    #         forcelist["LOCZ"],
    #         forcelist["TAG"],
    #         forcelist["MESHNODE"],
    #     ]  # forcelist[3:]
    #     node_list_fc, dir_fc = get_nodes_from_list(nodelist, Model.coord, Model.regions)
    #     force_value = float(forcelist["VAL"])
    #     force_dirc = forcelist["DOF"]
    #     inci = Model.inci
    #     coord = Model.coord
    #     tabmat = Model.tabmat
    #     tabgeo = Model.tabgeo
    #     intgauss = Model.intgauss
    #     fc_type = forcelist["DOF"]
    #     elmlist = get_elemen_from_nodelist(inci, node_list_fc)
    #     forcenodedof = np.zeros((1, 4))
    #     for ee in range(len(elmlist)):
    #         force_value_vector, nodelist, norm = (
    #             LoadStructural.__line_force_distribuition(
    #                 Model,
    #                 inci,
    #                 coord,
    #                 tabgeo,
    #                 intgauss,
    #                 node_list_fc,
    #                 elmlist[ee],
    #                 force_value,
    #                 fc_type,
    #             )
    #         )

    #         elem_set = Model.element.getElementSet()
    #         nodedof = len(elem_set["dofs"]["d"])

    #         if len(force_value_vector) > len(nodelist):
    #             nodelist = np.repeat(nodelist, nodedof)
    #             # fc_type_dof = np.tile([modelinfo["dofs"]["f"]["fx"], modelinfo["dofs"]["f"]["fy"]], int(len(nodelist)/nodedof))

    #         if forcelist["DOF"] == "pressure":
    #             shape_set = Model.shape.getShapeSet()
    #             nodesce = shape_set["nodesconecedge"]
    #             if int(norm[0]) == 1 and int(norm[1]) == 0:
    #                 fc_type_dof = Model.modelinfo["dofs"]["f"]["fx"] * np.ones_like(
    #                     nodelist
    #                 )
    #             elif int(norm[0]) == 0 and int(norm[1]) == 1:
    #                 fc_type_dof = Model.modelinfo["dofs"]["f"]["fy"] * np.ones_like(
    #                     nodelist
    #                 )
    #             else:
    #                 fc_type_dof = np.tile(
    #                     [
    #                         Model.modelinfo["dofs"]["f"]["fx"],
    #                         Model.modelinfo["dofs"]["f"]["fy"],
    #                     ],
    #                     nodesce,
    #                 )
    #         else:
    #             fc_type_dof = Model.modelinfo["dofs"]["f"][
    #                 forcelist["DOF"]
    #             ] * np.ones_like(nodelist)

    #         for j in range(len(nodelist)):
    #             fcdof = np.array(
    #                 [
    #                     [
    #                         int(nodelist[j]),
    #                         fc_type_dof[j],
    #                         force_value_vector[j],
    #                         int(forcelist["STEP"]),
    #                     ]
    #                 ]
    #             )
    #             forcenodedof = np.append(forcenodedof, fcdof, axis=0)
    #     forcenodedof = forcenodedof[1::][::]
    #     return forcenodedof

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
        inci = Model.inci
        coord = Model.coord
        tabgeo = Model.tabgeo
        intgauss = Model.intgauss
        fc_type = forcelist["DOF"]
        step = int(forcelist["STEP"])
        elmlist = get_elemen_from_nodelist(inci, node_list_fc)

        if len(elmlist) == 0:
            return np.empty((0, 4))

        # Tudo que é constante entre elementos é calculado UMA vez, fora do loop
        dofs_f = Model.modelinfo["dofs"]["f"]
        nodedof = len(Model.element.getElementSet()["dofs"]["d"])
        is_pressure = fc_type == "pressure"

        if is_pressure:
            nodescf = Model.shape.getShapeSet()["nodesconecface"]
            fx = dofs_f["fx"]
            fy = dofs_f["fy"]
            fz = dofs_f["fz"]
            # direção única da normal -> dof constante
            single_dof = {(1, 0, 0): fx, (0, 1, 0): fy, (0, 0, 1): fz}
            # duas componentes -> padrão alternado pré-montado
            multi_dof = {
                (1, 1, 0): np.tile([fx, fy], nodescf),
                (1, 0, 1): np.tile([fx, fz], nodescf),
                (0, 1, 1): np.tile([fy, fz], nodescf),
            }
            default_dof = np.tile([fx, fy, fz], nodescf)
        else:
            dof_const = dofs_f[fc_type]

        surf_force = LoadStructural.__surf_force_distribuition

        nodes_list = []
        dofs_list = []
        vals_list = []
        for elm in elmlist:
            force_value_vector, nodes, norm = surf_force(
                Model,
                inci,
                coord,
                tabgeo,
                intgauss,
                node_list_fc,
                elm,
                force_value,
                fc_type,
            )

            nvec = len(force_value_vector)
            if nvec > len(nodes):
                nodes = np.repeat(nodes, int(nvec / nodedof))
            else:
                nodes = np.asarray(nodes)
            n = nodes.size

            if is_pressure:
                key = (int(norm[0]), int(norm[1]), int(norm[2]))
                if key in single_dof:
                    dofs = np.full(n, single_dof[key])
                else:
                    dofs = multi_dof.get(key, default_dof)
            else:
                dofs = np.full(n, dof_const)

            # n linhas por elemento, como no loop original (range(len(nodelist)))
            nodes_list.append(nodes)
            dofs_list.append(dofs[:n])
            vals_list.append(np.asarray(force_value_vector, dtype=float).ravel()[:n])

        # Montagem única, sem np.append dentro do loop
        nodes_all = np.concatenate(nodes_list).astype(np.int64)

        forcenodedof = np.empty((nodes_all.size, 4))
        forcenodedof[:, 0] = nodes_all
        forcenodedof[:, 1] = np.concatenate(dofs_list)
        forcenodedof[:, 2] = np.concatenate(vals_list)
        forcenodedof[:, 3] = step
        return forcenodedof

    # def ForceSurfLoadApply(Model, forcelist):
    #     nodelist = [
    #         forcelist["DIR"],
    #         forcelist["LOCX"],
    #         forcelist["LOCY"],
    #         forcelist["LOCZ"],
    #         forcelist["TAG"],
    #         forcelist["MESHNODE"],
    #     ]  # forcelist[3:]
    #     node_list_fc, dir_fc = get_nodes_from_list(nodelist, Model.coord, Model.regions)
    #     force_value = float(forcelist["VAL"])
    #     force_dirc = forcelist["DOF"]
    #     inci = Model.inci
    #     coord = Model.coord
    #     tabmat = Model.tabmat
    #     tabgeo = Model.tabgeo
    #     intgauss = Model.intgauss
    #     fc_type = forcelist["DOF"]
    #     elmlist = get_elemen_from_nodelist(inci, node_list_fc)
    #     forcenodedof = np.zeros((1, 4))
    #     for ee in range(len(elmlist)):
    #         force_value_vector, nodelist, norm = (
    #             LoadStructural.__surf_force_distribuition(
    #                 Model,
    #                 inci,
    #                 coord,
    #                 tabgeo,
    #                 intgauss,
    #                 node_list_fc,
    #                 elmlist[ee],
    #                 force_value,
    #                 fc_type,
    #             )
    #         )

    #         elem_set = Model.element.getElementSet()
    #         nodedof = len(elem_set["dofs"]["d"])

    #         if len(force_value_vector) > len(nodelist):
    #             nodelist = np.repeat(nodelist, int(len(force_value_vector) / nodedof))
    #             # fc_type_dof = np.tile([modelinfo["dofs"]["f"]["fx"], modelinfo["dofs"]["f"]["fy"]], int(len(nodelist)/nodedof))

    #         if forcelist["DOF"] == "pressure":
    #             shape_set = Model.shape.getShapeSet()
    #             nodescf = shape_set["nodesconecface"]
    #             if int(norm[0]) == 1 and int(norm[1]) == 0 and int(norm[2]) == 0:
    #                 fc_type_dof = Model.modelinfo["dofs"]["f"]["fx"] * np.ones_like(
    #                     nodelist
    #                 )
    #             elif int(norm[0]) == 0 and int(norm[1]) == 1 and int(norm[2]) == 0:
    #                 fc_type_dof = Model.modelinfo["dofs"]["f"]["fy"] * np.ones_like(
    #                     nodelist
    #                 )
    #             elif int(norm[0]) == 0 and int(norm[1]) == 0 and int(norm[2]) == 1:
    #                 fc_type_dof = Model.modelinfo["dofs"]["f"]["fz"] * np.ones_like(
    #                     nodelist
    #                 )
    #             elif int(norm[0]) == 1 and int(norm[1]) == 1 and int(norm[2]) == 0:
    #                 fc_type_dof = np.tile(
    #                     [
    #                         Model.modelinfo["dofs"]["f"]["fx"],
    #                         Model.modelinfo["dofs"]["f"]["fy"],
    #                     ],
    #                     nodescf,
    #                 )
    #             elif int(norm[0]) == 1 and int(norm[1]) == 0 and int(norm[2]) == 1:
    #                 fc_type_dof = np.tile(
    #                     [
    #                         Model.modelinfo["dofs"]["f"]["fx"],
    #                         Model.modelinfo["dofs"]["f"]["fz"],
    #                     ],
    #                     nodescf,
    #                 )
    #             elif int(norm[0]) == 0 and int(norm[1]) == 1 and int(norm[2]) == 1:
    #                 fc_type_dof = np.tile(
    #                     [
    #                         Model.modelinfo["dofs"]["f"]["fy"],
    #                         Model.modelinfo["dofs"]["f"]["fz"],
    #                     ],
    #                     nodescf,
    #                 )
    #             else:
    #                 fc_type_dof = np.tile(
    #                     [
    #                         Model.modelinfo["dofs"]["f"]["fx"],
    #                         Model.modelinfo["dofs"]["f"]["fy"],
    #                         Model.modelinfo["dofs"]["f"]["fz"],
    #                     ],
    #                     nodescf,
    #                 )
    #         else:
    #             fc_type_dof = Model.modelinfo["dofs"]["f"][
    #                 forcelist["DOF"]
    #             ] * np.ones_like(nodelist)

    #         for j in range(len(nodelist)):
    #             fcdof = np.array(
    #                 [
    #                     [
    #                         int(nodelist[j]),
    #                         fc_type_dof[j],
    #                         force_value_vector[j],
    #                         int(forcelist["STEP"]),
    #                     ]
    #                 ]
    #             )
    #             forcenodedof = np.append(forcenodedof, fcdof, axis=0)
    #     forcenodedof = forcenodedof[1::][::]
    #     return forcenodedof

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
        inci = Model.inci
        coord = Model.coord
        tabgeo = Model.tabgeo
        intgauss = Model.intgauss
        fc_type = forcelist["DOF"]
        step = int(forcelist["STEP"])
        elmlist = get_elemen_from_nodelist(inci, node_list_fc)

        if len(elmlist) == 0:
            return np.empty((0, 4))

        # Tudo que é constante entre elementos é calculado UMA vez, fora do loop
        dofs_f = Model.modelinfo["dofs"]["f"]
        nodedof = Model.modelinfo["nodedof"]
        is_pressure = fc_type == "pressure"

        if is_pressure:
            fx = dofs_f["fx"]
            fy = dofs_f["fy"]
        else:
            patterns = {"fy": [2, 6], "fz": [3, 5], "fx": [1], "tx": [4]}
            if fc_type not in patterns:
                # no original: "pass" -> fc_type_dof indefinido (NameError)
                raise ValueError(f"DOF '{fc_type}' não suportado em ForceBeamLoadApply")
            pattern = patterns[fc_type]

        dof_cache = {}  # padrões tile reaproveitados entre elementos de mesmo tamanho
        beam_force = LoadStructural.__line_beam_distribuition

        nodes_list = []
        dofs_list = []
        vals_list = []
        for elm in elmlist:
            force_value_vector, nodeslist, norm = beam_force(
                Model,
                inci,
                coord,
                tabgeo,
                intgauss,
                node_list_fc,
                elm,
                force_value,
                fc_type,
            )

            nodes_orig = np.asarray(nodeslist)
            m = nodes_orig.size
            if len(force_value_vector) > m:
                nodes = np.repeat(nodes_orig, 2)
            else:
                nodes = nodes_orig
            n = nodes.size

            if is_pressure:
                n0 = int(norm[0])
                n1 = int(norm[1])
                if n0 == 1 and n1 == 0:
                    dofs = np.full(n, fx)
                elif n0 == 0 and n1 == 1:
                    dofs = np.full(n, fy)
                else:
                    key = int(m / nodedof)
                    dofs = dof_cache.get(key)
                    if dofs is None:
                        dofs = dof_cache[key] = np.tile([fx, fy], key)
            else:
                dofs = dof_cache.get(m)
                if dofs is None:
                    dofs = dof_cache[m] = np.tile(pattern, m)

            # n linhas por elemento, como no loop original (range(len(nodeslist)))
            nodes_list.append(nodes)
            dofs_list.append(dofs[:n])
            vals_list.append(np.asarray(force_value_vector, dtype=float).ravel()[:n])

        # Montagem única, sem np.append dentro do loop
        nodes_all = np.concatenate(nodes_list).astype(np.int64)

        forcenodedof = np.empty((nodes_all.size, 4))
        forcenodedof[:, 0] = nodes_all
        forcenodedof[:, 1] = np.concatenate(dofs_list)
        forcenodedof[:, 2] = np.concatenate(vals_list)
        forcenodedof[:, 3] = step
        return forcenodedof

    # def ForceBeamLoadApply(Model, forcelist):
    #     nodelist = [
    #         forcelist["DIR"],
    #         forcelist["LOCX"],
    #         forcelist["LOCY"],
    #         forcelist["LOCZ"],
    #         forcelist["TAG"],
    #         forcelist["MESHNODE"],
    #     ]  # forcelist[3:]
    #     node_list_fc, dir_fc = get_nodes_from_list(nodelist, Model.coord, Model.regions)
    #     force_value = float(forcelist["VAL"])
    #     force_dirc = forcelist["DOF"]
    #     inci = Model.inci
    #     coord = Model.coord
    #     tabmat = Model.tabmat
    #     tabgeo = Model.tabgeo
    #     intgauss = Model.intgauss
    #     fc_type = forcelist["DOF"]
    #     elmlist = get_elemen_from_nodelist(inci, node_list_fc)
    #     forcenodedof = np.zeros((1, 4))
    #     for ee in range(len(elmlist)):
    #         force_value_vector, nodeslist, norm = (
    #             LoadStructural.__line_beam_distribuition(
    #                 Model,
    #                 inci,
    #                 coord,
    #                 tabgeo,
    #                 intgauss,
    #                 node_list_fc,
    #                 elmlist[ee],
    #                 force_value,
    #                 fc_type,
    #             )
    #         )

    #         nodedof = Model.modelinfo["nodedof"]

    #         nodes = nodeslist
    #         if len(force_value_vector) > len(nodeslist):
    #             nodeslist = np.repeat(nodeslist, 2)
    #             # fc_type_dof = np.tile([modelinfo["dofs"]["f"]["fx"], modelinfo["dofs"]["f"]["fy"]], int(len(nodelist)/nodedof))

    #         if forcelist["DOF"] == "pressure":
    #             if int(norm[0]) == 1 and int(norm[1]) == 0:
    #                 fc_type_dof = Model.modelinfo["dofs"]["f"]["fx"] * np.ones_like(
    #                     nodeslist
    #                 )
    #             elif int(norm[0]) == 0 and int(norm[1]) == 1:
    #                 fc_type_dof = Model.modelinfo["dofs"]["f"]["fy"] * np.ones_like(
    #                     nodeslist
    #                 )
    #             else:
    #                 fc_type_dof = np.tile(
    #                     [
    #                         Model.modelinfo["dofs"]["f"]["fx"],
    #                         Model.modelinfo["dofs"]["f"]["fy"],
    #                     ],
    #                     int(len(nodes) / nodedof),
    #                 )  # modelinfo["dofs"]["f"]["fy"] * np.ones_like(nodelist)
    #         else:
    #             if force_dirc == "fy":
    #                 fc_type_dof = np.tile([2, 6], len(nodes))
    #             elif force_dirc == "fz":
    #                 fc_type_dof = np.tile([3, 5], len(nodes))
    #             elif force_dirc == "fx":
    #                 fc_type_dof = np.tile([1], len(nodes))
    #             elif force_dirc == "tx":
    #                 fc_type_dof = np.tile([4], len(nodes))
    #             else:
    #                 pass

    #         for j in range(len(nodeslist)):
    #             fcdof = np.array(
    #                 [
    #                     [
    #                         int(nodeslist[j]),
    #                         fc_type_dof[j],
    #                         force_value_vector[j],
    #                         int(forcelist["STEP"]),
    #                     ]
    #                 ]
    #             )
    #             forcenodedof = np.append(forcenodedof, fcdof, axis=0)
    #     forcenodedof = forcenodedof[1::][::]
    #     return forcenodedof

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
        shape_obj = Model.shape
        H = Model.element.getElementSet()["H"]
        nodedof = Model.modelinfo["nodedof"]
        edof = Model.modelinfo["elemdof"]

        nodelist = shape_obj.getNodeList(inci, element_number)
        elementcoord = shape_obj.getNodeCoord(coord, nodelist)
        C = Model.material.getElasticTensor(tabmat, inci, element_number)
        pt, wt = gauss_points(Model.modelinfo["shape"], intgauss)

        # B^T C e0 = B^T (C e0): o produto C @ e0 é constante no elemento
        # e vira matriz x vetor (bem mais barato que matriz x matriz)
        Ce0 = C @ np.asarray(strain_zero).ravel()

        # métodos ligados uma vez, fora do loop
        get_detJ = shape_obj.getdetJacobi
        get_diffN = shape_obj.getDiffShapeFuntion
        get_invJ = shape_obj.getinvJacobi
        get_B = shape_obj.getB

        force_value_vector = np.zeros(edof)
        for ip in range(intgauss):
            wi = wt[ip]
            for jp in range(intgauss):
                # ponto de Gauss criado uma única vez (antes eram 3 vezes)
                xi = np.array([pt[ip], pt[jp]])
                detJ = get_detJ(xi, elementcoord)
                diffN = get_diffN(xi, nodedof)
                invJ = get_invJ(xi, elementcoord, nodedof)
                B = get_B(H, invJ, diffN)
                force_value_vector += (B.T @ Ce0) * (abs(detJ) * wi * wt[jp])

        return force_value_vector, nodelist

    # def __plane_strain_zero(
    #     Model,
    #     inci,
    #     coord,
    #     tabmat,
    #     tabgeo,
    #     intgauss,
    #     element_number,
    #     strain_zero,
    # ):
    #     # internal balance forces
    #     elem_set = Model.element.getElementSet()
    #     H = elem_set["H"]
    #     nodedof = Model.modelinfo["nodedof"]
    #     shape = Model.modelinfo["shape"]
    #     edof = Model.modelinfo["elemdof"]
    #     nodelist = Model.shape.getNodeList(inci, element_number)
    #     elementcoord = Model.shape.getNodeCoord(coord, nodelist)
    #     C = Model.material.getElasticTensor(tabmat, inci, element_number)
    #     pt, wt = gauss_points(shape, intgauss)
    #     len_sigma = len(elem_set["tensor"])
    #     force_value_vector = np.zeros((edof, 1))
    #     for ip in range(intgauss):
    #         for jp in range(intgauss):
    #             detJ = Model.shape.getdetJacobi(
    #                 np.array([pt[ip], pt[jp]]), elementcoord
    #             )
    #             diffN = Model.shape.getDiffShapeFuntion(
    #                 np.array([pt[ip], pt[jp]]), nodedof
    #             )
    #             invJ = Model.shape.getinvJacobi(
    #                 np.array([pt[ip], pt[jp]]), elementcoord, nodedof
    #             )
    #             B = B = Model.shape.getB(H, invJ, diffN)
    #             force_value_vector += (
    #                 np.dot(np.dot(B.transpose(), C), strain_zero.reshape(-1, 1))
    #                 * abs(detJ)
    #                 * wt[ip]
    #                 * wt[jp]
    #             )
    #     force_value_vector = np.reshape(force_value_vector, (edof))
    #     return force_value_vector, nodelist

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
        shape_obj = Model.shape
        nodedof = Model.modelinfo["nodedof"]
        edof = Model.modelinfo["elemdof"]

        nodelist = shape_obj.getNodeList(inci, element_number)
        elementcoord = shape_obj.getNodeCoord(coord, nodelist)
        C = Model.material.getElasticTensor(tabmat, inci, element_number)
        pt, wt = gauss_points(Model.modelinfo["shape"], intgauss)

        # Bᵀ(C ε₀): produto constante no elemento, calculado uma vez
        Ce0 = C @ np.asarray(strain_zero).ravel()

        # O integrando não depende de kp (o ponto só usa pt[ip], pt[jp]),
        # então o laço em kp só multiplica pela soma dos pesos.
        wk_sum = np.sum(np.asarray(wt)[:intgauss])

        # métodos ligados uma vez, fora do loop
        get_detJ = shape_obj.getdetJacobi
        get_diffN = shape_obj.getDiffShapeFuntion
        get_invJ = shape_obj.getinvJacobi
        get_B = Model.element.getB

        force_value_vector = np.zeros(edof)
        for ip in range(intgauss):
            for jp in range(intgauss):
                wij = wt[ip] * wt[jp]
                for kp in range(intgauss):
                    xi = np.array([pt[ip], pt[jp], pt[kp]])
                    detJ = get_detJ(xi, elementcoord)
                    diffN = get_diffN(xi, nodedof)
                    invJ = get_invJ(xi, elementcoord, nodedof)
                    B = get_B(diffN, invJ)
                    force_value_vector += (B.T @ Ce0) * (abs(detJ) * wij * wt[kp])

        return force_value_vector, nodelist

    # def __solid_strain_zero(
    #     Model,
    #     inci,
    #     coord,
    #     tabmat,
    #     tabgeo,
    #     intgauss,
    #     element_number,
    #     strain_zero,
    # ):
    #     # internal balance forces
    #     elem_set = Model.element.getElementSet()
    #     nodedof = Model.modelinfo["nodedof"]
    #     shape = Model.modelinfo["shape"]
    #     edof = Model.modelinfo["elemdof"]
    #     nodelist = Model.shape.getNodeList(inci, element_number)
    #     elementcoord = Model.shape.getNodeCoord(coord, nodelist)
    #     C = Model.material.getElasticTensor(tabmat, inci, element_number)
    #     pt, wt = gauss_points(shape, intgauss)
    #     len_sigma = len(elem_set["tensor"])
    #     force_value_vector = np.zeros((edof, 1))
    #     for ip in range(intgauss):
    #         for jp in range(intgauss):
    #             for kp in range(intgauss):
    #                 detJ = Model.shape.getdetJacobi(
    #                     np.array([pt[ip], pt[jp]]), elementcoord
    #                 )
    #                 diffN = Model.shape.getDiffShapeFuntion(
    #                     np.array([pt[ip], pt[jp]]), nodedof
    #                 )
    #                 invJ = Model.shape.getinvJacobi(
    #                     np.array([pt[ip], pt[jp]]), elementcoord, nodedof
    #                 )
    #                 B = Model.element.getB(diffN, invJ)
    #                 force_value_vector += (
    #                     np.dot(np.dot(B.transpose(), C), strain_zero.reshape(-1, 1))
    #                     * abs(detJ)
    #                     * wt[ip]
    #                     * wt[jp]
    #                     * wt[kp]
    #                 )
    #     force_value_vector = np.reshape(force_value_vector, (edof))
    #     return force_value_vector, nodelist

    # def __line_body_force_volumetric(
    #     Model,
    #     inci,
    #     coord,
    #     tabmat,
    #     tabgeo,
    #     intgauss,
    #     element_number,
    #     gravity_value,
    #     fc_type_dof,
    # ):
    #     # body force line
    #     elem_set = Model.element.getElementSet()
    #     nodedof = len(elem_set["dofs"]["d"])
    #     shape_set = Model.shape.getShapeSet()
    #     nodecon = len(shape_set["nodes"])
    #     shape = shape_set["key"]
    #     edof = nodecon * nodedof
    #     nodelist = Model.shape.getNodeList(inci, element_number)
    #     elementcoord = Model.shape.getNodeCoord(coord, nodelist)
    #     R = tabmat[int(inci[element_number, 2]) - 1]["RHO"]
    #     AREA = tabgeo[int(inci[element_number, 3] - 1)]["AREACS"]
    #     pt, wt = gauss_points(shape, intgauss)
    #     G = gravity_value
    #     W = np.zeros((nodedof, 1))
    #     W[fc_type_dof - 1, 0] = R * G
    #     force_value_vector = np.zeros((edof, 1))
    #     idx_conec = np.array2string(idx_conec)
    #     get_side = Model.shape.getSideAxis(idx_conec[1:-1])
    #     for ip in range(intgauss):
    #         points = Model.shape.getIsoParaSide(get_side, pt[ip])
    #         N = Model.shape.getShapeFunctions(np.array(points), nodedof)
    #         J = Model.shape.getJacobian(np.array(points), elementcoord)
    #         detJ_e = Model.shape.getEdgeLength(J, get_side)
    #         force_value_vector += np.dot(np.array(N).transpose(), W) * AREA * abs(detJ_e) * wt[ip]
    #     force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
    #     return force_value_vector, nodelist

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
        shape_obj = Model.shape
        nodedof = Model.modelinfo["nodedof"]
        edof = Model.modelinfo["elemdof"]

        nodelist = shape_obj.getNodeList(inci, element_number)
        elementcoord = shape_obj.getNodeCoord(coord, nodelist)
        R = tabmat[int(inci[element_number, 2]) - 1]["RHO"]
        t = tabgeo[int(inci[element_number, 3] - 1)]["THICKN"]
        pt, wt = gauss_points(Model.modelinfo["shape"], intgauss)

        # W tem uma única componente não nula (R*G em fc_type_dof-1), então
        # Nᵀ·W = R*G * N[fc_type_dof-1, :]. Evita o produto matriz×vetor
        # e a criação de W a cada ponto de Gauss.
        idx = fc_type_dof - 1

        get_detJ = shape_obj.getdetJacobi
        get_N = shape_obj.getShapeFunctions

        acc = np.zeros(edof)
        for ip in range(intgauss):
            wi = wt[ip]
            for jp in range(intgauss):
                xi = np.array([pt[ip], pt[jp]])  # criado uma vez por ponto
                detJ = get_detJ(xi, elementcoord)
                N = get_N(xi, nodedof)
                acc += N[idx, :] * (abs(detJ) * wi * wt[jp])

        # constantes aplicadas uma única vez, no final
        force_value_vector = acc * (R * gravity_value * t)
        force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
        return force_value_vector, nodelist

    # def __plane_body_force_volumetric(
    #     Model,
    #     inci,
    #     coord,
    #     tabmat,
    #     tabgeo,
    #     intgauss,
    #     element_number,
    #     gravity_value,
    #     fc_type_dof,
    # ):
    #     # body force plane
    #     nodedof = Model.modelinfo["nodedof"]
    #     shape = Model.modelinfo["shape"]
    #     edof = Model.modelinfo["elemdof"]
    #     nodelist = Model.shape.getNodeList(inci, element_number)
    #     elementcoord = Model.shape.getNodeCoord(coord, nodelist)
    #     R = tabmat[int(inci[element_number, 2]) - 1]["RHO"]
    #     t = tabgeo[int(inci[element_number, 3] - 1)]["THICKN"]
    #     pt, wt = gauss_points(shape, intgauss)
    #     G = gravity_value
    #     W = np.zeros((nodedof, 1))
    #     W[fc_type_dof - 1, 0] = R * G
    #     force_value_vector = np.zeros((edof, 1))
    #     for ip in range(intgauss):
    #         for jp in range(intgauss):
    #             detJ = Model.shape.getdetJacobi(
    #                 np.array([pt[ip], pt[jp]]), elementcoord
    #             )
    #             N = Model.shape.getShapeFunctions(np.array([pt[ip], pt[jp]]), nodedof)
    #             force_value_vector += (
    #                 np.dot(N.transpose(), W) * t * abs(detJ) * wt[ip] * wt[jp]
    #             )
    #     force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
    #     return force_value_vector, nodelist

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
        shape_obj = Model.shape
        nodedof = Model.modelinfo["nodedof"]
        edof = Model.modelinfo["elemdof"]

        nodelist = shape_obj.getNodeList(inci, element_number)
        elementcoord = shape_obj.getNodeCoord(coord, nodelist)
        R = tabmat[int(inci[element_number, 2]) - 1]["RHO"]
        pt, wt = gauss_points(Model.modelinfo["shape"], intgauss)

        # W tem uma única componente não nula (R*G em fc_type_dof-1), então
        # Nᵀ·W = R*G * N[fc_type_dof-1, :]. Evita o produto matriz×vetor
        # e a criação de W.
        idx = fc_type_dof - 1

        get_detJ = shape_obj.getdetJacobi
        get_N = shape_obj.getShapeFunctions

        acc = np.zeros(edof)
        for ip in range(intgauss):
            wi = wt[ip]
            for jp in range(intgauss):
                wij = wi * wt[jp]
                for kp in range(intgauss):
                    xi = np.array([pt[ip], pt[jp], pt[kp]])  # criado uma vez por ponto
                    detJ = get_detJ(xi, elementcoord)
                    N = get_N(xi, nodedof)
                    acc += N[idx, :] * (abs(detJ) * wij * wt[kp])

        # constante aplicada uma única vez, no final
        force_value_vector = acc * (R * gravity_value)
        force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
        return force_value_vector, nodelist

    # def __solid_body_force_volumetric(
    #     Model,
    #     inci,
    #     coord,
    #     tabmat,
    #     tabgeo,
    #     intgauss,
    #     element_number,
    #     gravity_value,
    #     fc_type_dof,
    # ):
    #     # body force solid
    #     nodedof = Model.modelinfo["nodedof"]
    #     shape = Model.modelinfo["shape"]
    #     edof = Model.modelinfo["elemdof"]
    #     nodelist = Model.shape.getNodeList(inci, element_number)
    #     elementcoord = Model.shape.getNodeCoord(coord, nodelist)
    #     R = tabmat[int(inci[element_number, 2]) - 1]["RHO"]
    #     pt, wt = gauss_points(shape, intgauss)
    #     G = gravity_value
    #     W = np.zeros((nodedof, 1))
    #     W[fc_type_dof - 1, 0] = R * G
    #     force_value_vector = np.zeros((edof, 1))
    #     for ip in range(intgauss):
    #         for jp in range(intgauss):
    #             for kp in range(intgauss):
    #                 detJ = Model.shape.getdetJacobi(
    #                     np.array([pt[ip], pt[jp], pt[kp]]), elementcoord
    #                 )
    #                 N = Model.shape.getShapeFunctions(
    #                     np.array([pt[ip], pt[jp], pt[kp]]), nodedof
    #                 )
    #                 force_value_vector += (
    #                     np.dot(N.transpose(), W) * abs(detJ) * wt[ip] * wt[jp] * wt[kp]
    #                 )
    #     force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
    #     return force_value_vector, nodelist

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
        shape_obj = Model.shape
        nodedof = Model.modelinfo["nodedof"]
        edof = Model.modelinfo["elemdof"]

        nodelist = shape_obj.getNodeList(inci, element_number - 1)
        elementcoord = shape_obj.getNodeCoord(coord, nodelist)
        t = tabgeo[int(inci[element_number - 1, 3] - 1)]["THICKN"]

        nodelist_arr = np.asarray(nodelist)
        test = np.isin(nodelist_arr, node_list_fc, assume_unique=True)
        nodes = nodelist_arr[test]
        nodes_conec = np.flatnonzero(test)
        norm = np.zeros(2)

        if len(nodes_conec) < 2:
            nodes = np.repeat(nodes, 2)
            force_value_vector = np.zeros(len(nodes))
            return force_value_vector, nodes, norm

        idx_conec = np.array2string(nodes_conec)
        get_side = shape_obj.getSideAxis(idx_conec[1:-1])

        if fc_type == "fx":
            T = np.array([[force_value], [0.0]])
            norm[0] = 1
        elif fc_type == "fy":
            T = np.array([[0.0], [force_value]])
            norm[1] = 1
        elif fc_type == "pressure":  #  -->[+]<--
            normal = shape_obj.getNormalEdge(elementcoord, get_side)
            T = (force_value * normal).reshape(-1, 1)
            norm = (np.abs(np.asarray(normal)) >= 1e-6).astype(int)
        else:
            T = np.array([[0.0], [0.0]])

        # Carga nula (tipo não tratado ou força zero): o original integrava tudo
        # para no fim filtrar com np.nonzero e devolver um vetor vazio.
        if not T.any():
            return np.empty(0), nodes, norm

        pt, wt = gauss_points(Model.modelinfo["shape"], intgauss)

        get_iso = shape_obj.getIsoParaSide
        get_N = shape_obj.getShapeFunctions
        get_J = shape_obj.getJacobian
        get_len = shape_obj.getEdgeLength

        # T é constante: acumulo Σ N·(|detJ_e|·w) e faço Nᵀ·T uma única vez no final
        accN = 0.0
        for ip in range(intgauss):
            points = np.array(get_iso(get_side, pt[ip]))  # criado uma vez por ponto
            N = np.asarray(get_N(points, nodedof))
            J = get_J(points, elementcoord)
            accN = accN + N * (abs(get_len(J, get_side)) * wt[ip])

        force_value_vector = (accN.T @ T).ravel() * t
        force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
        return force_value_vector, nodes, norm

    # def __line_force_distribuition(
    #     Model,
    #     inci,
    #     coord,
    #     tabgeo,
    #     intgauss,
    #     node_list_fc,
    #     element_number,
    #     force_value,
    #     fc_type,
    # ):
    #     nodedof = Model.modelinfo["nodedof"]
    #     shape = Model.modelinfo["shape"]
    #     edof = Model.modelinfo["elemdof"]
    #     nodelist = Model.shape.getNodeList(inci, element_number - 1)
    #     elementcoord = Model.shape.getNodeCoord(coord, nodelist)
    #     t = tabgeo[int(inci[element_number - 1, 3] - 1)]["THICKN"]
    #     test = np.isin(nodelist, node_list_fc, assume_unique=True)
    #     nodes = np.array(nodelist)[test]
    #     nodes_conec = np.where(test == True)[0]
    #     norm = np.zeros((2))
    #     if len(nodes_conec) < 2:
    #         nodes = np.repeat(nodes, 2)
    #         force_value_vector = np.zeros((len(nodes)))
    #         pass
    #     else:
    #         idx_conec = np.array2string(nodes_conec)
    #         get_side = Model.shape.getSideAxis(idx_conec[1:-1])
    #         if fc_type == "fx":
    #             T = np.array([[force_value], [0.0]])  # force_value
    #             norm[0] = 1
    #         elif fc_type == "fy":
    #             T = np.array([[0.0], [force_value]])  # force_value
    #             norm[1] = 1
    #         elif fc_type == "pressure":  #  -->[+]<--
    #             normal = Model.shape.getNormalEdge(elementcoord, get_side)
    #             PN = force_value * normal
    #             T = np.array(PN).reshape(-1, 1)
    #             norm = np.abs(np.array(normal))
    #             norm[np.abs(norm) < 1e-6] = 0
    #             norm = (norm >= 1e-6).astype(int)  # norm/norm
    #         else:
    #             T = np.array([[0.0], [0.0]])
    #         pt, wt = gauss_points(shape, intgauss)
    #         force_value_vector = np.zeros((edof, 1))
    #         for ip in range(intgauss):
    #             points = Model.shape.getIsoParaSide(get_side, pt[ip])
    #             N = Model.shape.getShapeFunctions(np.array(points), nodedof)
    #             J = Model.shape.getJacobian(np.array(points), elementcoord)
    #             detJ_e = Model.shape.getEdgeLength(J, get_side)
    #             force_value_vector += (
    #                 np.dot(np.array(N).transpose(), T) * t * abs(detJ_e) * wt[ip]
    #             )
    #         force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
    #     return force_value_vector, nodes, norm

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

        # Carga nula (tipo não tratado ou força zero): o original integrava tudo
        # para no fim filtrar com np.nonzero e devolver um vetor vazio.
        if not T.any():
            return np.empty(0), nodes, norm

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
        force_value_vector[np.abs(force_value_vector) < 1e-6] = 0
        force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
        return force_value_vector, nodes, norm

    # def __surf_force_distribuition(
    #     Model,
    #     inci,
    #     coord,
    #     tabgeo,
    #     intgauss,
    #     node_list_fc,
    #     element_number,
    #     force_value,
    #     fc_type,
    # ):
    #     nodedof = Model.modelinfo["nodedof"]
    #     shape = Model.modelinfo["shape"]
    #     edof = Model.modelinfo["elemdof"]
    #     nodelist = Model.shape.getNodeList(inci, element_number - 1)
    #     elementcoord = Model.shape.getNodeCoord(coord, nodelist)
    #     test = np.isin(nodelist, node_list_fc, assume_unique=True)
    #     nodes = np.array(nodelist)[test]
    #     nodes_conec = np.where(test == True)[0]
    #     norm = np.zeros((3))
    #     if len(nodes_conec) < 3:
    #         nodes = np.repeat(nodes, 3)
    #         force_value_vector = np.zeros((len(nodes)))
    #         pass
    #     else:
    #         idx_conec = np.array2string(nodes_conec)
    #         get_side = Model.shape.getSideAxis(idx_conec[1:-1])
    #         if fc_type == "fx":
    #             T = np.array([[force_value], [0.0], [0.0]])  # force_value
    #             norm[0] = 1
    #         elif fc_type == "fy":
    #             T = np.array([[0.0], [force_value], [0.0]])  # force_value
    #             norm[1] = 1
    #         elif fc_type == "fz":
    #             T = np.array([[0.0], [0.0], [force_value]])  # force_value
    #             norm[2] = 1
    #         elif fc_type == "pressure":  #  -->[+]<--
    #             normal = Model.shape.getNormalFace(
    #                 elementcoord, get_side
    #             )  # [pt[ip], pt[jp]]
    #             PN = force_value * normal
    #             # T = np.sign(force_value)*(np.abs(np.array(PN).reshape(-1, 1)))
    #             T = np.array(PN).reshape(-1, 1)
    #             norm = np.abs(np.array(normal))
    #             norm[np.abs(norm) < 1e-6] = 0
    #             norm = (norm >= 1e-6).astype(int)  # norm/norm
    #         else:
    #             T = np.array([[0.0], [0.0], [0.0]])
    #         pt, wt = gauss_points(shape, intgauss)
    #         force_value_vector = np.zeros((edof, 1))
    #         for ip in range(intgauss):
    #             for jp in range(intgauss):
    #                 points = Model.shape.getIsoParaSide(
    #                     get_side, [pt[ip], pt[jp]]
    #                 )  # [pt[ip], pt[jp]]
    #                 N = Model.shape.getShapeFunctions(np.array(points), nodedof)
    #                 # J = Model.shape.getJacobian(np.array(points), elementcoord)
    #                 detJ_a = Model.shape.getAreaLength(get_side, elementcoord)
    #                 force_value_vector += (
    #                     np.dot(np.array(N).transpose(), T)
    #                     * abs(detJ_a)
    #                     * wt[ip]
    #                     * wt[jp]
    #                 )

    #         # if fc_type == "pressure":
    #         #     force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
    #         #     force_value_vector[np.abs(force_value_vector) < 1e-11] = 0
    #         # else:
    #         force_value_vector[np.abs(force_value_vector) < 1e-6] = 0
    #         force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
    #     return force_value_vector, nodes, norm

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
        shape_obj = Model.shape
        nodedof = Model.modelinfo["nodedof"]

        nodelist = shape_obj.getNodeList(inci, element_number - 1)
        elementcoord = shape_obj.getNodeCoord(coord, nodelist)

        nodelist_arr = np.asarray(nodelist)
        test = np.isin(nodelist_arr, node_list_fc, assume_unique=True)
        nodes = nodelist_arr[test]
        idx_conec = np.flatnonzero(test)
        norm = np.zeros(3)

        if len(idx_conec) < 2:
            nodes = np.repeat(nodes, 2)
            force_value_vector = np.zeros(len(nodes))
            return force_value_vector, nodes, norm

        if fc_type == "fx":
            T = np.array([[force_value], [0.0], [0.0], [0.0]])
            norm[0] = 1
        elif fc_type == "fy":
            T = np.array([[0.0], [force_value], [0.0], [0.0]])
            norm[1] = 1
        elif fc_type == "fz":
            T = np.array([[0.0], [0.0], [force_value], [0.0]])
            norm[2] = 1
        elif fc_type == "tx":
            T = np.array([[0.0], [0.0], [0.0], [force_value]])
            norm[0] = 1
        elif fc_type == "pressure":  #  -->[+]<--
            noi = idx_conec[0]
            noj = idx_conec[1]
            dx = abs(elementcoord[noj, 0] - elementcoord[noi, 0])
            dy = abs(elementcoord[noj, 1] - elementcoord[noi, 1])
            L = np.sqrt(dx**2 + dy**2)
            norm[0] = dy / L
            norm[1] = dx / L
            T = -force_value * norm.reshape(-1, 1)
        else:
            T = np.zeros((4, 1))

        # Carga nula (tipo não tratado ou força zero): o original integrava tudo
        # para no fim filtrar com np.nonzero e devolver um vetor vazio.
        if not T.any():
            return np.empty(0), nodes, norm

        get_side = shape_obj.getSideAxis(np.array2string(idx_conec)[1:-1])
        pt, wt = gauss_points(Model.modelinfo["shape"], intgauss)

        get_iso = shape_obj.getIsoParaSide
        get_N = shape_obj.getShapeFunctions
        get_J = shape_obj.getJacobian
        get_len = shape_obj.getEdgeLength

        # T é constante: acumulo Σ N·(|detJ_e|·w) e faço Nᵀ·T uma única vez no final
        accN = 0.0
        for ip in range(intgauss):
            points = np.array(get_iso(get_side, pt[ip]))  # criado uma vez por ponto
            N = np.asarray(get_N(points, nodedof))
            J = get_J(points, elementcoord)
            accN = accN + N * (abs(get_len(J, get_side)) * wt[ip])

        force_value_vector = (accN.T @ T).ravel()
        force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
        return force_value_vector, nodes, norm

    # def __line_beam_distribuition(
    #     Model,
    #     inci,
    #     coord,
    #     tabgeo,
    #     intgauss,
    #     node_list_fc,
    #     element_number,
    #     force_value,
    #     fc_type,
    # ):
    #     nodedof = Model.modelinfo["nodedof"]
    #     shape = Model.modelinfo["shape"]
    #     edof = Model.modelinfo["elemdof"]
    #     nodelist = Model.shape.getNodeList(inci, element_number - 1)
    #     elementcoord = Model.shape.getNodeCoord(coord, nodelist)
    #     test = np.isin(nodelist, node_list_fc, assume_unique=True)
    #     nodes = np.array(nodelist)[test]
    #     idx_conec = np.where(test == True)[0]
    #     norm = np.zeros((3))
    #     if len(idx_conec) < 2:
    #         nodes = np.repeat(nodes, 2)
    #         force_value_vector = np.zeros((len(nodes)))
    #         pass
    #     else:
    #         if fc_type == "fx":
    #             T = np.array([[force_value], [0.0], [0.0], [0.0]])  # force_value
    #             norm[0] = 1
    #         elif fc_type == "fy":
    #             T = np.array([[0.0], [force_value], [0.0], [0.0]])  # force_value
    #             norm[1] = 1
    #         elif fc_type == "fz":
    #             T = np.array([[0.0], [0.0], [force_value], [0.0]])  # force_value
    #             norm[2] = 1
    #         elif fc_type == "tx":
    #             T = np.array([[0.0], [0.0], [0.0], [force_value]])  # force_value
    #             norm[0] = 1
    #         elif fc_type == "pressure":  #  -->[+]<--
    #             noi = idx_conec[0]
    #             noj = idx_conec[1]
    #             dx = abs(elementcoord[noj, 0] - elementcoord[noi, 0])
    #             dy = abs(elementcoord[noj, 1] - elementcoord[noi, 1])
    #             L = np.sqrt(dx**2 + dy**2)
    #             # tx = (-dy / L) * force_value
    #             # ty = (dx / L) * force_value
    #             # T = np.array([[tx], [ty]])  # force_value
    #             norm[0] = dy / L
    #             norm[1] = dx / L
    #             T = -1 * (np.array([norm]).T) * force_value  # (np.array([norm]).T)*T
    #         else:
    #             T = np.array([[0.0], [0.0], [0.0], [0.0]])
    #         idx_conec = np.array2string(idx_conec)
    #         get_side = Model.shape.getSideAxis(idx_conec[1:-1])
    #         pt, wt = gauss_points(shape, intgauss)
    #         force_value_vector = np.zeros((edof, 1))
    #         for ip in range(intgauss):
    #             points = Model.shape.getIsoParaSide(get_side, pt[ip])
    #             N = Model.shape.getShapeFunctions(np.array(points), nodedof)
    #             J = Model.shape.getJacobian(np.array(points), elementcoord)
    #             detJ_e = Model.shape.getEdgeLength(J, get_side)
    #             force_value_vector += (
    #                 np.dot(np.array(N).transpose(), T) * abs(detJ_e) * wt[ip]
    #             )
    #         force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
    #     return force_value_vector, nodes, norm
