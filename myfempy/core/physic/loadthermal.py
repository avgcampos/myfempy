from __future__ import annotations

import numpy as np
from scipy.special import roots_legendre

from myfempy.core.physic.physics import Physics
from myfempy.core.utilities import (gauss_points, get_elemen_from_nodelist,
                                    get_nodes_from_list, poly_area)

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


class LoadThermal(Physics):
    """Thermal Load Class <ConcreteClassService>"""

    def getLoadApply(Model, forcelist):
        forcenodeaply = np.zeros((1, 4))
        if forcelist["TYPE"] == "thermalconvection":
            fapp = LoadThermal.ThermalConvecNode(Model, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        if forcelist["TYPE"] == "heatfluxedge":
            fapp = LoadThermal.HeatFluxEdge(Model, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        elif forcelist["TYPE"] == "heatfluxsurf":
            fapp = LoadThermal.HeatFluxSurf(Model, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        elif forcelist["TYPE"] == "heatgeneration":
            fapp = LoadThermal.HeatGeneration(Model, forcelist)
            forcenodeaply = np.append(forcenodeaply, fapp, axis=0)
        else:
            pass

        forcenodeaply = forcenodeaply[1::][::]
        return forcenodeaply

    def getUpdateMatrix(Model, matrix, loadaply):
        addConv = np.where(loadaply[:, 1] == 15)
        if addConv[0].size:
            addLoad = loadaply[addConv, :][0]
            addLoad = addLoad[~np.isclose(addLoad[:, 2], 0)]
            matrix["stiffness"] = Model.element.getUpdateMatrix(
                Model, matrix["stiffness"], addLoad
            )
        return matrix

    def getUpdateLoad(self):
        return None

    def ThermalConvecNode(Model, forcelist):
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

        # Montagem direta, sem loop e sem np.append
        forcenodedof = np.empty((nodes.size, 4))
        forcenodedof[:, 0] = nodes.astype(np.int64)
        forcenodedof[:, 1] = Model.modelinfo["dofs"]["f"][forcelist["DOF"]]
        forcenodedof[:, 2] = float(forcelist["VAL"])
        forcenodedof[:, 3] = int(forcelist["STEP"])
        return forcenodedof

    def HeatGeneration(Model, forcelist):
        type_shape = Model.modelinfo["type_shape"]

        # Decisão do tipo de elemento feita UMA vez, fora do loop
        if type_shape in ("hexa8", "tetr4"):
            elem_func = LoadThermal.__solid_body_force_volumetric
        elif type_shape in ("tria3", "tria6", "quad4", "quad8"):
            elem_func = LoadThermal.__plane_body_force_volumetric
        else:
            elem_func = None

        heatgen = float(forcelist["VAL"])
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
                heatgen,
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


    def HeatFluxEdge(Model, forcelist):
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

        fc_type_dof = Model.modelinfo["dofs"]["f"][fc_type]  # constante: uma vez
        line_force = LoadThermal.__line_force_distribuition

        nodes_list = []
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
            nodes = np.asarray(nodes)
            # n linhas por elemento, como no range(len(nodelist)) original
            nodes_list.append(nodes)
            vals_list.append(
                np.asarray(force_value_vector, dtype=float).ravel()[: nodes.size]
            )

        # Montagem única, sem np.append dentro do loop
        nodes_all = np.concatenate(nodes_list).astype(np.int64)

        forcenodedof = np.empty((nodes_all.size, 4))
        forcenodedof[:, 0] = nodes_all
        forcenodedof[:, 1] = fc_type_dof
        forcenodedof[:, 2] = np.concatenate(vals_list)
        forcenodedof[:, 3] = step
        return forcenodedof

    def HeatFluxSurf(Model, forcelist):
        nodelist = [
            forcelist["DIR"],
            forcelist["LOCX"],
            forcelist["LOCY"],
            forcelist["LOCZ"],
            forcelist["TAG"],
            forcelist["MESHNODE"],
        ]
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

        fc_type_dof = Model.modelinfo["dofs"]["f"][fc_type]  # constante: uma vez
        surf_force = LoadThermal.__surf_force_distribuition

        nodes_list = []
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
            nodes = np.asarray(nodes)
            nodes_list.append(nodes)
            vals_list.append(
                np.asarray(force_value_vector, dtype=float).ravel()[: nodes.size]
            )

        nodes_all = np.concatenate(nodes_list).astype(np.int64)

        forcenodedof = np.empty((nodes_all.size, 4))
        forcenodedof[:, 0] = nodes_all
        forcenodedof[:, 1] = fc_type_dof
        forcenodedof[:, 2] = np.concatenate(vals_list)
        forcenodedof[:, 3] = step
        return forcenodedof

    def __plane_body_force_volumetric(
        Model,
        inci,
        coord,
        tabmat,
        tabgeo,
        intgauss,
        element_number,
        heat_gen,
        fc_type_dof,
    ):
        shape_obj = Model.shape
        nodedof = Model.modelinfo["nodedof"]
        edof = Model.modelinfo["elemdof"]

        nodelist = shape_obj.getNodeList(inci, element_number)

        # Geração nula: o original integrava tudo para o np.nonzero devolver vazio
        if heat_gen == 0:
            return np.empty(0), nodelist

        elementcoord = shape_obj.getNodeCoord(coord, nodelist)
        t = tabgeo[int(inci[element_number, 3] - 1)]["THICKN"]
        pt, wt = gauss_points(Model.modelinfo["shape"], intgauss)

        get_detJ = shape_obj.getdetJacobi
        get_N = shape_obj.getShapeFunctions

        # Q e t são constantes: acumulo só Σ N·|detJ|·w e multiplico no final
        acc = np.zeros(edof)
        for ip in range(intgauss):
            wi = wt[ip]
            for jp in range(intgauss):
                xi = np.array([pt[ip], pt[jp]])  # criado uma vez por ponto
                detJ = get_detJ(xi, elementcoord)
                N = np.asarray(get_N(xi, nodedof))
                acc += N.ravel() * (abs(detJ) * wi * wt[jp])

        force_value_vector = acc * (heat_gen * t)
        force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
        return force_value_vector, nodelist

    def __solid_body_force_volumetric(
        Model,
        inci,
        coord,
        tabmat,
        tabgeo,
        intgauss,
        element_number,
        heat_gen,
        fc_type_dof,
    ):
        shape_obj = Model.shape
        nodedof = Model.modelinfo["nodedof"]
        edof = Model.modelinfo["elemdof"]

        nodelist = shape_obj.getNodeList(inci, element_number)

        if heat_gen == 0:
            return np.empty(0), nodelist

        elementcoord = shape_obj.getNodeCoord(coord, nodelist)
        pt, wt = gauss_points(Model.modelinfo["shape"], intgauss)

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
                    N = np.asarray(get_N(xi, nodedof))
                    acc += N.ravel() * (abs(detJ) * wij * wt[kp])

        force_value_vector = acc * heat_gen
        force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
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
        shape_obj = Model.shape
        nodedof = Model.modelinfo["nodedof"]
        edof = Model.modelinfo["elemdof"]

        nodelist = shape_obj.getNodeList(inci, element_number - 1)
        elementcoord = shape_obj.getNodeCoord(coord, nodelist)
        t = tabgeo[int(inci[element_number - 1, 3] - 1)]["THICKN"]

        nodelist_arr = np.asarray(nodelist)
        test = np.isin(nodelist_arr, node_list_fc, assume_unique=True)
        nodes = nodelist_arr[test]
        idx_conec = np.flatnonzero(test)
        norm = np.zeros(2)

        if len(idx_conec) < 2:
            nodes = np.repeat(nodes, 2)
            return np.zeros(len(nodes)), nodes, norm

        if force_value == 0:
            return np.zeros(len(nodes)), nodes, norm

        get_side = shape_obj.getSideAxis(np.array2string(idx_conec)[1:-1])
        pt, wt = gauss_points(Model.modelinfo["shape"], intgauss)

        get_iso = shape_obj.getIsoParaSide
        get_N = shape_obj.getShapeFunctions
        get_J = shape_obj.getJacobian
        get_len = shape_obj.getEdgeLength

        # q e t são constantes: acumulo Σ N·|detJ_e|·w e multiplico no final
        acc = np.zeros(edof)
        for ip in range(intgauss):
            points = np.array(get_iso(get_side, pt[ip]))  # criado uma vez por ponto
            N = np.asarray(get_N(points, nodedof))
            J = get_J(points, elementcoord)
            acc += N.ravel() * (abs(get_len(J, get_side)) * wt[ip])

        force_value_vector = acc * (force_value * t)
        force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
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
        edof = Model.modelinfo["elemdof"]

        nodelist = shape_obj.getNodeList(inci, element_number - 1)
        elementcoord = shape_obj.getNodeCoord(coord, nodelist)

        nodelist_arr = np.asarray(nodelist)
        test = np.isin(nodelist_arr, node_list_fc, assume_unique=True)
        nodes = nodelist_arr[test]
        nodes_conec = np.flatnonzero(test)
        norm = np.zeros(3)

        if len(nodes_conec) < 3:
            nodes = np.repeat(nodes, 2)
            return np.zeros(len(nodes)), nodes, norm

        if force_value == 0:
            return np.zeros(len(nodes)), nodes, norm

        get_side = shape_obj.getSideAxis(np.array2string(nodes_conec)[1:-1])
        pt, wt = gauss_points(Model.modelinfo["shape"], intgauss)

        # A área da face não depende do ponto de Gauss: calculada uma vez
        detJ_a = abs(shape_obj.getAreaLength(get_side, elementcoord))

        get_iso = shape_obj.getIsoParaSide
        get_N = shape_obj.getShapeFunctions

        acc = np.zeros(edof)
        for ip in range(intgauss):
            wi = wt[ip]
            for jp in range(intgauss):
                points = np.array(get_iso(get_side, [pt[ip], pt[jp]]))
                N = np.asarray(get_N(points, nodedof))
                acc += N.ravel() * (wi * wt[jp])

        force_value_vector = acc * (force_value * detJ_a)
        force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
        return force_value_vector, nodes, norm


    # def ThermalConvecNode(Model, forcelist):
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

    # def HeatGeneration(Model, forcelist):
    #     forcenodedof = np.zeros((1, 4))
    #     heatgen = float(forcelist["VAL"])
    #     inci = Model.inci
    #     coord = Model.coord
    #     tabmat = Model.tabmat
    #     tabgeo = Model.tabgeo
    #     intgauss = Model.intgauss
    #     fc_type_dof = Model.modelinfo["dofs"]["f"][forcelist["DOF"]]
    #     for ee in range(inci.shape[0]):
    #         if (
    #             Model.modelinfo["type_shape"] == "hexa8"
    #             or Model.modelinfo["type_shape"] == "tetr4"
    #         ):
    #             force_value_vector, nodelist = (
    #                 LoadThermal.__solid_body_force_volumetric(
    #                     Model,
    #                     inci,
    #                     coord,
    #                     tabmat,
    #                     tabgeo,
    #                     intgauss,
    #                     ee,
    #                     heatgen,
    #                     fc_type_dof,
    #                 )
    #             )
    #         elif (
    #             Model.modelinfo["type_shape"] == "tria3"
    #             or Model.modelinfo["type_shape"] == "tria6"
    #             or Model.modelinfo["type_shape"] == "quad4"
    #             or Model.modelinfo["type_shape"] == "quad8"
    #         ):
    #             force_value_vector, nodelist = (
    #                 LoadThermal.__plane_body_force_volumetric(
    #                     Model,
    #                     inci,
    #                     coord,
    #                     tabmat,
    #                     tabgeo,
    #                     intgauss,
    #                     ee,
    #                     heatgen,
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

    # def HeatFluxEdge(Model, forcelist):
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

    #         (
    #             force_value_vector,
    #             nodelist,
    #             norm,
    #         ) = LoadThermal.__line_force_distribuition(
    #             Model,
    #             inci,
    #             coord,
    #             tabgeo,
    #             intgauss,
    #             node_list_fc,
    #             elmlist[ee],
    #             force_value,
    #             fc_type,
    #         )

    #         fc_type_dof = Model.modelinfo["dofs"]["f"][forcelist["DOF"]] * np.ones_like(
    #             nodelist
    #         )

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

    # def HeatFluxSurf(Model, forcelist):
    #     nodelist = [
    #         forcelist["DIR"],
    #         forcelist["LOCX"],
    #         forcelist["LOCY"],
    #         forcelist["LOCZ"],
    #         forcelist["TAG"],
    #         forcelist["MESHNODE"],
    #     ]
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
    #         (
    #             force_value_vector,
    #             nodelist,
    #             norm,
    #         ) = LoadThermal.__surf_force_distribuition(
    #             Model,
    #             inci,
    #             coord,
    #             tabgeo,
    #             intgauss,
    #             node_list_fc,
    #             elmlist[ee],
    #             force_value,
    #             fc_type,
    #         )

    #         fc_type_dof = Model.modelinfo["dofs"]["f"][forcelist["DOF"]] * np.ones_like(
    #             nodelist
    #         )

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

    # def __plane_body_force_volumetric(
    #     Model,
    #     inci,
    #     coord,
    #     tabmat,
    #     tabgeo,
    #     intgauss,
    #     element_number,
    #     heat_gen,
    #     fc_type_dof,
    # ):
    #     nodedof = Model.modelinfo["nodedof"]
    #     shape = Model.modelinfo["shape"]
    #     edof = Model.modelinfo["elemdof"]
    #     nodelist = Model.shape.getNodeList(inci, element_number)
    #     elementcoord = Model.shape.getNodeCoord(coord, nodelist)
    #     t = tabgeo[int(inci[element_number, 3] - 1)]["THICKN"]
    #     pt, wt = gauss_points(shape, intgauss)
    #     Q = heat_gen
    #     force_value_vector = np.zeros((edof, 1))
    #     for ip in range(intgauss):
    #         for jp in range(intgauss):
    #             detJ = Model.shape.getdetJacobi(
    #                 np.array([pt[ip], pt[jp]]), elementcoord
    #             )
    #             N = Model.shape.getShapeFunctions(np.array([pt[ip], pt[jp]]), nodedof)
    #             force_value_vector += (
    #                 np.dot(N.transpose(), Q) * t * abs(detJ) * wt[ip] * wt[jp]
    #             )
    #     force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
    #     return force_value_vector, nodelist

    # def __solid_body_force_volumetric(
    #     Model,
    #     inci,
    #     coord,
    #     tabmat,
    #     tabgeo,
    #     intgauss,
    #     element_number,
    #     heat_gen,
    #     fc_type_dof,
    # ):
    #     nodedof = Model.modelinfo["nodedof"]
    #     shape = Model.modelinfo["shape"]
    #     edof = Model.modelinfo["elemdof"]
    #     nodelist = Model.shape.getNodeList(inci, element_number)
    #     elementcoord = Model.shape.getNodeCoord(coord, nodelist)
    #     pt, wt = gauss_points(shape, intgauss)
    #     Q = heat_gen
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
    #                     np.dot(N.transpose(), Q) * abs(detJ) * wt[ip] * wt[jp] * wt[kp]
    #                 )
    #     force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
    #     return force_value_vector, nodelist

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
    #     idx_conec = np.where(test == True)[0]
    #     norm = np.zeros((2))
    #     q = force_value
    #     if len(idx_conec) < 2:
    #         nodes = np.repeat(nodes, 2)
    #         force_value_vector = np.zeros((len(nodes)))
    #         pass
    #     else:
    #         force_value_vector = np.zeros((len(nodes)))
    #         if force_value == 0:
    #             pass
    #         else:
    #             idx_conec = np.array2string(idx_conec)
    #             get_side = Model.shape.getSideAxis(idx_conec[1:-1])
    #             pt, wt = gauss_points(shape, intgauss)
    #             force_value_vector = np.zeros((edof, 1))
    #             for ip in range(intgauss):
    #                 points = Model.shape.getIsoParaSide(get_side, pt[ip])
    #                 N = Model.shape.getShapeFunctions(np.array(points), nodedof)
    #                 J = Model.shape.getJacobian(np.array(points), elementcoord)
    #                 detJ_e = Model.shape.getEdgeLength(J, get_side)
    #                 force_value_vector += (
    #                     np.dot(N.transpose(), q) * t * abs(detJ_e) * wt[ip]
    #                 )
    #             force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
    #     return force_value_vector, nodes, norm

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
    #     q = force_value
    #     if len(nodes_conec) < 3:
    #         nodes = np.repeat(nodes, 2)
    #         force_value_vector = np.zeros((len(nodes)))
    #         pass
    #     else:
    #         force_value_vector = np.zeros((len(nodes)))
    #         if force_value == 0:
    #             pass
    #         else:
    #             idx_conec = np.array2string(nodes_conec)
    #             get_side = Model.shape.getSideAxis(idx_conec[1:-1])
    #             pt, wt = gauss_points(shape, intgauss)
    #             force_value_vector = np.zeros((edof, 1))
    #             for ip in range(intgauss):
    #                 for jp in range(intgauss):
    #                     points = Model.shape.getIsoParaSide(
    #                         get_side, [pt[ip], pt[jp]]
    #                     )  # [pt[ip], pt[jp]]
    #                     N = Model.shape.getShapeFunctions(np.array(points), nodedof)
    #                     detJ_a = Model.shape.getAreaLength(get_side, elementcoord)
    #                     force_value_vector += (
    #                         np.dot(N.transpose(), q) * abs(detJ_a) * wt[ip] * wt[jp]
    #                     )
    #             force_value_vector = force_value_vector[np.nonzero(force_value_vector)]
    #     return force_value_vector, nodes, norm
