#!/usr/bin/env python
__doc__ = """
calculate loads and constrains
"""
import numpy as np


def get_forces(modelinfo, flist):
    """_summary_

    Arguments:
        modelinfo -- _description_
        flist -- _description_

    Returns:
        _description_
    """
    forcenodeaply = np.zeros((1, 4))
    for k in range(int(len(np.unique([flist[:, 7]])))):
        forcelist = flist[np.where(flist[:, 7] == str(k + 1))[0], :]
        node_list_fc = np.array([])
        for i in range(len(forcelist)):
            # ----- SEEKERS IN LOC -----
            if forcelist[i][3] == "lengthx":
                coord_0 = float(forcelist[i, 5])
                coord_1 = float(forcelist[i, 6])
                node_list_fc = modelinfo["coord"][
                    np.where(
                        (modelinfo["coord"][:, 1] >= coord_0)
                        & (modelinfo["coord"][:, 1] <= coord_1)
                    ),
                    0,
                ][0]
                dir_fc = "x"
            elif forcelist[i][3] == "lengthy":
                coord_0 = float(forcelist[i, 5])
                coord_1 = float(forcelist[i, 6])
                node_list_fc = modelinfo["coord"][
                    np.where(
                        (modelinfo["coord"][:, 2] >= coord_0)
                        & (modelinfo["coord"][:, 2] <= coord_1)
                    ),
                    0,
                ][0]
                dir_fc = "y"
            elif forcelist[i][3] == "lengthz":
                coord_0 = float(forcelist[i, 5])
                coord_1 = float(forcelist[i, 6])
                node_list_fc = modelinfo["coord"][
                    np.where(
                        (modelinfo["coord"][:, 3] >= coord_0)
                        & (modelinfo["coord"][:, 3] <= coord_1)
                    ),
                    0,
                ][0]
                dir_fc = "z"
            elif forcelist[i][3] == "edgex":
                # sys.path.append('../lib')
                from myfempy.felib.physics.getnode import search_edgex

                edge_coordX = float(forcelist[i, 4])
                if float(forcelist[i, 5]) == 999:
                    dir_fc = "x_y"
                    coord_fc = (
                        modelinfo["coord"][
                            np.where(
                                modelinfo["coord"][:, 3] == float(forcelist[i, 6])
                            ),
                            :,
                        ]
                    )[0]
                elif float(forcelist[i, 6]) == 999:
                    dir_fc = "x_z"
                    coord_fc = (
                        modelinfo["coord"][
                            np.where(
                                modelinfo["coord"][:, 2] == float(forcelist[i, 5])
                            ),
                            :,
                        ]
                    )[0]
                node_list_fc = search_edgex(edge_coordX, coord_fc, 2e-3)
            elif forcelist[i][3] == "edgey":
                from myfempy.felib.physics.getnode import search_edgey

                edge_coordY = float(forcelist[i, 5])
                if float(forcelist[i, 4]) == 999:
                    dir_fc = "y_x"
                    coord_fc = (
                        modelinfo["coord"][
                            np.where(
                                modelinfo["coord"][:, 3] == float(forcelist[i, 6])
                            ),
                            :,
                        ]
                    )[0]
                elif float(forcelist[i, 6]) == 999:
                    dir_fc = "y_z"
                    coord_fc = (
                        modelinfo["coord"][
                            np.where(
                                modelinfo["coord"][:, 1] == float(forcelist[i, 4])
                            ),
                            :,
                        ]
                    )[0]
                node_list_fc = search_edgey(edge_coordY, coord_fc, 2e-3)
            elif forcelist[i][3] == "edgez":
                from myfempy.felib.physics.getnode import search_edgez

                edge_coordZ = float(forcelist[i, 6])
                if float(forcelist[i, 4]) == 999:
                    dir_fc = "z_x"
                    coord_fc = (
                        modelinfo["coord"][
                            np.where(
                                modelinfo["coord"][:, 2] == float(forcelist[i, 5])
                            ),
                            :,
                        ]
                    )[0]
                elif float(forcelist[i, 5]) == 999:
                    dir_fc = "z_x"
                    coord_fc = (
                        modelinfo["coord"][
                            np.where(
                                modelinfo["coord"][:, 1] == float(forcelist[i, 4])
                            ),
                            :,
                        ]
                    )[0]
                node_list_fc = search_edgez(edge_coordZ, coord_fc, 2e-3)
            elif forcelist[i][3] == "surfxy":
                from myfempy.felib.physics.getnode import search_surfxy

                orthg_coordZ = float(forcelist[i, 6])
                node_list_fc = search_surfxy(orthg_coordZ, modelinfo["coord"], 2e-3)
                dir_fc = "z"
            elif forcelist[i][3] == "surfyz":
                from myfempy.felib.physics.getnode import search_surfyz

                orthg_coordX = float(forcelist[i, 4])
                node_list_fc = search_surfyz(orthg_coordX, modelinfo["coord"], 2e-3)
                dir_fc = "x"
            elif forcelist[i][3] == "surfzx":
                from myfempy.felib.physics.getnode import search_surfzx

                orthg_coordY = float(forcelist[i, 5])
                node_list_fc = search_surfzx(orthg_coordY, modelinfo["coord"], 2e-3)
                dir_fc = "y"
            elif forcelist[i][3] == "node":
                from myfempy.felib.physics.getnode import search_nodexyz

                node_coordX = float(forcelist[i, 4])
                node_coordY = float(forcelist[i, 5])
                node_coordZ = float(forcelist[i, 6])
                node_list_fc = search_nodexyz(
                    node_coordX, node_coordY, node_coordZ, modelinfo["coord"], 2e-3
                )
            # ----- SEEKERS IN TAG -----
            elif forcelist[i][3] == "point":
                if forcelist[i][1] == "fx":
                    dir_fc = "x"
                elif forcelist[i][1] == "fy":
                    dir_fc = "y"
                elif forcelist[i][1] == "fz":
                    dir_fc = "z"
                elif forcelist[i][1] == "tx":
                    dir_fc = "x"
                elif forcelist[i][1] == "ty":
                    dir_fc = "y"
                elif forcelist[i][1] == "tz":
                    dir_fc = "z"
                else:
                    pass
                node_list_fc = modelinfo["regions"][0][1][int(forcelist[i][8]) - 1][1][
                    :
                ]
            elif forcelist[i][3] == "edge":
                if forcelist[i][1] == "fx":
                    dir_fc = "x"
                elif forcelist[i][1] == "fy":
                    dir_fc = "y"
                elif forcelist[i][1] == "fz":
                    dir_fc = "z"
                elif forcelist[i][1] == "tx":
                    dir_fc = "x"
                elif forcelist[i][1] == "ty":
                    dir_fc = "y"
                elif forcelist[i][1] == "tz":
                    dir_fc = "z"
                else:
                    pass
                node_list_fc = modelinfo["regions"][1][1][int(forcelist[i][8]) - 1][1][
                    :
                ]
            elif forcelist[i][3] == "surf":
                if forcelist[i][1] == "fx":
                    dir_fc = "x"
                elif forcelist[i][1] == "fy":
                    dir_fc = "y"
                elif forcelist[i][1] == "fz":
                    dir_fc = "z"
                elif forcelist[i][1] == "tx":
                    dir_fc = "x"
                elif forcelist[i][1] == "ty":
                    dir_fc = "y"
                elif forcelist[i][1] == "tz":
                    dir_fc = "z"
                else:
                    pass
                node_list_fc = modelinfo["regions"][2][1][int(forcelist[i][8]) - 1][1][
                    :
                ]
            else:
                print("input erro: force_opt don't defined")
            # ----- TYPE OF LOADS -----
            if forcelist[i][0] == "forcenode":
                force_value_vector = np.ones_like(node_list_fc) * float(forcelist[i, 2])
                if forcelist[i][1] == "fx":
                    fc_type_dof = 1 * np.ones_like(node_list_fc)
                elif forcelist[i][1] == "fy":
                    fc_type_dof = 2 * np.ones_like(node_list_fc)
                elif forcelist[i][1] == "fz":
                    fc_type_dof = 3 * np.ones_like(node_list_fc)
                elif forcelist[i][1] == "tx":
                    fc_type_dof = 4 * np.ones_like(node_list_fc)
                elif forcelist[i][1] == "ty":
                    fc_type_dof = 5 * np.ones_like(node_list_fc)
                elif forcelist[i][1] == "tz":
                    fc_type_dof = 6 * np.ones_like(node_list_fc)
                elif forcelist[i][1] == "masspoint":
                    fc_type_dof = 15 * np.ones_like(node_list_fc)
                elif forcelist[i][1] == "spring2ground":
                    fc_type_dof = 16 * np.ones_like(node_list_fc)
                elif forcelist[i][1] == "damper2ground":
                    fc_type_dof = 17 * np.ones_like(node_list_fc)
                else:
                    print("input erro: force_opt_dir don't defined")
                for j in range(len(node_list_fc)):
                    fcdef = np.array(
                        [
                            [
                                int(node_list_fc[j]),
                                fc_type_dof[j],
                                force_value_vector[j],
                                int(forcelist[i][7]),
                            ]
                        ]
                    )
                    forcenodeaply = np.append(forcenodeaply, fcdef, axis=0)
            elif forcelist[i][0] == "forcebeam":
                from myfempy.felib.physics.force2node import force_beam

                line_fc = int(float(forcelist[i, 8]))
                force_value = float(forcelist[i, 2])
                force_dirc = forcelist[i, 1]
                force_value_vector, fc_type_dof = force_beam(
                    modelinfo, force_value, force_dirc, dir_fc, node_list_fc, line_fc
                )
                fc_type_dof = np.tile(fc_type_dof, len(node_list_fc))
                node_list_fc = np.repeat(node_list_fc, modelinfo["nodedof"][0], axis=0)
                for j in range(len(node_list_fc)):
                    fcdef = np.array(
                        [
                            [
                                int(node_list_fc[j]),
                                fc_type_dof[j],
                                force_value_vector[j, 0],
                                int(forcelist[i][7]),
                            ]
                        ]
                    )
                    forcenodeaply = np.append(forcenodeaply, fcdef, axis=0)
            elif forcelist[i][0] == "forceedge":
                from myfempy.felib.physics.force2node import force_edge

                force_value = float(forcelist[i, 2])
                force_dirc = forcelist[i, 1]
                force_value_vector, fc_type_dof = force_edge(
                    modelinfo, force_value, force_dirc, node_list_fc, dir_fc
                )
                fc_type_dof = np.tile(fc_type_dof, len(node_list_fc))
                node_list_fc = np.repeat(node_list_fc, modelinfo["nodedof"][0], axis=0)
                for j in range(len(node_list_fc)):
                    fcdef = np.array(
                        [
                            [
                                int(node_list_fc[j]),
                                fc_type_dof[j],
                                force_value_vector[j, 0],
                                int(forcelist[i][7]),
                            ]
                        ]
                    )
                    forcenodeaply = np.append(forcenodeaply, fcdef, axis=0)
            elif forcelist[i][0] == "forcesurf":
                from myfempy.felib.physics.force2node import force_surf

                force_value = float(forcelist[i, 2])
                force_dirc = forcelist[i, 1]
                force_value_vector, fc_type_dof = force_surf(
                    modelinfo, force_value, force_dirc, node_list_fc, dir_fc
                )
                fc_type_dof = np.tile(fc_type_dof, len(node_list_fc))
                node_list_fc = np.repeat(node_list_fc, modelinfo["nodedof"][0], axis=0)
                for j in range(len(node_list_fc)):
                    fcdef = np.array(
                        [
                            [
                                int(node_list_fc[j]),
                                fc_type_dof[j],
                                force_value_vector[j, 0],
                                int(forcelist[i][7]),
                            ]
                        ]
                    )
                    forcenodeaply = np.append(forcenodeaply, fcdef, axis=0)
            # elif eval('usrlog.force_typ_'+str(i)) == "gravity":
            # elif eval('usrlog.force_typ_'+str(i)) == "pressure":
            else:
                print("input erro: force_typ don't defined")
    forces = forcenodeaply[1::][::]
    return forces


# BOUNDARY CONDITIONS


def get_constrain(modelinfo, blist):
    """_summary_

    Arguments:
        modelinfo -- _description_
        blist -- _description_

    Returns:
        _description_
    """
    boncdnodeaply = np.zeros((1, 2))
    node_list_bc = np.array([])
    for i in range(len(blist)):
        if blist[i][0] == "fixed":
            # ----- SEEKERS IN LOC -----
            if blist[i][2] == "edgex":
                from myfempy.felib.physics.getnode import search_edgex

                edge_coordX = float(blist[i, 3])
                if float(blist[i, 4]) == 999:
                    coord_bc = (
                        modelinfo["coord"][
                            np.where(modelinfo["coord"][:, 3] == float(blist[i, 5])), :
                        ]
                    )[0]
                elif float(blist[i, 5]) == 999:
                    coord_bc = (
                        modelinfo["coord"][
                            np.where(modelinfo["coord"][:, 2] == float(blist[i, 4])), :
                        ]
                    )[0]
                node_list_bc = search_edgex(edge_coordX, coord_bc, 2e-3)
            elif blist[i][2] == "edgey":
                from myfempy.felib.physics.getnode import search_edgey

                edge_coordY = float(blist[i, 4])
                if float(blist[i, 3]) == 999:
                    coord_bc = (
                        modelinfo["coord"][
                            np.where(modelinfo["coord"][:, 3] == float(blist[i, 5])), :
                        ]
                    )[0]
                elif float(blist[i, 5]) == 999:
                    coord_bc = (
                        modelinfo["coord"][
                            np.where(modelinfo["coord"][:, 1] == float(blist[i, 3])), :
                        ]
                    )[0]
                node_list_bc = search_edgey(edge_coordY, coord_bc, 2e-3)
            elif blist[i][2] == "edgez":
                from myfempy.felib.physics.getnode import search_edgez

                edge_coordZ = float(blist[i, 5])
                if float(blist[i, 3]) == 999:
                    coord_bc = (
                        modelinfo["coord"][
                            np.where(modelinfo["coord"][:, 2] == float(blist[i, 4])), :
                        ]
                    )[0]
                elif float(blist[i, 4]) == 999:
                    coord_bc = (
                        modelinfo["coord"][
                            np.where(modelinfo["coord"][:, 1] == float(blist[i, 3])), :
                        ]
                    )[0]
                node_list_bc = search_edgez(edge_coordZ, coord_bc, 2e-3)
            elif blist[i][2] == "surfxy":
                from myfempy.felib.physics.getnode import search_surfxy

                orthg_coordZ = float(blist[i, 5])
                node_list_bc = search_surfxy(orthg_coordZ, modelinfo["coord"], 2e-3)
            elif blist[i][2] == "surfyz":
                from myfempy.felib.physics.getnode import search_surfyz

                orthg_coordX = float(blist[i, 3])
                node_list_bc = search_surfyz(orthg_coordX, modelinfo["coord"], 2e-3)
            elif blist[i][2] == "surfzx":
                from myfempy.felib.physics.getnode import search_surfzx

                orthg_coordY = float(blist[i, 4])
                node_list_bc = search_surfzx(orthg_coordY, modelinfo["coord"], 2e-3)
            elif blist[i][2] == "node":
                from myfempy.felib.physics.getnode import search_nodexyz

                node_coordX = float(blist[i, 3])
                node_coordY = float(blist[i, 4])
                node_coordZ = float(blist[i, 5])
                node_list_bc = search_nodexyz(
                    node_coordX, node_coordY, node_coordZ, modelinfo["coord"], 2e-3
                )
            # ----- SEEKERS IN TAG -----    
            elif blist[i][2] == "point":
                node_list_bc = modelinfo["regions"][0][1][int(blist[i][6]) - 1][1][:]
            elif blist[i][2] == "edge":
                node_list_bc = modelinfo["regions"][1][1][int(blist[i][6]) - 1][1][:]
            elif blist[i][2] == "surf":
                node_list_bc = modelinfo["regions"][2][1][int(blist[i][6]) - 1][1][:]
            else:
                print("input erro: bc_opt don't defined")
            for j in range(len(node_list_bc)):
                if blist[i][1] == "ux":
                    bcdef = np.array([[1, int(node_list_bc[j])]])
                    boncdnodeaply = np.append(boncdnodeaply, bcdef, axis=0)
                elif blist[i][1] == "uy":
                    bcdef = np.array([[2, int(node_list_bc[j])]])
                    boncdnodeaply = np.append(boncdnodeaply, bcdef, axis=0)
                elif blist[i][1] == "uz":
                    bcdef = np.array([[3, int(node_list_bc[j])]])
                    boncdnodeaply = np.append(boncdnodeaply, bcdef, axis=0)
                elif blist[i][1] == "rx":
                    bcdef = np.array([[4, int(node_list_bc[j])]])
                    boncdnodeaply = np.append(boncdnodeaply, bcdef, axis=0)
                elif blist[i][1] == "ry":
                    bcdef = np.array([[5, int(node_list_bc[j])]])
                    boncdnodeaply = np.append(boncdnodeaply, bcdef, axis=0)
                elif blist[i][1] == "rz":
                    bcdef = np.array([[6, int(node_list_bc[j])]])
                    boncdnodeaply = np.append(boncdnodeaply, bcdef, axis=0)
                elif blist[i][1] == "all":
                    bcdef = np.array([[0, int(node_list_bc[j])]])
                    boncdnodeaply = np.append(boncdnodeaply, bcdef, axis=0)
                else:
                    print("input erro: bc_opt_dir don't defined")
        else:
            print("input erro: bc_typ don't defined")
    constrains = boncdnodeaply[1::][::]
    return constrains
