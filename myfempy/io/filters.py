#!/usr/bin/env python
__doc__ = """
Filters I/O
"""
import numpy as np


def input_filter(path_user, file_name):
    pointlist = np.zeros((1, 4))
    linelist = np.zeros((1, 3))
    propgeonewlist = np.zeros((1, 6))
    propgeobiblist = np.zeros((1, 6))
    compmaterial = np.zeros((1, 4))
    propmatlist = np.zeros((1, 8))
    solutconfig = np.zeros((1, 3))
    meshconfig = np.zeros((1, 5))
    forcelist = np.zeros((1, 9))
    boncdlist = np.zeros((1, 7))
    solvercfg = np.zeros((1, 7))
    outputcfg = np.zeros((1, 4))
    file_input = str(path_user + "/" + file_name)
    with open(file_input, "r") as file_object:
        line = file_object.readline()
        line = file_object.readline()
        line = file_object.readline()
        line = file_object.readline()
        line = file_object.readline()
        line = file_object.readline()
        line = file_object.readline()
        if line == "#SOLVER_CFG\n":
            line = file_object.readline()
            while line != "\n":
                lineaux = line.split()
                linearray = np.array(
                    [
                        lineaux[0],
                        lineaux[2],
                        lineaux[3],
                        lineaux[4],
                        lineaux[5],
                        lineaux[6],
                        lineaux[7],
                    ]
                )
                solvercfg = np.append(solvercfg, [linearray], axis=0)
                line = file_object.readline()
        line = file_object.readline()
        if line == "#IMPORT_GEO\n":
            line = file_object.readline()
        line = file_object.readline()
        if line == "#NEW_GEO\n":
            line = file_object.readline()
            if line == "##POINT_COORD\n":
                line = file_object.readline()
                while line != "\n":
                    lineaux = line.split()
                    linearray = np.array(
                        [lineaux[0], lineaux[2], lineaux[3], lineaux[4]]
                    )
                    pointlist = np.append(pointlist, [linearray], axis=0)
                    line = file_object.readline()
            line = file_object.readline()
            if line == "##LINE_CONEC\n":
                line = file_object.readline()
                while line != "\n":
                    lineaux = line.split()
                    linearray = np.array([lineaux[0], lineaux[2], lineaux[3]])
                    linelist = np.append(linelist, [linearray], axis=0)
                    line = file_object.readline()
            line = file_object.readline()
            planelist = []
            if line == "##PLANE_CONEC\n":
                line = file_object.readline()
                while line != "\n":
                    lineaux = line.split()
                    planelist.append(lineaux[2:])
                    line = file_object.readline()
        line = file_object.readline()
        if line == "#PROP_GEOMETRIA\n":
            line = file_object.readline()
            if line == "##GEO_INP\n":
                line = file_object.readline()
                while line != "\n":
                    lineaux = line.split()
                    linearray = np.array(
                        [
                            lineaux[0],
                            lineaux[2],
                            lineaux[3],
                            lineaux[4],
                            lineaux[5],
                            lineaux[6],
                        ]
                    )
                    propgeonewlist = np.append(propgeonewlist, [linearray], axis=0)
                    line = file_object.readline()
            line = file_object.readline()
            if line == "##GEO_BIB\n":
                line = file_object.readline()
                while line != "\n":
                    lineaux = line.split()
                    linearray = np.array(
                        [
                            lineaux[0],
                            lineaux[2],
                            lineaux[3],
                            lineaux[4],
                            lineaux[5],
                            lineaux[6],
                        ]
                    )
                    propgeobiblist = np.append(propgeobiblist, [linearray], axis=0)
                    line = file_object.readline()
        line = file_object.readline()
        if line == "#PROP_MATERIAL\n":
            line = file_object.readline()
            while line != "\n":
                lineaux = line.split()
                linearray = np.array([lineaux[0], lineaux[2], lineaux[3], lineaux[4]])
                compmaterial = np.append(compmaterial, [linearray], axis=0)
                line = file_object.readline()
        line = file_object.readline()
        if line == "##MAT_CFG\n":
            line = file_object.readline()
            while line != "\n":
                lineaux = line.split()
                linearray = np.array(
                    [
                        lineaux[0],
                        lineaux[2],
                        lineaux[3],
                        lineaux[4],
                        lineaux[5],
                        lineaux[6],
                        lineaux[7],
                        lineaux[8],
                    ]
                )
                propmatlist = np.append(propmatlist, [linearray], axis=0)
                line = file_object.readline()
        line = file_object.readline()
        if line == "#MODEL_CFG\n":
            line = file_object.readline()
            while line != "\n":
                lineaux = line.split()
                linearray = np.array([lineaux[0], lineaux[2], lineaux[3]])
                solutconfig = np.append(solutconfig, [linearray], axis=0)
                line = file_object.readline()
        line = file_object.readline()
        if line == "#MESH_CFG\n":
            line = file_object.readline()
            while line != "\n":
                lineaux = line.split()
                linearray = np.array(
                    [lineaux[0], lineaux[2], lineaux[3], lineaux[4], lineaux[5]]
                )
                meshconfig = np.append(meshconfig, [linearray], axis=0)
                line = file_object.readline()
        line = file_object.readline()
        if line == "#SOLUTION_CFG\n":
            line = file_object.readline()
            if line == "##LOADS\n":
                line = file_object.readline()
                while line != "\n":
                    lineaux = line.split()
                    linearray = np.array(
                        [
                            lineaux[0],
                            lineaux[2],
                            lineaux[3],
                            lineaux[4],
                            lineaux[5],
                            lineaux[6],
                            lineaux[7],
                            lineaux[8],
                            lineaux[9],
                        ]
                    )
                    forcelist = np.append(forcelist, [linearray], axis=0)
                    line = file_object.readline()
            line = file_object.readline()
            if line == "##BOND_COND\n":
                line = file_object.readline()
                while line != "\n":
                    lineaux = line.split()
                    linearray = np.array(
                        [
                            lineaux[0],
                            lineaux[2],
                            lineaux[3],
                            lineaux[4],
                            lineaux[5],
                            lineaux[6],
                            lineaux[7],
                        ]
                    )
                    boncdlist = np.append(boncdlist, [linearray], axis=0)
                    line = file_object.readline()
        line = file_object.readline()
        if line == "#VIEWSOLUTION_CFG\n":
            line = file_object.readline()
            while line != "\n":
                lineaux = line.split()
                linearray = np.array([lineaux[0], lineaux[2], lineaux[3], lineaux[4]])
                outputcfg = np.append(outputcfg, [linearray], axis=0)
                line = file_object.readline()
        line = file_object.readline()
        if line == "#GRAPH_OUT\n":
            line = file_object.readline()
            while line != "\n":
                lineaux = line.split()
                graphout = np.zeros((1, 1))
                graphout = np.append(graphout, np.array([[lineaux[0]]]), axis=1)
                graphout = np.append(graphout, np.array([lineaux[2:]]), axis=1)
                line = file_object.readline()
        line = file_object.readline()
        if line == "#FILE_OUT\n":
            line = file_object.readline()
            while line != "\n":
                lineaux = line.split()
                fileout = np.zeros((1, 1))
                fileout = np.append(fileout, np.array([[lineaux[0]]]), axis=1)
                fileout = np.append(fileout, np.array([lineaux[2:]]), axis=1)
                line = file_object.readline()
    pointlist = pointlist[1::][::]
    linelist = linelist[1::][::]
    propgeonewlist = propgeonewlist[1::][::]
    propgeobiblist = propgeobiblist[1::][::]
    propmatlist = propmatlist[1::][::]
    compmaterial = compmaterial[1::][::]
    solutconfig = solutconfig[1::][::]
    meshconfig = meshconfig[1::][::]
    forcelist = forcelist[1::][::]
    boncdlist = boncdlist[1::][::]
    solvercfg = solvercfg[1::][::]
    outputcfg = outputcfg[1::][::]
    graphout = graphout[0, 1:]
    fileout = fileout[0, 1:]
    return (
        solutconfig,
        solvercfg,
        pointlist,
        linelist,
        planelist,
        propgeonewlist,
        propgeobiblist,
        compmaterial,
        propmatlist,
        meshconfig,
        forcelist,
        boncdlist,
        solvercfg,
        outputcfg,
        graphout,
        fileout,
    )


def gen_usrlog(
    usrlog_file,
    solutconfig,
    meshconfig,
    compmaterial,
    forcelist,
    boncdlist,
    solvercfg,
    outputcfg,
    graphout,
    fileout,
):
    with open(usrlog_file, "w") as file_object:
        for i in range(len(solvercfg)):
            file_object.write("solver_typ = " + "'" + solvercfg[i, 1] + "'" + "\n")
            file_object.write("solver_def = " + "'" + solvercfg[i, 2] + "'" + "\n")
            file_object.write("solver_opt = " + "'" + solvercfg[i, 3] + "'" + "\n")
            file_object.write("solver_start = " + "'" + solvercfg[i, 4] + "'" + "\n")
            file_object.write("solver_end = " + "'" + solvercfg[i, 5] + "'" + "\n")
            file_object.write("solver_step = " + "'" + solvercfg[i, 6] + "'" + "\n")
        for i in range(len(solutconfig)):
            file_object.write("mod_typ = " + "'" + solutconfig[i, 1] + "'" + "\n")
            file_object.write("mod_opt = " + "'" + solutconfig[i, 2] + "'" + "\n")
        for i in range(len(meshconfig)):
            file_object.write("mesh_typ = " + "'" + meshconfig[i, 1] + "'" + "\n")
        for i in range(len(compmaterial)):
            file_object.write("mat_typ = " + "'" + compmaterial[i, 1] + "'" + "\n")
            file_object.write("mat_opt = " + "'" + compmaterial[i, 2] + "'" + "\n")
            file_object.write("mat_def = " + "'" + compmaterial[i, 3] + "'" + "\n")
        for j in range(len(np.unique(forcelist[:, 8]))):
            forceliststep = forcelist[np.where(forcelist[:, 8] == str(j + 1)), :][0]
            for i in range(len(forceliststep)):
                file_object.write(
                    "force_typ_" + str(i) + " "
                    "= " + "'" + forceliststep[i, 1] + "'" + "\n"
                )
                file_object.write(
                    "force_opt_" + str(i) + " "
                    "= " + "'" + forceliststep[i, 4] + "'" + "\n"
                )
                file_object.write(
                    "force_opt_dir" + str(i) + " "
                    "= " + "'" + forceliststep[i, 2] + "'" + "\n"
                )
                file_object.write(
                    "force_step" + str(j) + " "
                    "= " + "'" + forceliststep[i, 8] + "'" + "\n"
                )
        for i in range(len(boncdlist)):
            file_object.write(
                "bc_typ_" + str(i) + " " "= " + "'" + boncdlist[i, 1] + "'" + "\n"
            )
            file_object.write(
                "bc_opt_" + str(i) + " " "= " + "'" + boncdlist[i, 3] + "'" + "\n"
            )
            file_object.write(
                "bc_opt_dir" + str(i) + " " "= " + "'" + boncdlist[i, 2] + "'" + "\n"
            )
        file_object.write("prev_sol = " + "'" + outputcfg[0, 1] + "'" + "\n")
        file_object.write("graph_out = " + "'" + outputcfg[0, 2] + "'" + "\n")
        file_object.write("save_fig = " + "'" + outputcfg[0, 3] + "'" + "\n")
        for i in range(len(graphout)):
            file_object.write(f"graph_out_{i} = " + "'" + graphout[i] + "'" + "\n")
        for i in range(len(fileout)):
            file_object.write(f"file_out_{i} = " + "'" + fileout[i] + "'" + "\n")


def import_mesh_line(file_imp):
    with open(file_imp, "r") as file_object:
        file_object.readline()
        NNOD = int(file_object.readline())
        coord = np.zeros((NNOD, 4))
        for ii in range(0, NNOD):
            line = file_object.readline()
            lineaux = line.split()
            contstr = lineaux[0:4]
            coord[ii, :] = [
                float(contstr[0]),
                float(contstr[1]),
                float(contstr[2]),
                float(contstr[3]),
            ]
        file_object.readline()
        file_object.readline()
        NELM = int(file_object.readline())
        conec_elm = np.zeros((NELM, 9))
        prop_elm = np.zeros((NELM, 3))
        contelm = 0
        for kk in range(0, NELM):
            line = file_object.readline()
            lineaux = line.split()
            type_elm = int(lineaux[1])
            if type_elm == 1:
                contelm += 1
                prop_elm[contelm - 1, :] = [
                    type_elm,
                    float(lineaux[3]),
                    float(lineaux[3]),
                ]
                conec_elm[contelm - 1, :] = [
                    contelm,
                    float(lineaux[5]),
                    float(lineaux[6]),
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                ]
    conec_elm = conec_elm[:contelm, :]
    prop_elm = prop_elm[:contelm, :]
    inci = np.concatenate(
        (conec_elm[:, 0][:, np.newaxis], prop_elm, conec_elm[:, 1:]), axis=1
    )
    return coord, inci, type_elm


def import_mesh_tria(file_imp):
    with open(file_imp, "r") as file_object:
        file_object.readline()
        NNOD = int(file_object.readline())
        coord = np.zeros((NNOD, 4))
        for ii in range(0, NNOD):
            line = file_object.readline()
            lineaux = line.split()
            contstr = lineaux[0:4]
            coord[ii, :] = [
                float(contstr[0]),
                float(contstr[1]),
                float(contstr[2]),
                float(contstr[3]),
            ]
        file_object.readline()
        file_object.readline()
        NELM = int(file_object.readline())
        conec_elm = np.zeros((NELM, 9))
        prop_elm = np.zeros((NELM, 3))
        contelm = 0
        for kk in range(0, NELM):
            line = file_object.readline()
            lineaux = line.split()
            type_elm = int(lineaux[1])
            if type_elm == 2:
                contelm += 1
                prop_elm[contelm - 1, :] = [
                    type_elm,
                    float(lineaux[3]),
                    float(lineaux[3]),
                ]
                conec_elm[contelm - 1, :] = [
                    contelm,
                    float(lineaux[5]),
                    float(lineaux[6]),
                    float(lineaux[7]),
                    0,
                    0,
                    0,
                    0,
                    0,
                ]
    conec_elm = conec_elm[:contelm, :]
    prop_elm = prop_elm[:contelm, :]
    inci = np.concatenate(
        (conec_elm[:, 0][:, np.newaxis], prop_elm, conec_elm[:, 1:]), axis=1
    )
    return coord, inci, type_elm


def import_mesh_quad(file_imp):
    with open(file_imp, "r") as file_object:
        file_object.readline()
        NNOD = int(file_object.readline())
        coord = np.zeros((NNOD, 4))
        for ii in range(0, NNOD):
            line = file_object.readline()
            lineaux = line.split()
            contstr = lineaux[0:4]
            coord[ii, :] = [
                float(contstr[0]),
                float(contstr[1]),
                float(contstr[2]),
                float(contstr[3]),
            ]
        file_object.readline()
        file_object.readline()
        NELM = int(file_object.readline())
        conec_elm = np.zeros((NELM, 9))
        prop_elm = np.zeros((NELM, 3))
        contelm = 0
        for kk in range(0, NELM):
            line = file_object.readline()
            lineaux = line.split()
            type_elm = int(lineaux[1])
            if type_elm == 3:
                contelm += 1
                prop_elm[contelm - 1, :] = [
                    type_elm,
                    float(lineaux[3]),
                    float(lineaux[3]),
                ]
                conec_elm[contelm - 1, :] = [
                    contelm,
                    float(lineaux[5]),
                    float(lineaux[6]),
                    float(lineaux[7]),
                    float(lineaux[8]),
                    0,
                    0,
                    0,
                    0,
                ]
    conec_elm = conec_elm[:contelm, :]
    prop_elm = prop_elm[:contelm, :]
    inci = np.concatenate(
        (conec_elm[:, 0][:, np.newaxis], prop_elm, conec_elm[:, 1:]), axis=1
    )
    return coord, inci, type_elm


def import_mesh_hexa(file_imp):
    with open(file_imp, "r") as file_object:
        file_object.readline()
        NNOD = int(file_object.readline())
        coord = np.zeros((NNOD, 4))
        for ii in range(0, NNOD):
            line = file_object.readline()
            lineaux = line.split()
            contstr = lineaux[0:4]
            coord[ii, :] = [
                float(contstr[0]),
                float(contstr[1]),
                float(contstr[2]),
                float(contstr[3]),
            ]
        file_object.readline()
        file_object.readline()
        NELM = int(file_object.readline())
        conec_elm = np.zeros((NELM, 9))
        prop_elm = np.zeros((NELM, 3))
        contelm = 0
        for kk in range(0, NELM):
            line = file_object.readline()
            lineaux = line.split()
            type_elm = int(lineaux[1])
            if type_elm == 5:
                contelm += 1
                prop_elm[contelm - 1, :] = [
                    type_elm,
                    float(lineaux[3]),
                    float(lineaux[3]),
                ]
                conec_elm[contelm - 1, :] = [
                    contelm,
                    float(lineaux[5]),
                    float(lineaux[6]),
                    float(lineaux[7]),
                    float(lineaux[8]),
                    float(lineaux[9]),
                    float(lineaux[10]),
                    float(lineaux[11]),
                    float(lineaux[12]),
                ]
    conec_elm = conec_elm[:contelm, :]
    prop_elm = prop_elm[:contelm, :]
    inci = np.concatenate(
        (conec_elm[:, 0][:, np.newaxis], prop_elm, conec_elm[:, 1:]), axis=1
    )
    return coord, inci, type_elm
