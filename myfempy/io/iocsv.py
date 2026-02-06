#!/usr/bin/env python
import csv

import numpy as np


def write2log(Model, Physic, log_data, solstatus, log_file):
    with open(log_file, "w") as file_object:
        file_object.write(
            "===============================================================================\n"
        )
        file_object.write(
            "                        M Y F E M P Y   -   R E P O R T                        \n"
        )
        file_object.write(
            "===============================================================================\n"
        )
        if "get" in log_data.keys():
            if "nelem" in log_data["get"].keys():
                file_object.write("\n")
                file_object.write("GET NELEM " + str(len(Model.inci)) + "\n")
            if "nnode" in log_data["get"].keys():
                file_object.write("\n")
                file_object.write("GET NNODE " + str(len(Model.coord)) + "\n")
            if "inci" in log_data["get"].keys():
                file_object.write("\n")
                file_object.write("LIST OF ELEMENTS\n")
                file_object.write(
                    "{0:<7}{1:<7}{2:<7}{3:<7}{4:<25}\n".format(
                        "ELEM", "KEY", "MAT", "GEO", "NODES"
                    )
                )
                for row in range(len(Model.inci)):
                    node_list = Model.inci[row][4:].astype(
                        int
                    )  # node_list[np.nonzero(node_list)].tolist()
                    node_list = node_list[np.nonzero(node_list)].tolist()
                    node_list = " ".join(str(x) for x in node_list)
                    file_object.write(
                        "{0:<7}{1:<7}{2:<7}{3:<7}{4:<25}\n".format(
                            Model.inci[row][0].astype(int),
                            Model.inci[row][1].astype(int),
                            Model.inci[row][2].astype(int),
                            Model.inci[row][3].astype(int),
                            node_list,
                        )
                    )
            if "coord" in log_data["get"].keys():
                file_object.write("\n")
                file_object.write("LIST OF NODES COORDINATE\n")
                file_object.write(
                    "{0:<7}{1:<10}{2:<10}{3:<10}\n".format("NODE", "X", "Y", "Z")
                )
                for row in range(len(Model.coord)):
                    file_object.write(
                        "{0:<7}{1:<20}{2:<20}{3:<20}\n".format(
                            Model.coord[row][0].astype(int),
                            Model.coord[row][1],
                            Model.coord[row][2],
                            Model.coord[row][3],
                        )
                    )
            if "boundcond_list" in log_data["get"].keys():
                file_object.write("\n")
                file_object.write("LIST OF CONSTRAINTS\n")
                file_object.write("{0:<7}{1:<10}\n".format("DOF", "NODE"))
                for row in range(len(Physic.constrains)):
                    bc_type = Physic.constrains[row][0]
                    if bc_type == 0:
                        file_object.write(
                            "{0:<7}{1:<10}\n".format(
                                "FULL", Physic.constrains[row][1].astype(int)
                            )
                        )
                    elif bc_type == 1:
                        file_object.write(
                            "{0:<7}{1:<10}\n".format(
                                "UX", Physic.constrains[row][1].astype(int)
                            )
                        )
                    elif bc_type == 2:
                        file_object.write(
                            "{0:<7}{1:<10}\n".format(
                                "UY", Physic.constrains[row][1].astype(int)
                            )
                        )
                    elif bc_type == 3:
                        file_object.write(
                            "{0:<7}{1:<10}\n".format(
                                "UZ", Physic.constrains[row][1].astype(int)
                            )
                        )
                    elif bc_type == 4:
                        file_object.write(
                            "{0:<7}{1:<10}\n".format(
                                "RX", Physic.constrains[row][1].astype(int)
                            )
                        )
                    elif bc_type == 5:
                        file_object.write(
                            "{0:<7}{1:<10}\n".format(
                                "RY", Physic.constrains[row][1].astype(int)
                            )
                        )
                    elif bc_type == 6:
                        file_object.write(
                            "{0:<7}{1:<10}\n".format(
                                "RZ", Physic.constrains[row][1].astype(int)
                            )
                        )
                    else:
                        pass
            if "forces_list" in log_data["get"].keys():
                file_object.write("\n")
                file_object.write("LIST OF FORCES\n")
                file_object.write(
                    "{0:<7}{1:<10}{2:<25}{3:<10}\n".format(
                        "NODE", "TYPE", "VALUE", "STEP"
                    )
                )
                for row in range(len(Physic.forces)):
                    fc_type = Physic.forces[row][1]
                    if fc_type == 1:
                        file_object.write(
                            "{0:<7}{1:<10}{2:<25}{3:<10}\n".format(
                               Physic.forces[row][0].astype(int),
                                "FX",
                                Physic.forces[row][2],
                                Physic.forces[row][3].astype(int),
                            )
                        )
                    elif fc_type == 2:
                        file_object.write(
                            "{0:<7}{1:<10}{2:<25}{3:<10}\n".format(
                                Physic.forces[row][0].astype(int),
                                "FY",
                                Physic.forces[row][2],
                                Physic.forces[row][3].astype(int),
                            )
                        )
                    elif fc_type == 3:
                        file_object.write(
                            "{0:<7}{1:<10}{2:<25}{3:<10}\n".format(
                                Physic.forces[row][0].astype(int),
                                "FZ",
                                Physic.forces[row][2],
                                Physic.forces[row][3].astype(int),
                            )
                        )
                    elif fc_type == 4:
                        file_object.write(
                            "{0:<7}{1:<10}{2:<25}{3:<10}\n".format(
                                Physic.forces[row][0].astype(int),
                                "TX",
                                Physic.forces[row][2],
                                Physic.forces[row][3].astype(int),
                            )
                        )
                    elif fc_type == 5:
                        file_object.write(
                            "{0:<7}{1:<10}{2:<25}{3:<10}\n".format(
                                Physic.forces[row][0].astype(int),
                                "TY",
                                Physic.forces[row][2],
                                Physic.forces[row][3].astype(int),
                            )
                        )
                    elif fc_type == 6:
                        file_object.write(
                            "{0:<7}{1:<10}{2:<25}{3:<10}\n".format(
                                Physic.forces[row][0].astype(int),
                                "TZ",
                                Physic.forces[row][2],
                                Physic.forces[row][3].astype(int),
                            )
                        )
                    else:
                        pass
            if "tabmat" in log_data["get"].keys():
                file_object.write("\n")
                file_object.write("LIST OF MATERIAL PROPERTY\n")
                file_object.write(
                    "{0:<7}{1:<10}{2:<10}{3:<10}{4:<10}{5:<10}{6:<10}{7:<10}{8:<10}{9:<10}{10:<10}{11:<10}{12:<10}{13:<10}{14:<10}{15:<10}{16:<10}\n".format(
                        "EXX",
                        "VXY",
                        "GXY",
                        "EYY",
                        "VYZ",
                        "GYZ",
                        "EZZ",
                        "VZX",
                        "GZX",
                        "RHO",
                        "KXX",
                        "KYY",
                        "KZZ",
                        "CTE",
                        "VIS",
                        "STIF",
                        "DAMP",
                    )
                )
                for row in range(len(Model.tabmat)):
                    file_object.write(
                        "{0:<7}{1:<10}{2:<10}{3:<10}{4:<10}{5:<10}{6:<10}{7:<10}{8:<10}{9:<10}{10:<10}{11:<10}\n".format(
                            str(row + 1),
                            Model.tabmat[row][0],
                            Model.tabmat[row][1],
                            Model.tabmat[row][2],
                            Model.tabmat[row][3],
                            Model.tabmat[row][4],
                            Model.tabmat[row][5],
                            Model.tabmat[row][6],
                            Model.tabmat[row][7],
                            Model.tabmat[row][8],
                            Model.tabmat[row][9].astype(int),
                            Model.tabmat[row][10].astype(int),
                        )
                    )
            if "tabgeo" in log_data["get"].keys():
                file_object.write("\n")
                file_object.write("LIST OF MATERIAL PROPERTY\n")
                file_object.write(
                    "{0:<7}{1:<10}{2:<10}{3:<10}{4:<10}{5:<10}{6:<10}{7:<10}{8:<10}{9:<10}}\n".format(
                        "AREACS",
                        "INERYY",
                        "INERZZ",
                        "INERXX",
                        "THICKN",
                        "b",
                        "h",
                        "t",
                        "d",
                        "ID",
                    )
                )
                for row in range(len(Model.tabgeo)):
                    file_object.write(
                        "{0:<7}{1:<10}{2:<10}{3:<10}{4:<10}{5:<10}{6:<10}{7:<10}{8:<10}{9:<10}{10:<10}\n".format(
                            str(row + 1),
                            Model.tabgeo[row][0],
                            Model.tabgeo[row][1],
                            Model.tabgeo[row][2],
                            Model.tabgeo[row][3],
                            Model.tabgeo[row][4],
                            Model.tabgeo[row][5],
                            Model.tabgeo[row][6],
                            Model.tabgeo[row][7],
                            Model.tabgeo[row][8],
                            Model.tabgeo[row][9].astype(int),
                        )
                    )
        # =============================================================================================================
        if "log" in log_data.keys():
            file_object.write("\n")
            file_object.write(
                "+---------------------------- S O L V E R   L O G ----------------------------+\n"
            )
            file_object.write(
                "{0:<30} : {1:<10} SEC\n".format(
                    "ASSEMBLY FULL TIME SPEND ",
                    str(solstatus["solverstatus"]["timeasb"]),
                )
            )
            file_object.write(
                "{0:<30} : {1:<10} SEC\n".format(
                    "SOLVE FULL TIME SPEND ", str(solstatus["solverstatus"]["timesim"])
                )
            )
            file_object.write(
                "{0:<30} : {1:<10} DOF\n".format(
                    "NUMBER OF EQUATION ", str(Model.modelinfo["fulldofs"])
                )
            )
            file_object.write(
                "{0:<30} : {1:<10} MB\n".format(
                    "STIFFNESS SIZE ", str(solstatus["solverstatus"]["memorysize"])
                )
            )
            file_object.write(
                "{0:<30} : {1:<10} SET\n".format(
                    "TYPE ASSEMBLER ", str(solstatus["solverstatus"]["typeasmb"])
                )
            )
            file_object.write(
                "{0:<30} : {1:<10} INT\n".format(
                    "SOLVER CORE ", str(solstatus["solverstatus"]["ncpu"])
                )
            )


def writer2csv(csv2write_file, csv2write_data, label):
    with open(csv2write_file, "w") as file_object:
        writer = csv.writer(file_object)
        writer.writerow([label[0], label[1]])
        for tt in range(len(csv2write_data[0])):
            writer.writerow([csv2write_data[0][tt], csv2write_data[1][tt]])


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
