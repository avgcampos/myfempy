#!/usr/bin/env python
__doc__ = """
Filters I/O
"""
import csv

import numpy as np


def write2log(log_file, log_data, modelinfo, solstatus):
    with open(log_file, "w") as file_object:
        file_object.write(
            "===============================================================================\n"
        )
        file_object.write(
            "                                   MYFEMPY                                     \n"
        )
        file_object.write(
            "===============================================================================\n"
        )
        if "get" in log_data.keys():
            if "nelem" in log_data["get"].keys():
                file_object.write("\n")
                file_object.write("GET NELEM " + str(len(modelinfo["inci"])) + "\n")
            if "nnode" in log_data["get"].keys():
                file_object.write("\n")
                file_object.write("GET NNODE " + str(len(modelinfo["coord"])) + "\n")
            if "inci" in log_data["get"].keys():
                file_object.write("\n")
                file_object.write("LIST OF ELEMENTS\n")
                file_object.write(
                    "{0:<7}{1:<7}{2:<7}{3:<7}{4:<25}\n".format(
                        "ELEM", "KEY", "MAT", "GEO", "NODES"
                    )
                )
                for row in range(len(modelinfo["inci"])):
                    node_list = modelinfo["inci"][row][4:].astype(
                        int
                    )  # node_list[np.nonzero(node_list)].tolist()
                    node_list = node_list[np.nonzero(node_list)].tolist()
                    node_list = " ".join(str(x) for x in node_list)
                    file_object.write(
                        "{0:<7}{1:<7}{2:<7}{3:<7}{4:<25}\n".format(
                            modelinfo["inci"][row][0].astype(int),
                            modelinfo["inci"][row][1].astype(int),
                            modelinfo["inci"][row][2].astype(int),
                            modelinfo["inci"][row][3].astype(int),
                            node_list,
                        )
                    )
            if "coord" in log_data["get"].keys():
                file_object.write("\n")
                file_object.write("LIST OF NODES COORDINATE\n")
                file_object.write(
                    "{0:<7}{1:<10}{2:<10}{3:<10}\n".format("NODE", "X", "Y", "Z")
                )
                for row in range(len(modelinfo["coord"])):
                    file_object.write(
                        "{0:<7}{1:<10}{2:<10}{3:<10}\n".format(
                            modelinfo["coord"][row][0].astype(int),
                            modelinfo["coord"][row][1],
                            modelinfo["coord"][row][2],
                            modelinfo["coord"][row][3],
                        )
                    )
            if "boundcond_list" in log_data["get"].keys():
                file_object.write("\n")
                file_object.write("LIST OF CONSTRAINTS\n")
                file_object.write("{0:<7}{1:<10}\n".format("BC", "NODE"))
                for row in range(len(modelinfo["constrains"])):
                    bc_type = modelinfo["constrains"][row][0]
                    if bc_type == 0:
                        file_object.write(
                            "{0:<7}{1:<10}\n".format(
                                "FULL", modelinfo["constrains"][row][1].astype(int)
                            )
                        )
                    elif bc_type == 1:
                        file_object.write(
                            "{0:<7}{1:<10}\n".format(
                                "UX", modelinfo["constrains"][row][1].astype(int)
                            )
                        )
                    elif bc_type == 2:
                        file_object.write(
                            "{0:<7}{1:<10}\n".format(
                                "UY", modelinfo["constrains"][row][1].astype(int)
                            )
                        )
                    elif bc_type == 3:
                        file_object.write(
                            "{0:<7}{1:<10}\n".format(
                                "UZ", modelinfo["constrains"][row][1].astype(int)
                            )
                        )
                    elif bc_type == 4:
                        file_object.write(
                            "{0:<7}{1:<10}\n".format(
                                "RX", modelinfo["constrains"][row][1].astype(int)
                            )
                        )
                    elif bc_type == 5:
                        file_object.write(
                            "{0:<7}{1:<10}\n".format(
                                "RY", modelinfo["constrains"][row][1].astype(int)
                            )
                        )
                    elif bc_type == 6:
                        file_object.write(
                            "{0:<7}{1:<10}\n".format(
                                "RZ", modelinfo["constrains"][row][1].astype(int)
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
                for row in range(len(modelinfo["forces"])):
                    fc_type = modelinfo["forces"][row][1]
                    if fc_type == 1:
                        file_object.write(
                            "{0:<7}{1:<10}{2:<25}{3:<10}\n".format(
                                modelinfo["forces"][row][0].astype(int),
                                "FX",
                                modelinfo["forces"][row][2],
                                modelinfo["forces"][row][3].astype(int),
                            )
                        )
                    elif fc_type == 2:
                        file_object.write(
                            "{0:<7}{1:<10}{2:<25}{3:<10}\n".format(
                                modelinfo["forces"][row][0].astype(int),
                                "FY",
                                modelinfo["forces"][row][2],
                                modelinfo["forces"][row][3].astype(int),
                            )
                        )
                    elif fc_type == 3:
                        file_object.write(
                            "{0:<7}{1:<10}{2:<25}{3:<10}\n".format(
                                modelinfo["forces"][row][0].astype(int),
                                "FZ",
                                modelinfo["forces"][row][2],
                                modelinfo["forces"][row][3].astype(int),
                            )
                        )
                    elif fc_type == 4:
                        file_object.write(
                            "{0:<7}{1:<10}{2:<25}{3:<10}\n".format(
                                modelinfo["forces"][row][0].astype(int),
                                "TX",
                                modelinfo["forces"][row][2],
                                modelinfo["forces"][row][3].astype(int),
                            )
                        )
                    elif fc_type == 5:
                        file_object.write(
                            "{0:<7}{1:<10}{2:<25}{3:<10}\n".format(
                                modelinfo["forces"][row][0].astype(int),
                                "TY",
                                modelinfo["forces"][row][2],
                                modelinfo["forces"][row][3].astype(int),
                            )
                        )
                    elif fc_type == 6:
                        file_object.write(
                            "{0:<7}{1:<10}{2:<25}{3:<10}\n".format(
                                modelinfo["forces"][row][0].astype(int),
                                "TZ",
                                modelinfo["forces"][row][2],
                                modelinfo["forces"][row][3].astype(int),
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
                for row in range(len(modelinfo["tabmat"])):
                    file_object.write(
                        "{0:<7}{1:<10}{2:<10}{3:<10}{4:<10}{5:<10}{6:<10}{7:<10}{8:<10}{9:<10}{10:<10}{11:<10}\n".format(
                            str(row + 1),
                            modelinfo["tabmat"][row][0],
                            modelinfo["tabmat"][row][1],
                            modelinfo["tabmat"][row][2],
                            modelinfo["tabmat"][row][3],
                            modelinfo["tabmat"][row][4],
                            modelinfo["tabmat"][row][5],
                            modelinfo["tabmat"][row][6],
                            modelinfo["tabmat"][row][7],
                            modelinfo["tabmat"][row][8],
                            modelinfo["tabmat"][row][9].astype(int),
                            modelinfo["tabmat"][row][10].astype(int),
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
                for row in range(len(modelinfo["tabgeo"])):
                    file_object.write(
                        "{0:<7}{1:<10}{2:<10}{3:<10}{4:<10}{5:<10}{6:<10}{7:<10}{8:<10}{9:<10}{10:<10}\n".format(
                            str(row + 1),
                            modelinfo["tabgeo"][row][0],
                            modelinfo["tabgeo"][row][1],
                            modelinfo["tabgeo"][row][2],
                            modelinfo["tabgeo"][row][3],
                            modelinfo["tabgeo"][row][4],
                            modelinfo["tabgeo"][row][5],
                            modelinfo["tabgeo"][row][6],
                            modelinfo["tabgeo"][row][7],
                            modelinfo["tabgeo"][row][8],
                            modelinfo["tabgeo"][row][9].astype(int),
                        )
                    )
        # =============================================================================================================
        if "log" in log_data.keys():
            file_object.write("\n")
            file_object.write(
                "+------------------S O L V E R   S T A T U S------------------+\n"
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
                    "NUMBER OF EQUATION ", str(modelinfo["fulldofs"])
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
