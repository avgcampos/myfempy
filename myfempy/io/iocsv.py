#!/usr/bin/env python
import csv

import numpy as np


__docformat__ = "google"

__doc__ = """

==========================================================================
                            __                                
         _ __ ___   _   _  / _|  ___  _ __ ___   _ __   _   _ 
        | '_ ` _ \ | | | || |_  / _ \| '_ ` _ \ | '_ \ | | | |
        | | | | | || |_| ||  _||  __/| | | | | || |_) || |_| |
        |_| |_| |_| \__, ||_|   \___||_| |_| |_|| .__/  \__, |
                    |___/                       |_|     |___/ 
        myfempy -- MultiphYsics Finite Element Module to PYthon    
                    COMPUTATIONAL ANALYSIS PROGRAM                   
        Copyright (C) 2022-2026 Antonio Vinicius Garcia Campos        
==========================================================================
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

def write2log(Model, Physic, log_data, solstatus, log_file):
    
    if "numpy_decimals" in log_data["get"].keys():
        np_decimals = log_data["get"]["numpy_decimals"]
    else:
        np_decimals = 4

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
                        "{0:<7}{1:<10}{2:<10}{3:<10}\n".format(
                            Model.coord[row][0].astype(int),
                            Model.coord[row][1].round(decimals=np_decimals),
                            Model.coord[row][2].round(decimals=np_decimals),
                            Model.coord[row][3].round(decimals=np_decimals),
                        )
                    )


            if "tabmat" in log_data["get"].keys():
                file_object.write("\n")
                file_object.write("LIST OF MATERIAL PROPERTY\n")

                for row in range(len(Model.tabmat)):
                    file_object.write("------------------------\n")
                    file_object.write("{0:<12}{1:<3}\n".format("MATERIAL", row + 1))
                    file_object.write("{0:<10}{1:<10}\n".format("PROPMAT", "VALUE"))
                    for key in Model.tabmat[row]:
                        value_propmat = Model.tabmat[row][key]
                        if value_propmat == "NULL":
                            pass
                        else:
                            file_object.write("{0:<10}{1:<10}\n".format(key, value_propmat))

            if "tabgeo" in log_data["get"].keys():
                file_object.write("\n")
                file_object.write("LIST OF GEOMETRY PROPERTY\n")

                for row in range(len(Model.tabgeo)):
                    file_object.write("------------------------\n")
                    file_object.write("{0:<12}{1:<3}\n".format("GEOMETRY", row + 1))
                    file_object.write("{0:<10}{1:<10}\n".format("PROPGEO", "VALUE"))
                    for key in Model.tabgeo[row]:
                        value_propmat = Model.tabgeo[row][key]
                        if value_propmat == "NULL":
                            pass
                        else:
                            file_object.write("{0:<10}{1:<10}\n".format(key, value_propmat))

            if "bc_list" in log_data["get"].keys():
                file_object.write("\n")
                file_object.write("LIST OF BOUNDARY CONDITIONS\n")
                file_object.write(
                    "{0:<7}{1:<10}{2:<16}{3:<10}\n".format(
                        "NODE", "DOF", "VALUE", "STEP"
                    )
                )
                model_dict = Model.modelinfo["dofs"]["d"]
                for row in range(len(Physic.constrains)):
                    bc_type = Physic.constrains[row][1]
                    if bc_type == 0:
                        file_object.write(
                            "{0:<7}{1:<10}{2:<16}{3:<10}\n".format(
                                Physic.constrains[row][0].astype(int),
                                "full",
                                Physic.constrains[row][2].round(decimals=np_decimals),
                                Physic.constrains[row][3].astype(int),
                            )
                        )
                    else:
                        bc_dof = next((k for k, v in model_dict.items() if v == bc_type), None)
                        file_object.write(
                            "{0:<7}{1:<10}{2:<16}{3:<10}\n".format(
                                Physic.constrains[row][0].astype(int),
                                bc_dof,
                                Physic.constrains[row][2],
                                Physic.constrains[row][3].astype(int),
                            )
                        )

            if "lo_list" in log_data["get"].keys():
                file_object.write("\n")
                file_object.write("LIST OF LOADS\n")
                file_object.write(
                    "{0:<7}{1:<20}{2:<16}{3:<10}\n".format(
                        "NODE", "TYPE", "VALUE", "STEP"
                    )
                )
                model_dict = Model.modelinfo["dofs"]["f"]
                for row in range(len(Physic.forces)):
                    fc_type = Physic.forces[row][1]

                    fc_dof = next((k for k, v in model_dict.items() if v == fc_type), None)

                    file_object.write(
                        "{0:<7}{1:<20}{2:<16}{3:<10}\n".format(
                            Physic.forces[row][0].astype(int),
                            fc_dof,
                            Physic.forces[row][2].round(decimals=np_decimals),
                            Physic.forces[row][3].astype(int),
                        )
                    )

            if "u_list" in log_data["get"].keys():
                file_object.write("\n")
                file_object.write("LIST OF SOLUTIONS\n")
                model_dict = list(Model.modelinfo["dofs"]["d"].keys()) 
                for row in range(len(solstatus['SOLUTION'])):
                    file_object.write("------------------------\n")
                    file_object.write("{0:<7}{1:<4}\n".format("STEP", row + 1))
                    file_object.write('dof, ' + str(model_dict)[1:-1])
                    file_object.write("\n")
                    array_sol = solstatus['SOLUTION'][row][Model.element.setTitleDeformation()][:, :len(model_dict)]
                    line_dof = 0
                    for line in array_sol:
                        line_str = ", ".join(map(str, line.round(decimals=np_decimals)))
                        file_object.write(str(line_dof) + ', ' + line_str + "\n")
                        line_dof+=1
                
            else:
                file_object.write("\n")
                file_object.write("log_data['get']: key erro\n")
                pass

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
