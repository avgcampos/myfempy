#!/usr/bin/env python
import csv
import numpy as np

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

def write2log(Model, Physic, log_data, solstatus, log_file):
    lines = [
        "===============================================================================",
        "                                 M Y F E M P Y                                 ",
        "version "+str(solstatus["solverstatus"]["myfempyversion"]),
        "===============================================================================",
    ]

    get_keys = log_data.get("get", {})
    np_decimals = get_keys.get("numpy_decimals", 4)

    if "nelem" in get_keys:
        lines.append("")
        lines.append(f"GET NELEM {len(Model.inci)}")

    if "nnode" in get_keys:
        lines.append("")
        lines.append(f"GET NNODE {len(Model.coord)}")

    if "inci" in get_keys:
        lines.append("")
        lines.append("LIST OF ELEMENTS")
        lines.append("{0:<7}{1:<7}{2:<7}{3:<7}{4:<25}".format("ELEM", "KEY", "MAT", "GEO", "NODES"))
        for row in range(len(Model.inci)):
            node_list = Model.inci[row][4:].astype(int)
            node_list = node_list[np.nonzero(node_list)].tolist()
            node_str = " ".join(str(x) for x in node_list)
            lines.append("{0:<7}{1:<7}{2:<7}{3:<7}{4:<25}".format(
                int(Model.inci[row][0]),
                int(Model.inci[row][1]),
                int(Model.inci[row][2]),
                int(Model.inci[row][3]),
                node_str
            ))

    if "coord" in get_keys:
        lines.append("")
        lines.append("LIST OF NODES COORDINATE")
        lines.append("{0:<7}{1:<10}{2:<10}{3:<10}".format("NODE", "X", "Y", "Z"))
        for row in range(len(Model.coord)):
            lines.append("{0:<7}{1:<10}{2:<10}{3:<10}".format(
                int(Model.coord[row][0]),
                round(float(Model.coord[row][1]), np_decimals),
                round(float(Model.coord[row][2]), np_decimals),
                round(float(Model.coord[row][3]), np_decimals)
            ))

    if "tabmat" in get_keys:
        lines.append("")
        lines.append("LIST OF MATERIAL PROPERTY")
        for row in range(len(Model.tabmat)):
            lines.append("------------------------")
            lines.append("{0:<12}{1:<3}".format("MATERIAL", row + 1))
            lines.append("{0:<10}{1:<10}".format("PROPMAT", "VALUE"))
            for key, value_propmat in Model.tabmat[row].items():
                if value_propmat != "NULL":
                    lines.append("{0:<10}{1:<10}".format(key, value_propmat))

    if "tabgeo" in get_keys:
        lines.append("")
        lines.append("LIST OF GEOMETRY PROPERTY")
        for row in range(len(Model.tabgeo)):
            lines.append("------------------------")
            lines.append("{0:<12}{1:<3}".format("GEOMETRY", row + 1))
            lines.append("{0:<10}{1:<10}".format("PROPGEO", "VALUE"))
            for key, value_propmat in Model.tabgeo[row].items():
                if value_propmat != "NULL":
                    lines.append("{0:<10}{1:<10}".format(key, value_propmat))

    if "bc_list" in get_keys:
        lines.append("")
        lines.append("LIST OF BOUNDARY CONDITIONS")
        lines.append("{0:<7}{1:<10}{2:<16}{3:<10}".format("NODE", "DOF", "VALUE", "STEP"))
        model_dict = Model.modelinfo["dofs"]["d"]
        for row in range(len(Physic.constrains)):
            bc_type = Physic.constrains[row][1]
            node_id = int(Physic.constrains[row][0])
            step_id = int(Physic.constrains[row][3])
            if bc_type == 0:
                val = round(float(Physic.constrains[row][2]), np_decimals)
                lines.append("{0:<7}{1:<10}{2:<16}{3:<10}".format(node_id, "full", val, step_id))
            else:
                bc_dof = next((k for k, v in model_dict.items() if v == bc_type), None)
                val = Physic.constrains[row][2]
                lines.append("{0:<7}{1:<10}{2:<16}{3:<10}".format(node_id, bc_dof, val, step_id))

    if "lo_list" in get_keys:
        lines.append("")
        lines.append("LIST OF LOADS")
        lines.append("{0:<7}{1:<20}{2:<16}{3:<10}".format("NODE", "TYPE", "VALUE", "STEP"))
        model_dict = Model.modelinfo["dofs"]["f"]
        for row in range(len(Physic.forces)):
            fc_type = Physic.forces[row][1]
            fc_dof = next((k for k, v in model_dict.items() if v == fc_type), None)
            val = round(float(Physic.forces[row][2]), np_decimals)
            lines.append("{0:<7}{1:<20}{2:<16}{3:<10}".format(
                int(Physic.forces[row][0]),
                str(fc_dof),
                val,
                int(Physic.forces[row][3])
            ))

    if "u_list" in get_keys:
        lines.append("")
        lines.append("LIST OF SOLUTIONS")
        model_dict = list(Model.modelinfo["dofs"]["d"].keys())
        dof_header = str(model_dict)[1:-1]
        for row in range(len(solstatus['SOLUTION'])):
            lines.append("------------------------")
            lines.append("{0:<7}{1:<4}".format("STEP", row + 1))
            lines.append(f"dof, {dof_header}")
            array_sol = solstatus['SOLUTION'][row][Model.element.setTitleDeformation()][:, :len(model_dict)]
            for line_dof, line in enumerate(array_sol):
                line_str = ", ".join(map(str, line.round(decimals=np_decimals)))
                lines.append(f"{line_dof}, {line_str}")

    if "log" in log_data.keys():
        lines.append("")
        lines.append("+---------------------------- S O L V E R   L O G ----------------------------+")
        lines.append("{0:<30} : {1:<10}".format("ANALYZED ON ", str(log_data["timesolver"])))
        lines.append("{0:<30} : {1:<10} SEC".format("ASSEMBLY FULL TIME SPEND ", str(solstatus["solverstatus"]["timeasb"])))
        lines.append("{0:<30} : {1:<10} SEC".format("SOLVE FULL TIME SPEND ", str(solstatus["solverstatus"]["timesim"])))
        lines.append("{0:<30} : {1:<10} DOF".format("NUMBER OF EQUATION ", str(Model.modelinfo["fulldofs"])))
        lines.append("{0:<30} : {1:<10} MB".format("STIFFNESS SIZE ", str(solstatus["solverstatus"]["memorysize"])))
        lines.append("{0:<30} : {1:<10} ".format("TYPE ASSEMBLER ", str(solstatus["solverstatus"]["typeasmb"])))
        lines.append("{0:<30} : {1:<10} ".format("SOLVER CORE ", str(solstatus["solverstatus"]["solvercore"])))

    with open(log_file, "w") as file_object:
        file_object.write("\n".join(lines) + "\n")


def writer2csv(csv2write_file, csv2write_data, label):
    with open(csv2write_file, "w", newline="") as file_object:
        writer = csv.writer(file_object)
        writer.writerow([label[0], label[1]])
        writer.writerows(zip(csv2write_data[0], csv2write_data[1]))


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
    lines = []
    for i in range(len(solvercfg)):
        lines.append(f"solver_typ = '{solvercfg[i, 1]}'")
        lines.append(f"solver_def = '{solvercfg[i, 2]}'")
        lines.append(f"solver_opt = '{solvercfg[i, 3]}'")
        lines.append(f"solver_start = '{solvercfg[i, 4]}'")
        lines.append(f"solver_end = '{solvercfg[i, 5]}'")
        lines.append(f"solver_step = '{solvercfg[i, 6]}'")

    for i in range(len(solutconfig)):
        lines.append(f"mod_typ = '{solutconfig[i, 1]}'")
        lines.append(f"mod_opt = '{solutconfig[i, 2]}'")

    for i in range(len(meshconfig)):
        lines.append(f"mesh_typ = '{meshconfig[i, 1]}'")

    for i in range(len(compmaterial)):
        lines.append(f"mat_typ = '{compmaterial[i, 1]}'")
        lines.append(f"mat_opt = '{compmaterial[i, 2]}'")
        lines.append(f"mat_def = '{compmaterial[i, 3]}'")

    unique_steps = np.unique(forcelist[:, 8])
    for j, step_val in enumerate(unique_steps):
        forceliststep = forcelist[np.where(forcelist[:, 8] == step_val), :][0]
        for i in range(len(forceliststep)):
            lines.append(f"force_typ_{i} = '{forceliststep[i, 1]}'")
            lines.append(f"force_opt_{i} = '{forceliststep[i, 4]}'")
            lines.append(f"force_opt_dir{i} = '{forceliststep[i, 2]}'")
            lines.append(f"force_step{j} = '{forceliststep[i, 8]}'")

    for i in range(len(boncdlist)):
        lines.append(f"bc_typ_{i} = '{boncdlist[i, 1]}'")
        lines.append(f"bc_opt_{i} = '{boncdlist[i, 3]}'")
        lines.append(f"bc_opt_dir{i} = '{boncdlist[i, 2]}'")

    lines.append(f"prev_sol = '{outputcfg[0, 1]}'")
    lines.append(f"graph_out = '{outputcfg[0, 2]}'")
    lines.append(f"save_fig = '{outputcfg[0, 3]}'")

    for i in range(len(graphout)):
        lines.append(f"graph_out_{i} = '{graphout[i]}'")

    for i in range(len(fileout)):
        lines.append(f"file_out_{i} = '{fileout[i]}'")

    with open(usrlog_file, "w") as file_object:
        file_object.write("\n".join(lines) + "\n")