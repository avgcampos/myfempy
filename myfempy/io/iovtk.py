#!/usr/bin/env python
"""
I/O VTK
"""
import numpy as np


def meshid2vtkid(elemid):
    # https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    # space + dofnode + numbconecelem + firstorder(1)/secondorder(2)
    CellType = {
        "1621": 3,
        "1632": 21,
        "2231": 5,
        "2331": 5,
        "2262": 22,
        "2241": 9,
        "2341": 9,
        "2282": 23,
        "3341": 10,
        "33102": 24,
        "3381": 12,
        "33202": 25,
    }
    vtkCT = CellType[elemid]
    return vtkCT


def convert_to_vtk(plotdata):
    """_summary_

    Arguments:
        plotdata -- _description_
    """
    # ii: cython.int
    # jj: cython.int
    # md: cython.int
    # numnodes: cython.int
    # numelem: cython.int
    # numdsp_pntdata: cython.int
    # numstr_pntdata: cython.int
    # nummds_pntdata: cython.int
    # num_modes: cython.int
    # numstr_clldata: cython.int

    numnodes = int(len(plotdata["coord"]))
    numelem = int(len(plotdata["inci"]))
    numdsp_pntdata = int(len(plotdata["displ_POINT_DATA_val"]))
    numstr_pntdata = int(len(plotdata["stress_POINT_DATA_title"]))
    nummds_pntdata = int(len(plotdata["modes_POINT_DATA"]))

    numstr_clldata = int(len(plotdata["stress_CELL_DATA_title"]))
    # numstreng_clldata = int(len(plotdata["strain_energy_CELL_DATA_title"]))

    with open(plotdata["filename"] + ".vtk", "w") as file_object:
        file_object.write("# vtk DataFile Version 4.0\n")
        file_object.write("vtk output from myfempy solver\n")
        file_object.write("ASCII\n")
        file_object.write("DATASET UNSTRUCTURED_GRID\n")
        file_object.write("POINTS " + str(int(len(plotdata["coord"]))) + " double\n")
        for ii in range(0, numnodes):
            list2write = plotdata["coord"][ii, 1:].astype(str).tolist()
            file_object.write(" ".join(list2write) + "\n")
        file_object.write("\n")
        file_object.write(
            "CELLS"
            + " "
            + str(int(len(plotdata["inci"])))
            + " "
            + str((plotdata["nodecon"] + 1) * int(len(plotdata["inci"])))
            + "\n"
        )
        for ii in range(0, numelem):
            listnodes = plotdata["inci"][ii, 4:]
            list2write = (
                np.array(list(map(lambda x: x - 1, listnodes[listnodes.nonzero()])))
                .astype(int)
                .astype(str)
                .tolist()
            )
            file_object.write(str(len(list2write)) + " " + " ".join(list2write) + "\n")
        file_object.write("\n")
        file_object.write("CELL_TYPES" + " " + str(int(len(plotdata["inci"]))) + "\n")
        for ii in range(0, numelem):
            vtkCT = meshid2vtkid(str(int(plotdata["inci"][ii, 1])))
            file_object.write(str(vtkCT) + "\n")
        file_object.write("\n")
        # -----------------------------------------------------
        # POINT DATA
        file_object.write("POINT_DATA" + " " + str(int(len(plotdata["coord"]))) + "\n")
        if len(plotdata["stress_POINT_DATA_val"]) > 0:
            for jj in range(0, numstr_pntdata):
                file_object.write(
                    "SCALARS " + plotdata["stress_POINT_DATA_title"][jj] + " float 1\n"
                )
                file_object.write("LOOKUP_TABLE default\n")
                list2write = (
                    plotdata["stress_POINT_DATA_val"][:, jj].astype(str).tolist()
                )
                file_object.write("\n".join(list2write) + "\n")
                file_object.write("\n")
        file_object.write("\n")

        if len(plotdata["displ_POINT_DATA_val"]) > 0:
            file_object.write("FIELD FieldData 1\n")
            file_object.write(
                plotdata["displ_POINT_DATA_title"]
                + " 3 "
                + str(int(len(plotdata["displ_POINT_DATA_val"])))
                + " float\n"
            )
            for ii in range(0, numdsp_pntdata):
                list2write = (
                    plotdata["displ_POINT_DATA_val"][ii, :].astype(str).tolist()
                )
                file_object.write(" ".join(list2write) + "\n")

        if len(plotdata["modes_POINT_DATA"]) > 0:
            for md in range(0, nummds_pntdata):
                file_object.write("\n")
                file_object.write("FIELD FieldData 1\n")
                file_object.write(
                    plotdata["modes_POINT_DATA"][md]["title"]
                    + " 3 "
                    + str(int(len(plotdata["modes_POINT_DATA"][md]["val"][:, 1:])))
                    + " float\n"
                )
                for ii in range(
                    0, int(len(plotdata["modes_POINT_DATA"][md]["val"][:, 1:]))
                ):
                    list2write = (
                        plotdata["modes_POINT_DATA"][md]["val"][:, 1:][ii, :]
                        .astype(str)
                        .tolist()
                    )
                    file_object.write(" ".join(list2write) + "\n")

        # POINT DATA FUTURE...
        file_object.write("\n")
        # -----------------------------------------------------
        # CELL DATA
        file_object.write("CELL_DATA" + " " + str(int(len(plotdata["inci"]))) + "\n")
        if len(plotdata["stress_CELL_DATA_val"]) > 0:
            for jj in range(0, numstr_clldata):
                file_object.write(
                    "SCALARS " + plotdata["stress_CELL_DATA_title"][jj] + " float 1\n"
                )
                file_object.write("LOOKUP_TABLE default\n")
                list2write = (
                    plotdata["stress_CELL_DATA_val"][:, jj].astype(str).tolist()
                )
                file_object.write("\n".join(list2write) + "\n")
        # if len(plotdata["strain_energy_CELL_DATA_val"]) > 0:
        #     for jj in range(0, numstreng_clldata):
        #         file_object.write(
        #             "SCALARS " + plotdata["strain_energy_CELL_DATA_title"][jj] + " float 1\n"
        #         )
        #         file_object.write("LOOKUP_TABLE default\n")
        #         list2write = (
        #             plotdata["strain_energy_CELL_DATA_val"][:, jj].astype(str).tolist()
        #         )
        #         file_object.write("\n".join(list2write) + "\n")

        # CELL DATA FUTURE...


def convert_from_vtk(filename):
    """_summary_

    Arguments:
        filename -- _description_

    Returns:
        _description_
    """

    file_imp = filename + ".vtk"
    with open(file_imp, "r") as file_object:
        file_object.readline()
        file_object.readline()
        file_object.readline()
        file_object.readline()
        line = file_object.readline()
        lineaux = line.split()
        nnod = int(lineaux[1])
        nodelist = [[None] * 4]
        for ii in range(0, nnod):
            line = file_object.readline()
            lineaux = line.split()
            contstr = lineaux[0:3]
            nodelist.append(
                [int(ii + 1), float(contstr[0]), float(contstr[1]), float(contstr[2])]
            )
        nodelist = nodelist[1::][::]
        file_object.readline()
        line = file_object.readline()
        lineaux = line.split()
        nelm = int(lineaux[1])
        conec_elm = []
        for kk in range(0, nelm):
            line = file_object.readline()
            lineaux = line.split()
            conec_elm.append(list(map(float, lineaux[:])))
        conec_elm = conec_elm[1::][::]
    return conec_elm, nodelist
