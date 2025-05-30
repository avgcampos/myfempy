#!/usr/bin/env python
__doc__ = """
GMSH GEN MESH
"""
import os

# import numpy as np
from numpy import abs


def meshid2gmshid(elemid):
    # https://gmsh.info/dev/doc/texinfo/gmsh.pdf
    # space + dofnode + numbconecelem + firstorder(1)/secondorder(2)
    celltype = {
        "1121": 1,
        "1132": 8,
        "1621": 1,
        "1632": 8,
        "2131": 2,
        "2162": 9,
        "2141": 3,
        "2182": 16,
        "2231": 2,
        "2262": 9,
        "2241": 3,
        "2282": 16,
        "3141": 4,
        "3341": 4,
        "33102": 11,
        "3181": 5,
        "3381": 5,
        "33202": 17,
    }
    gmshtyp = celltype[elemid]
    return gmshtyp


def gmsh_key(meshtype: str):
    l = {
        "line2": "-1",
        "line3": "-1",
        "tria3": "-2",
        "tria6": "-2",
        "quad4": "-2",
        "quad8": "-2",
        "tetr4": "-3",
        "tetr10": "-3",
        "hexa8": "-3",
        "hexa20": "-3",
    }
    return l[meshtype]


def get_gmsh_msh(filename, meshdata):
    os.system("echo MESHING...")
    cmd = (
        "gmsh"
        + " "
        + (filename + ".geo")
        + " "
        + gmsh_key(meshdata["meshconfig"]["mesh"])
        + " -o "
        + (filename + ".msh1")
    )
    os.system("echo GENERATING MESH FROM EXTERNAL GMSH")
    os.system(cmd)
    os.system("echo MESH IS DONE")
    # os.system("echo SAVING")


def get_gmsh_geo(filename, meshdata):
    with open((filename + ".geo"), "w") as file_object:
        file_object.write("// GMSH GEOMETRY FILE FROM MYFEMPY\n")
        file_object.write('SetFactory("OpenCASCADE");\n')
        if "pointlist" in meshdata.keys():
            numlinlist = len(meshdata["linelist"])
            line_list = ""
            for i in range(numlinlist):
                line_list += str(i + 1) + ","
            line_list = line_list[0:-1]
            if "planelist" in meshdata.keys():
                numplalistP = len(meshdata["planelist"])
                planes = ""
                for i in range(numplalistP):
                    planes += str(i + 1) + ","
                planes = planes[0:-1]
            else:
                pass
            numpnt = len(meshdata["pointlist"])
            for i in range(0, numpnt):
                file_object.write(
                    "Point("
                    + str(i + 1)
                    + ") = {"
                    + str(meshdata["pointlist"][i][0])
                    + ","
                    + str(meshdata["pointlist"][i][1])
                    + ","
                    + str(meshdata["pointlist"][i][2])
                    # + ","
                    # + str(meshdata["meshconfig"]["sizeelement"])
                    + "};"
                    + "\n"
                )
            if "circle" in meshdata.keys():
                numincl = len(meshdata["circle"])
                for inl in range(numincl):
                    d = meshdata["circle"][inl][0]
                    cx = meshdata["circle"][inl][1][0]
                    cy = meshdata["circle"][inl][1][1]
                    cz = meshdata["circle"][inl][1][2]
                    arc0 = meshdata["circle"][inl][2][0]
                    arc1 = meshdata["circle"][inl][2][1]
                    file_object.write(
                        "Circle("
                        + str(numlinlist + inl + 1)
                        + ") = {"
                        + str(cx)
                        + ","
                        + str(cy)
                        + ","
                        + str(cz)
                        + ","
                        + str(d)
                        + ","
                        + arc0
                        + ","
                        + arc1
                        + "};\n"
                    )
            else:
                numincl = 0
            if "arc" in meshdata.keys():
                numarcs = len(meshdata["arc"])
                for iarc in range(numarcs):
                    file_object.write(
                        "Circle("
                        + str(numlinlist + numincl + iarc + 1)
                        + ") = {"
                        + str(meshdata["arc"][iarc][0])
                        + ","
                        + str(meshdata["arc"][iarc][1])
                        + ","
                        + str(meshdata["arc"][iarc][2])
                        + "};\n"
                    )
                else:
                    numarcs = 0
            for i in range(0, numlinlist):
                file_object.write(
                    "Line("
                    + str(i + 1)
                    + ") = {"
                    + str(meshdata["linelist"][i][0])
                    + ","
                    + str(meshdata["linelist"][i][1])
                    + "};\n"
                )
        if meshdata["meshconfig"]["mesh"] == "line2":
            if "numbernodes" in meshdata["meshconfig"].keys():
                file_object.write(
                    "Transfinite Curve {"
                    + line_list
                    + "} = "
                    + str(meshdata["meshconfig"]["numbernodes"])
                    + " Using Progression 1;\n"
                )
            elif "sizeelement" in meshdata["meshconfig"].keys():
                for i in range(0, numpnt):
                    file_object.write(
                        "Point("
                        + str(i + 1)
                        + ") = {"
                        + str(meshdata["pointlist"][i][0])
                        + ","
                        + str(meshdata["pointlist"][i][1])
                        + ","
                        + str(meshdata["pointlist"][i][2])
                        + ","
                        + str(meshdata["meshconfig"]["sizeelement"])
                        + "};"
                        + "\n"
                    )
            file_object.write("// MESH CONFIGURATION\n")
            file_object.write("Mesh.CharacteristicLengthExtendFromBoundary = 1;\n")
            file_object.write("Mesh.CharacteristicLengthMin = 0;\n")
            file_object.write("Mesh.CharacteristicLengthFromPoints = 1;\n")
            file_object.write("Mesh.Optimize = 1;\n")
            file_object.write("Mesh.HighOrderOptimize = 0;\n")
            file_object.write("Mesh.Algorithm = 8;\n")
            file_object.write("Mesh.ElementOrder = 1;\n")

        elif meshdata["meshconfig"]["mesh"] == "line3":
            if "numbernodes" in meshdata["meshconfig"].keys():
                file_object.write(
                    "Transfinite Curve {"
                    + line_list
                    + "} = "
                    + str(meshdata["meshconfig"]["numbernodes"])
                    + " Using Progression 1;\n"
                )
            elif "sizeelement" in meshdata["meshconfig"].keys():
                for i in range(0, numpnt):
                    file_object.write(
                        "Point("
                        + str(i + 1)
                        + ") = {"
                        + str(meshdata["pointlist"][i][0])
                        + ","
                        + str(meshdata["pointlist"][i][1])
                        + ","
                        + str(meshdata["pointlist"][i][2])
                        + ","
                        + str(meshdata["meshconfig"]["sizeelement"])
                        + "};"
                        + "\n"
                    )
            file_object.write("// MESH CONFIGURATION\n")
            file_object.write("Mesh.CharacteristicLengthExtendFromBoundary = 1;\n")
            file_object.write("Mesh.CharacteristicLengthMin = 0;\n")
            file_object.write("Mesh.CharacteristicLengthFromPoints = 1;\n")
            file_object.write("Mesh.Optimize = 1;\n")
            file_object.write("Mesh.HighOrderOptimize = 0;\n")
            file_object.write("Mesh.Algorithm = 8;\n")
            file_object.write("Mesh.SecondOrderIncomplete = 1;\n")
            file_object.write("Mesh.ElementOrder = 2;\n")

        else:
            if "cadimport" in meshdata.keys():
                file_object.write('Merge "' + meshdata["cadimport"] + '";\n')
            else:
                npl = 0
                phl = 0
                for i in range(0, numplalistP):
                    npl += 1
                    file_object.write(
                        "Curve Loop("
                        + str(npl)
                        + ") = {"
                        + (str(abs(meshdata["planelist"][i][:]).tolist()))[1:-1]
                        + "};\n"
                    )

                npladd = 0
                lplrm = []
                nplrm = 0
                for i in range(0, numplalistP):
                    if any(jj < 0 for jj in meshdata["planelist"][i][:]):
                        nplrm += 1
                        lplrm.append(npladd + i)
                    else:
                        npladd += 1
                        lplrm.append(npladd)

                addpl = 0
                if nplrm > 0:
                    addpl = 1
                    file_object.write(
                        "Plane Surface("
                        + str(addpl)
                        + ") = {"
                        + str(list(set(lplrm)))[1:-1]
                        + "};\n"
                    )
                if npladd >= 1:
                    for iap in range(addpl, npladd):
                        file_object.write(
                            "Plane Surface("
                            + str(iap + 1)
                            + ") = {"
                            + str(addpl + iap + 1)
                            + "};\n"
                        )

                file_object.write(
                    "Characteristic Length {:} = "
                    + str(meshdata["meshconfig"]["sizeelement"])
                    + ";\n"
                )

            if "meshmap" in meshdata["meshconfig"].keys():
                if meshdata["meshconfig"]["meshmap"]["on"] == True:
                    file_object.write("//FACE MAPPING \n")
                    if "numbernodes" in meshdata["meshconfig"]["meshmap"].keys():
                        if meshdata["meshconfig"]["meshmap"]["edge"] == "all":
                            file_object.write(
                                "Transfinite Curve {:} = "
                                + str(meshdata["meshconfig"]["meshmap"]["numbernodes"])
                                + " Using Progression 1;\n"
                            )
                        else:
                            for ed in range(
                                len(meshdata["meshconfig"]["meshmap"]["edge"])
                            ):
                                file_object.write(
                                    "Transfinite Curve {"
                                    + str(
                                        meshdata["meshconfig"]["meshmap"]["edge"][ed]
                                    )[1:-1]
                                    + "} = "
                                    + str(
                                        meshdata["meshconfig"]["meshmap"][
                                            "numbernodes"
                                        ][ed]
                                    )
                                    + " Using Progression 1;\n"
                                )
                    file_object.write("Transfinite Surface {:};\n")
                else:
                    pass
            else:
                pass
            file_object.write("// MESH "+meshdata["meshconfig"]["mesh"]+" CONFIGURATION\n")
            if meshdata["meshconfig"]["mesh"] == "tria3":
                file_object.write("Mesh.CharacteristicLengthExtendFromBoundary = 1;\n")
                file_object.write("Mesh.CharacteristicLengthMin = 0;\n")
                file_object.write(
                    "Mesh.CharacteristicLengthMax = "
                    + str(meshdata["meshconfig"]["sizeelement"])
                    + ";\n"
                )
                file_object.write("Mesh.CharacteristicLengthFromPoints = 1;\n")
                file_object.write("Mesh.Optimize = 1;\n")
                file_object.write("Mesh.HighOrderOptimize = 0;\n")
                file_object.write("Mesh.Algorithm = 8;\n")
                file_object.write("Mesh.ElementOrder = 1;\n")

            elif meshdata["meshconfig"]["mesh"] == "tria6":
                file_object.write("Mesh.CharacteristicLengthExtendFromBoundary = 1;\n")
                file_object.write("Mesh.CharacteristicLengthMin = 0;\n")
                file_object.write(
                    "Mesh.CharacteristicLengthMax = "
                    + str(meshdata["meshconfig"]["sizeelement"])
                    + ";\n"
                )
                file_object.write("Mesh.CharacteristicLengthFromPoints = 1;\n")
                file_object.write("Mesh.Optimize = 1;\n")
                file_object.write("Mesh.HighOrderOptimize = 0;\n")
                file_object.write("Mesh.Algorithm = 8;\n")
                file_object.write("Mesh.SecondOrderIncomplete = 1;\n")
                file_object.write("Mesh.ElementOrder = 2;\n")

            elif meshdata["meshconfig"]["mesh"] == "quad4":
                file_object.write("Recombine Surface {:};\n")
                file_object.write("Mesh.RecombinationAlgorithm = 1;\n")
                file_object.write("Mesh.RecombineAll = 1;\n")
                file_object.write("Mesh.SubdivisionAlgorithm = 1;\n")
                file_object.write("Mesh.CharacteristicLengthExtendFromBoundary = 1;\n")
                file_object.write("Mesh.CharacteristicLengthMin = 0;\n")
                file_object.write(
                    "Mesh.CharacteristicLengthMax = "
                    + str(meshdata["meshconfig"]["sizeelement"])
                    + ";\n"
                )
                file_object.write("Mesh.CharacteristicLengthFromPoints = 1;\n")
                file_object.write("Mesh.Optimize = 1;\n")
                file_object.write("Mesh.HighOrderOptimize = 0;\n")
                file_object.write("Mesh.Algorithm = 8;\n")
                file_object.write("Mesh.ElementOrder = 1;\n")

            elif meshdata["meshconfig"]["mesh"] == "quad8":
                file_object.write("Recombine Surface {:};\n")
                file_object.write("Mesh.RecombinationAlgorithm = 1;\n")
                file_object.write("Mesh.RecombineAll = 1;\n")
                file_object.write("Mesh.SubdivisionAlgorithm = 1;\n")
                file_object.write("Mesh.CharacteristicLengthExtendFromBoundary = 1;\n")
                file_object.write("Mesh.CharacteristicLengthMin = 0;\n")
                file_object.write(
                    "Mesh.CharacteristicLengthMax = "
                    + str(meshdata["meshconfig"]["sizeelement"])
                    + ";\n"
                )
                file_object.write("Mesh.CharacteristicLengthFromPoints = 1;\n")
                file_object.write("Mesh.Optimize = 1;\n")
                file_object.write("Mesh.HighOrderOptimize = 0;\n")
                file_object.write("Mesh.Algorithm = 8;\n")
                file_object.write("Mesh.SecondOrderIncomplete = 1;\n")
                file_object.write("Mesh.ElementOrder = 2;\n")

            elif meshdata["meshconfig"]["mesh"] == "tetr4":
                if "extrude" in meshdata["meshconfig"].keys():
                    thck = meshdata["meshconfig"]["extrude"]
                    file_object.write(
                        "Extrude {0, 0, " + str(float(thck)) + "} {Surface{:};}\n"
                    )
                file_object.write("Mesh.Algorithm = 2;\n")  # 4
                file_object.write("Mesh.Algorithm3D = 4;\n")  # 7
                file_object.write("Mesh.CharacteristicLengthExtendFromBoundary = 1;\n")
                file_object.write("Mesh.CharacteristicLengthMin = 0;\n")
                file_object.write(
                    "Mesh.CharacteristicLengthMax = "
                    + str(meshdata["meshconfig"]["sizeelement"])
                    + ";\n"
                )
                file_object.write("Mesh.ElementOrder = 1;\n")
                file_object.write("Mesh.Optimize = 1;\n")
                file_object.write("Mesh.HighOrderOptimize = 0;\n")

            elif meshdata["meshconfig"]["mesh"] == "hexa8":
                if "extrude" in meshdata["meshconfig"].keys():
                    thck = meshdata["meshconfig"]["extrude"]
                    # file_object.write('Extrude {0, 0, '+str(float(thck))+'} {Surface{:};Layers{'+str(int(float(thck)/float(meshdata['GMSH']['meshconfig']['sizeelement'])))+'};Recombine;};\n')
                    file_object.write(
                        "Extrude {0, 0, "
                        + str(float(thck))
                        + "} {Surface{:};Layers{"
                        + str(
                            int(
                                float(thck)
                                / float(meshdata["meshconfig"]["sizeelement"])
                            )
                        )
                        + "};Recombine;};\n"
                    )
                file_object.write("Recombine Surface {:};\n")
                file_object.write("Mesh.Algorithm = 2;\n")  # 8
                file_object.write("Mesh.Algorithm3D = 4;\n")  # 7
                file_object.write("Mesh.CharacteristicLengthExtendFromBoundary = 1;\n")
                file_object.write("Mesh.CharacteristicLengthMin = 0;\n")
                file_object.write(
                    "Mesh.CharacteristicLengthMax = "
                    + str(meshdata["meshconfig"]["sizeelement"])
                    + ";\n"
                )
                file_object.write("Mesh.ElementOrder = 1;\n")
                file_object.write("Mesh.Optimize = 1;\n")
                file_object.write("Mesh.HighOrderOptimize = 0;\n")
                file_object.write("Mesh.RecombinationAlgorithm = 0;\n")
                file_object.write("Mesh.SubdivisionAlgorithm = 2;\n")
                file_object.write("Mesh.RecombineAll = 1;\n")
