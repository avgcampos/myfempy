#!/usr/bin/env python
__doc__ = """
GMSH GEN MESH
"""
import os


def gmsh_key(meshtype: str):
    l = {
        "line2": "-1",
        "tria3": "-2",
        "quad4": "-2",
        "hexa8": "-3",
        "tetr4": "-3",
    }
    return l[meshtype]


def get_gmsh_msh(meshdata: dict):
    cmd = (
        "gmsh"
        + " "
        + (meshdata["GMSH"]["filename"] + ".geo")
        + " "
        + gmsh_key(meshdata["GMSH"]["meshconfig"]["mesh"])
        + " -o "
        + (meshdata["GMSH"]["filename"] + ".msh1")
    )
    # os.system("echo GENERATING MESH FROM EXTERNAL GMSH")
    os.system(cmd)
    os.system("echo MESHING IS DONE")
    # os.system("echo SAVING AND EXIT")


def get_gmsh_geo(meshdata: dict):
    with open((meshdata["GMSH"]["filename"] + ".geo"), "w") as file_object:
        file_object.write("// GMSH GEOMETRY FILE FROM MYFEMPY\n")
        file_object.write('SetFactory("OpenCASCADE");\n')
        if "pointlist" in meshdata["GMSH"].keys():
            numlinlist = len(meshdata["GMSH"]["linelist"])
            line_list = ""
            for i in range(numlinlist):
                line_list += str(i + 1) + ","
            line_list = line_list[0:-1]
            if "planelist" in meshdata["GMSH"].keys():
                numplalistP = len(meshdata["GMSH"]["planelist"])
                planes = ""
                for i in range(numplalistP):
                    planes += str(i + 1) + ","
                planes = planes[0:-1]
            else:
                pass
            numpnt = len(meshdata["GMSH"]["pointlist"])
            for i in range(0, numpnt):
                file_object.write(
                    "Point("
                    + str(i + 1)
                    + ") = {"
                    + str(meshdata["GMSH"]["pointlist"][i][0])
                    + ","
                    + str(meshdata["GMSH"]["pointlist"][i][1])
                    + ","
                    + str(meshdata["GMSH"]["pointlist"][i][2])
                    + ","
                    + str(meshdata["GMSH"]["meshconfig"]["sizeelement"])
                    + "};"
                    + "\n"
                )
            if "arc" in meshdata["GMSH"].keys():
                numincl = len(meshdata["GMSH"]["arc"])
                for inl in range(numincl):
                    d = meshdata["GMSH"]["arc"][inl][0]
                    cx = meshdata["GMSH"]["arc"][inl][1][0]
                    cy = meshdata["GMSH"]["arc"][inl][1][1]
                    cz = meshdata["GMSH"]["arc"][inl][1][2]
                    arc0 = meshdata["GMSH"]["arc"][inl][2][0]
                    arc1 = meshdata["GMSH"]["arc"][inl][2][1]
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
            for i in range(0, numlinlist):
                file_object.write(
                    "Line("
                    + str(i + 1)
                    + ") = {"
                    + str(meshdata["GMSH"]["linelist"][i][0])
                    + ","
                    + str(meshdata["GMSH"]["linelist"][i][1])
                    + "};\n"
                )
        if meshdata["GMSH"]["meshconfig"]["mesh"] == "line2":
            if "numbernodes" in meshdata["GMSH"]["meshconfig"].keys():
                file_object.write(
                    "Transfinite Curve {"
                    + line_list
                    + "} = "
                    + str(meshdata["GMSH"]["meshconfig"]["numbernodes"])
                    + " Using Progression 1;\n"
                )
            elif "sizeelement" in meshdata["GMSH"]["meshconfig"].keys():
                for i in range(0, numpnt):
                    file_object.write(
                        "Point("
                        + str(i + 1)
                        + ") = {"
                        + str(meshdata["GMSH"]["pointlist"][i][0])
                        + ","
                        + str(meshdata["GMSH"]["pointlist"][i][1])
                        + ","
                        + str(meshdata["GMSH"]["pointlist"][i][2])
                        + ","
                        + str(meshdata["GMSH"]["meshconfig"]["sizeelement"])
                        + "};"
                        + "\n"
                    )
            else:
                pass
        elif meshdata["GMSH"]["meshconfig"]["mesh"] != "line2":
            if "cadimport" in meshdata["GMSH"].keys():
                file_object.write(
                    'Merge "' + meshdata["GMSH"]["cadimport"]["object"] + '";\n'
                )
            else:
                npl = 0
                phl = 0
                for i in range(0, numplalistP):
                    npl += 1
                    file_object.write(
                        "Curve Loop("
                        + str(npl)
                        + ") = {"
                        + (str(meshdata["GMSH"]["planelist"][i][:]))[1:-1]
                        + "};\n"
                    )
                    plane_hole = ""
                    plane_hole += str(npl) + ","
                    for inl in range(numincl):
                        if meshdata["GMSH"]["arc"][inl][0] == i + 1:
                            npl += 1
                            plane_hole += str(npl) + ","
                            phl = meshdata["GMSH"]["arc"][inl][1]
                            file_object.write(
                                "Curve Loop("
                                + str(i + 1 + inl + 1)
                                + ") = {"
                                + str(numlinlist + inl + 1)
                                + "};\n"
                            )
                    if "arc" in meshdata["GMSH"].keys():
                        file_object.write(
                            "Plane Surface("
                            + str(i + 1)
                            + ") = {"
                            + plane_hole[0:-1]
                            + "};\n"
                        )
                    else:
                        file_object.write(
                            "Plane Surface(" + str(i + 1) + ") = {" + str(npl) + "};\n"
                        )
                file_object.write(
                    "Characteristic Length {:} = "
                    + str(meshdata["GMSH"]["meshconfig"]["sizeelement"])
                    + ";\n"
                )
            if meshdata["GMSH"]["meshconfig"]["meshmap"]["on"] == True:
                file_object.write("//FACE MAPPING \n")
                if "numbernodes" in meshdata["GMSH"]["meshconfig"]["meshmap"].keys():
                    if meshdata["GMSH"]["meshconfig"]["meshmap"]["edge"] == "all":
                        file_object.write(
                            "Transfinite Curve {:} = "
                            + str(
                                meshdata["GMSH"]["meshconfig"]["meshmap"]["numbernodes"]
                            )
                            + " Using Progression 1;\n"
                        )
                    else:
                        file_object.write(
                            "Transfinite Curve {"
                            + str(meshdata["GMSH"]["meshconfig"]["meshmap"]["edge"])[
                                1:-1
                            ]
                            + "} = "
                            + str(
                                meshdata["GMSH"]["meshconfig"]["meshmap"]["numbernodes"]
                            )
                            + " Using Progression 1;\n"
                        )
                file_object.write("Transfinite Surface {:};\n")
            else:
                pass
            if meshdata["GMSH"]["meshconfig"]["mesh"] == "tria3":
                file_object.write("// MESH CONFIGURATION\n")
                file_object.write("Mesh.CharacteristicLengthExtendFromBoundary = 1;\n")
                file_object.write("Mesh.CharacteristicLengthMin = 0;\n")
                file_object.write(
                    "Mesh.CharacteristicLengthMax = "
                    + str(meshdata["GMSH"]["meshconfig"]["sizeelement"])
                    + ";\n"
                )
                file_object.write("Mesh.CharacteristicLengthFromPoints = 1;\n")
                file_object.write("Mesh.Optimize = 1;\n")
                file_object.write("Mesh.HighOrderOptimize = 0;\n")
                file_object.write("Mesh.Algorithm = 8;\n")
                file_object.write("Mesh.ElementOrder = 1;\n")
            elif meshdata["GMSH"]["meshconfig"]["mesh"] == "quad4":
                file_object.write("Recombine Surface {:};\n")
                file_object.write("// MESH CONFIGURATION\n")
                file_object.write("Mesh.RecombinationAlgorithm = 1;\n")
                file_object.write("Mesh.RecombineAll = 1;\n")
                file_object.write("Mesh.SubdivisionAlgorithm = 1;\n")
                file_object.write("Mesh.CharacteristicLengthExtendFromBoundary = 1;\n")
                file_object.write("Mesh.CharacteristicLengthMin = 0;\n")
                file_object.write(
                    "Mesh.CharacteristicLengthMax = "
                    + str(meshdata["GMSH"]["meshconfig"]["sizeelement"])
                    + ";\n"
                )
                file_object.write("Mesh.CharacteristicLengthFromPoints = 1;\n")
                file_object.write("Mesh.Optimize = 1;\n")
                file_object.write("Mesh.HighOrderOptimize = 0;\n")
                file_object.write("Mesh.Algorithm = 8;\n")
                file_object.write("Mesh.ElementOrder = 1;\n")
            elif meshdata["GMSH"]["meshconfig"]["mesh"] == "tetr4":
                if "extrude" in meshdata["GMSH"]["meshconfig"].keys():
                    thck = meshdata["GMSH"]["meshconfig"]["extrude"]
                    file_object.write(
                        "Extrude {0, 0, " + str(float(thck)) + "} {Surface{:};}\n"
                    )
                file_object.write("// MESH CONFIGURATION\n")
                file_object.write("Mesh.Algorithm = 2;\n")  # 4
                file_object.write("Mesh.Algorithm3D = 4;\n")  # 7
                file_object.write("Mesh.CharacteristicLengthExtendFromBoundary = 1;\n")
                file_object.write("Mesh.CharacteristicLengthMin = 0;\n")
                file_object.write(
                    "Mesh.CharacteristicLengthMax = "
                    + str(meshdata["GMSH"]["meshconfig"]["sizeelement"])
                    + ";\n"
                )
                file_object.write("Mesh.ElementOrder = 1;\n")
                file_object.write("Mesh.Optimize = 1;\n")
                file_object.write("Mesh.HighOrderOptimize = 0;\n")
            elif meshdata["GMSH"]["meshconfig"]["mesh"] == "hexa8":
                if "extrude" in meshdata["GMSH"]["meshconfig"].keys():
                    thck = meshdata["GMSH"]["meshconfig"]["extrude"]
                    # file_object.write('Extrude {0, 0, '+str(float(thck))+'} {Surface{:};Layers{'+str(int(float(thck)/float(meshdata['GMSH']['meshconfig']['sizeelement'])))+'};Recombine;};\n')
                    file_object.write(
                        "Extrude {0, 0, "
                        + str(float(thck))
                        + "} {Surface{:};Layers{"
                        + str(
                            int(
                                float(thck)
                                / float(meshdata["GMSH"]["meshconfig"]["sizeelement"])
                            )
                        )
                        + "};Recombine;};\n"
                    )
                file_object.write("Recombine Surface {:};\n")
                file_object.write("// MESH CONFIGURATION\n")
                file_object.write("Mesh.Algorithm = 2;\n")  # 8
                file_object.write("Mesh.Algorithm3D = 4;\n")  # 7
                file_object.write("Mesh.CharacteristicLengthExtendFromBoundary = 1;\n")
                file_object.write("Mesh.CharacteristicLengthMin = 0;\n")
                file_object.write(
                    "Mesh.CharacteristicLengthMax = "
                    + str(meshdata["GMSH"]["meshconfig"]["sizeelement"])
                    + ";\n"
                )
                file_object.write("Mesh.ElementOrder = 1;\n")
                file_object.write("Mesh.Optimize = 1;\n")
                file_object.write("Mesh.HighOrderOptimize = 0;\n")
                file_object.write("Mesh.RecombinationAlgorithm = 0;\n")
                file_object.write("Mesh.SubdivisionAlgorithm = 2;\n")
                file_object.write("Mesh.RecombineAll = 1;\n")
        else:
            print("input erro: mesh_cfg don't defined")
