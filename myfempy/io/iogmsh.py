#!/usr/bin/env python
import os

import gmsh
import sys
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


def get_mesh_gmsh(filename, meshdata):
    # os.system("echo MESHING...")
    cmd = (
        "gmsh"
        + " "
        + "-v 0"
        + " "
        + (filename + ".geo")
        + " "
        + gmsh_key(meshdata["meshconfig"]["mesh"])
        + " -o "
        + (filename + ".msh2")
    )
    # os.system("echo GENERATING MESH FROM EXTERNAL GMSH")
    os.system(cmd)
    # os.system("echo MESH IS DONE")

def get_reorder_mesh(filename, meshdata):
        

    gmsh.initialize()
    # --- Hide terminal output ---
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(filename + ".geo")

    gmsh.model.mesh.generate(abs(int(gmsh_key(meshdata["meshconfig"]["mesh"]))))

    # original_model = gmsh.model.getCurrent()

    # 1. Coletar dados e calcular ordenação por PESO GEOMÉTRICO (Z -> Y -> X)
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    coords = coords.reshape(-1, 3)
    c_round = np.round(coords, 6)

    # Score: Z (1e12) + Y (1e6) + X (1)
    sort_scores = c_round[:, 2] * 1e12 + c_round[:, 1] * 1e6 + c_round[:, 0]
    indices = np.argsort(sort_scores)

    # Mapeamento: ID_Antigo -> ID_Novo
    old_to_new = {int(node_tags[idx]): i + 1 for i, idx in enumerate(indices)}

    # 2. Coletar estrutura COMPLETA do modelo original
    physicals = []
    for dim, p_tag in gmsh.model.getPhysicalGroups():
        physicals.append((dim, p_tag, gmsh.model.getPhysicalName(dim, p_tag), gmsh.model.getEntitiesForPhysicalGroup(dim, p_tag)))

    entities_data = []
    node_ownership = {}
    for dim in range(4):
        for _, e_tag in gmsh.model.getEntities(dim):
            e_node_tags, _, _ = gmsh.model.mesh.getNodes(dim, e_tag)
            for nt in e_node_tags:
                nt_int = int(nt)
                if nt_int not in node_ownership:
                    node_ownership[nt_int] = (dim, e_tag)
            e_types, e_tags, e_conn = gmsh.model.mesh.getElements(dim, e_tag)
            entities_data.append((dim, e_tag, e_node_tags, e_types, e_tags, e_conn))

    # 3. CRIAR NOVO MODELO E RECONSTRUIR
    gmsh.model.add("MalhaReordenada")

    # IMPORTANTE: No MSH2, a ordem dos nós no arquivo segue a ordem das entidades.
    # Para que a numeração 1, 2, 3... apareça em ordem no arquivo, 
    # vamos colocar TODOS os nós em uma única entidade discreta de maior dimensão.
    max_dim = 3 if any(d == 3 for d,_,_,_,_,_ in entities_data) else 2
    main_entity = gmsh.model.addDiscreteEntity(max_dim, 9999)

    # Adicionamos TODOS os nós de uma vez na entidade 9999, na ordem correta (1, 2, 3...)
    sorted_tags = [i + 1 for i in range(len(indices))]
    sorted_coords = []
    for idx in indices:
        sorted_coords.extend(coords[idx])
    gmsh.model.mesh.addNodes(max_dim, 9999, sorted_tags, sorted_coords)

    # Agora adicionamos os elementos em suas entidades originais
    for dim, e_tag, e_node_tags, e_types, e_tags, e_conn in entities_data:
        gmsh.model.addDiscreteEntity(dim, e_tag)
        for i in range(len(e_types)):
            new_conn = [old_to_new[int(n)] for n in e_conn[i]]
            gmsh.model.mesh.addElements(dim, e_tag, [e_types[i]], [e_tags[i]], [new_conn])

    # Re-criar Grupos Físicos
    for dim, p_tag, name, ents in physicals:
        gmsh.model.addPhysicalGroup(dim, ents, p_tag)
        if name: gmsh.model.setPhysicalName(dim, p_tag, name)

    # gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    # gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(filename + ".msh2")
    gmsh.finalize()
    print("--- SUCESSO: MALHA REORDENADA E SALVA ---")

# def get_reorder_mesh(filename, meshdata):
#     gmsh.initialize()
#     gmsh.open(filename + ".geo")

#     # 1. MUDANÇA: Gera malha 3D
#     gmsh.model.mesh.generate(abs(int(gmsh_key(meshdata["meshconfig"]["mesh"]))))

#     node_tags, coords, _ = gmsh.model.mesh.getNodes()
#     coords = coords.reshape(-1, 3)
#     elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements()

#     coords_rounded = np.round(coords, 6)
#     # Ordenação léxica já contempla Z, Y e X (Z é o critério primário aqui)
#     indices = np.lexsort((coords_rounded[:, 0], coords_rounded[:, 1], coords_rounded[:, 2]))
#     old_to_new = {int(old): i + 1 for i, old in enumerate(node_tags[indices])}

#     # NOVO MODELO
#     gmsh.model.add("MalhaReordenada")

#     # 2. MUDANÇA: Adicionamos uma entidade para cada dimensão possível (0 a 3)
#     entities = {
#         0: gmsh.model.addDiscreteEntity(0, 1),
#         1: gmsh.model.addDiscreteEntity(1, 1),
#         2: gmsh.model.addDiscreteEntity(2, 1),
#         3: gmsh.model.addDiscreteEntity(3, 1)
#     }

#     # Adicionamos todos os nós na entidade de maior dimensão (3 se for 3D, ou 2 se for 2D)
#     max_dim = max(entities.keys())
#     for idx in indices:
#         new_tag = old_to_new[int(node_tags[idx])]
#         c = coords[idx]
#         gmsh.model.mesh.addNodes(max_dim, entities[max_dim], [new_tag], [c[0], c[1], c[2]])

#     # 3. MUDANÇA: O loop agora direciona para a entidade correta (0, 1, 2 ou 3)
#     for i in range(len(elem_types)):
#         e_type = elem_types[i]
#         dim = gmsh.model.mesh.getElementProperties(e_type)[1]
        
#         target_tag = entities[dim]
#         new_connectivity = [old_to_new[int(node)] for node in elem_node_tags[i]]
#         gmsh.model.mesh.addElements(dim, target_tag, [e_type], [elem_tags[i].astype(int).tolist()], [new_connectivity])

#     # gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
#     # gmsh.option.setNumber("Mesh.Binary", 0)
#     gmsh.write(filename + ".msh2")
#     gmsh.finalize()
#     print("\n--- SUCESSO: MALHA SALVA E REORDENADA ---")


def set_gmsh_geo(filename, meshdata):
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
                file_object.write('Merge "' + meshdata["cadimport"]['object'] + '";\n')
            else:
                npl = 0
                phl = 0
                for i in range(0, numplalistP):
                    npl += 1
                    file_object.write(
                        "Curve Loop("
                        + str(npl)
                        + ") = {"
                        + (str(np.abs(meshdata["planelist"][i][:]).tolist()))[1:-1]
                        + "};\n"
                    )

                npladd = 0
                lplrm = []
                nplrm = 0
                for i in range(0, numplalistP): 
                    if any(jj < 0 for jj in meshdata["planelist"][i][:]):
                        nplrm += 1
                        lplrm.append(i + 1)
                    else:
                        npladd += 1

                addpl = 0
                if nplrm > 0:
                    lplrm.insert(0, lplrm[0] - 1)
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
