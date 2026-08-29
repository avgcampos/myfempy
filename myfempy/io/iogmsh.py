#!/usr/bin/env python
import os
import gmsh
import sys
import numpy as np

__docformat__ = "google"

__doc__ = """

==========================================================================
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
    return celltype[elemid]


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
    cmd = f"gmsh -v 0 {filename}.geo {gmsh_key(meshdata['meshconfig']['mesh'])} -o {filename}.msh2"
    os.system(cmd)


def get_reorder_mesh(filename, meshdata):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(filename + ".geo")

    gmsh.model.mesh.generate(abs(int(gmsh_key(meshdata["meshconfig"]["mesh"]))))

    # 1. Coleta e ordenação vetorial por PESO GEOMÉTRICO (Z -> Y -> X)
    node_tags, coords, _ = gmsh.model.mesh.getNodes()
    coords = coords.reshape(-1, 3)
    c_round = np.round(coords, 6)

    sort_scores = c_round[:, 2] * 1e12 + c_round[:, 1] * 1e6 + c_round[:, 0]
    indices = np.argsort(sort_scores)

    old_to_new = {int(node_tags[idx]): i + 1 for i, idx in enumerate(indices)}

    # 2. Coleta estruturada de entidades e grupos físicos
    physicals = [
        (dim, p_tag, gmsh.model.getPhysicalName(dim, p_tag), gmsh.model.getEntitiesForPhysicalGroup(dim, p_tag))
        for dim, p_tag in gmsh.model.getPhysicalGroups()
    ]

    entities_data = []
    for dim in range(4):
        for _, e_tag in gmsh.model.getEntities(dim):
            e_node_tags, _, _ = gmsh.model.mesh.getNodes(dim, e_tag)
            e_types, e_tags, e_conn = gmsh.model.mesh.getElements(dim, e_tag)
            entities_data.append((dim, e_tag, e_node_tags, e_types, e_tags, e_conn))

    # 3. Reconstrução Otimizada no Gmsh
    gmsh.model.add("MalhaReordenada")

    max_dim = 3 if any(d == 3 for d, _, _, _, _, _ in entities_data) else 2
    gmsh.model.addDiscreteEntity(max_dim, 9999)

    sorted_tags = np.arange(1, len(indices) + 1).tolist()
    sorted_coords = coords[indices].flatten().tolist()
    gmsh.model.mesh.addNodes(max_dim, 9999, sorted_tags, sorted_coords)

    for dim, e_tag, _, e_types, e_tags, e_conn in entities_data:
        gmsh.model.addDiscreteEntity(dim, e_tag)
        for i in range(len(e_types)):
            new_conn = [old_to_new[int(n)] for n in e_conn[i]]
            gmsh.model.mesh.addElements(dim, e_tag, [e_types[i]], [e_tags[i]], [new_conn])

    for dim, p_tag, name, ents in physicals:
        gmsh.model.addPhysicalGroup(dim, ents, p_tag)
        if name:
            gmsh.model.setPhysicalName(dim, p_tag, name)

    gmsh.write(filename + ".msh2")
    gmsh.finalize()
    print("--- SUCESSO: MALHA REORDENADA E SALVA ---")


def set_gmsh_geo(filename, meshdata):
    lines = ["// GMSH GEOMETRY FILE FROM MYFEMPY", 'SetFactory("OpenCASCADE");']
    
    has_points = "pointlist" in meshdata.keys()
    if has_points:
        numlinlist = len(meshdata.get("linelist", []))
        line_list = ",".join(str(i + 1) for i in range(numlinlist))
        
        for i, pt in enumerate(meshdata["pointlist"]):
            lines.append(f"Point({i + 1}) = {{{pt[0]}, {pt[1]}, {pt[2]}}};")

        if "circle" in meshdata.keys():
            for inl, circ in enumerate(meshdata["circle"]):
                d, (cx, cy, cz), (arc0, arc1) = circ[0], circ[1], circ[2]
                lines.append(f"Circle({numlinlist + inl + 1}) = {{{cx}, {cy}, {cz}, {d}, {arc0}, {arc1}}};")

        if "arc" in meshdata.keys():
            numincl = len(meshdata.get("circle", []))
            for iarc, arc in enumerate(meshdata["arc"]):
                lines.append(f"Circle({numlinlist + numincl + iarc + 1}) = {{{arc[0]}, {arc[1]}, {arc[2]}}};")

        for i, ln in enumerate(meshdata.get("linelist", [])):
            lines.append(f"Line({i + 1}) = {{{ln[0]}, {ln[1]}}};")

    mesh_type = meshdata["meshconfig"]["mesh"]

    if mesh_type in ["line2", "line3"]:
        if "numbernodes" in meshdata["meshconfig"]:
            lines.append(f"Transfinite Curve {{{line_list}}} = {meshdata['meshconfig']['numbernodes']} Using Progression 1;")
        elif "sizeelement" in meshdata["meshconfig"]:
            for i, pt in enumerate(meshdata["pointlist"]):
                lines[i + 2] = f"Point({i + 1}) = {{{pt[0]}, {pt[1]}, {pt[2]}, {meshdata['meshconfig']['sizeelement']}}};"
        
        lines.extend([
            "// MESH CONFIGURATION",
            "Mesh.CharacteristicLengthExtendFromBoundary = 1;",
            "Mesh.CharacteristicLengthMin = 0;",
            "Mesh.CharacteristicLengthFromPoints = 1;",
            "Mesh.Optimize = 1;",
            "Mesh.HighOrderOptimize = 0;",
            "Mesh.Algorithm = 8;"
        ])
        if mesh_type == "line2":
            lines.append("Mesh.ElementOrder = 1;")
        else:
            lines.extend(["Mesh.SecondOrderIncomplete = 1;", "Mesh.ElementOrder = 2;"])
    else:
        if "cadimport" in meshdata.keys():
            lines.append(f'Merge "{meshdata["cadimport"]["object"]}";')
        else:
            planelist = meshdata.get("planelist", [])
            numplalistP = len(planelist)
            for i, plane in enumerate(planelist):
                pl_str = str(np.abs(plane).tolist())[1:-1]
                lines.append(f"Curve Loop({i + 1}) = {{{pl_str}}};")

            lplrm = []
            npladd = 0
            for i, plane in enumerate(planelist):
                if any(jj < 0 for jj in plane):
                    lplrm.append(i + 1)
                else:
                    npladd += 1

            addpl = 0
            if lplrm:
                lplrm.insert(0, lplrm[0] - 1)
                addpl = 1
                lines.append(f"Plane Surface({addpl}) = {{{str(list(set(lplrm)))[1:-1]}}};")
            
            for iap in range(addpl, npladd):
                lines.append(f"Plane Surface({iap + 1}) = {{{addpl + iap + 1}}};")

            lines.append(f"Characteristic Length {{:}} = {meshdata['meshconfig']['sizeelement']};")

        if "meshmap" in meshdata["meshconfig"] and meshdata["meshconfig"]["meshmap"].get("on"):
                    lines.append("//FACE MAPPING")
                    mmap = meshdata["meshconfig"]["meshmap"]
                    if "numbernodes" in mmap:
                        if mmap["edge"] == "all":
                            lines.append(f"Transfinite Curve {{:}} = {mmap['numbernodes']} Using Progression 1;")
                        else:
                            for ed, edge_val in enumerate(mmap["edge"]):
                                # Mantém a formatação original removendo os colchetes das extremidades
                                edge_str = str(edge_val)[1:-1]
                                lines.append(f"Transfinite Curve {{{edge_str}}} = {mmap['numbernodes'][ed]} Using Progression 1;")
                    elif mmap["edge"] == "all":
                        lines.append("Transfinite Surface {:};")
                    else:
                        lines.append("Transfinite Surface {:};")

        lines.append(f"// MESH {mesh_type} CONFIGURATION")
        sz = meshdata["meshconfig"]["sizeelement"]

        if mesh_type in ["tria3", "tria6"]:
            lines.extend([
                "Mesh.CharacteristicLengthExtendFromBoundary = 1;",
                "Mesh.CharacteristicLengthMin = 0;",
                f"Mesh.CharacteristicLengthMax = {sz};",
                "Mesh.CharacteristicLengthFromPoints = 1;",
                "Mesh.Optimize = 1;",
                "Mesh.HighOrderOptimize = 0;",
                "Mesh.Algorithm = 8;"
            ])
            if mesh_type == "tria3":
                lines.append("Mesh.ElementOrder = 1;")
            else:
                lines.extend(["Mesh.SecondOrderIncomplete = 1;", "Mesh.ElementOrder = 2;"])

        elif mesh_type in ["quad4", "quad8"]:
            lines.extend([
                "Recombine Surface {:};",
                "Mesh.RecombinationAlgorithm = 1;",
                "Mesh.RecombineAll = 1;",
                "Mesh.SubdivisionAlgorithm = 1;",
                "Mesh.CharacteristicLengthExtendFromBoundary = 1;",
                "Mesh.CharacteristicLengthMin = 0;",
                f"Mesh.CharacteristicLengthMax = {sz};",
                "Mesh.CharacteristicLengthFromPoints = 1;",
                "Mesh.Optimize = 1;",
                "Mesh.HighOrderOptimize = 0;",
                "Mesh.Algorithm = 8;"
            ])
            if mesh_type == "quad4":
                lines.append("Mesh.ElementOrder = 1;")
            else:
                lines.extend(["Mesh.SecondOrderIncomplete = 1;", "Mesh.ElementOrder = 2;"])

        elif mesh_type == "tetr4":
            if "extrude" in meshdata["meshconfig"]:
                thck = meshdata["meshconfig"]["extrude"]
                lines.append(f"Extrude {{0, 0, {float(thck)}}} {{Surface{{:}};}}")
            lines.extend([
                "Mesh.Algorithm = 2;",
                "Mesh.Algorithm3D = 4;",
                "Mesh.CharacteristicLengthExtendFromBoundary = 1;",
                "Mesh.CharacteristicLengthMin = 0;",
                f"Mesh.CharacteristicLengthMax = {sz};",
                "Mesh.ElementOrder = 1;",
                "Mesh.Optimize = 1;",
                "Mesh.HighOrderOptimize = 0;"
            ])

        elif mesh_type == "hexa8":
            if "extrude" in meshdata["meshconfig"]:
                thck = float(meshdata["meshconfig"]["extrude"])
                layers = int(thck / float(sz))
                lines.append(f"Extrude {{0, 0, {thck}}} {{Surface{{:}};Layers{{{layers}}};Recombine;}};")
            lines.extend([
                "Recombine Surface {:};",
                "Mesh.Algorithm = 2;",
                "Mesh.Algorithm3D = 4;",
                "Mesh.CharacteristicLengthExtendFromBoundary = 1;",
                "Mesh.CharacteristicLengthMin = 0;",
                f"Mesh.CharacteristicLengthMax = {sz};",
                "Mesh.ElementOrder = 1;",
                "Mesh.Optimize = 1;",
                "Mesh.HighOrderOptimize = 0;",
                "Mesh.RecombinationAlgorithm = 0;",
                "Mesh.SubdivisionAlgorithm = 2;",
                "Mesh.RecombineAll = 1;"
            ])

    with open(filename + ".geo", "w") as file_object:
        file_object.write("\n".join(lines) + "\n")