#!/usr/bin/env python
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

# Tabela estática fora da função para acesso mais rápido
_CELL_TYPE_MAP = {
    1121: 3, 1132: 21, 1621: 3, 1632: 21,
    2131: 5, 2162: 22, 2141: 9, 2182: 23,
    2231: 5, 2262: 22, 2241: 9, 2282: 23,
    3141: 10, 3341: 10, 33102: 24, 3181: 12,
    3381: 12, 33202: 25,
}

def meshid2vtkid(elemid):
    """Retorna o ID do VTK baseado no ID da malha."""
    return _CELL_TYPE_MAP.get(int(elemid), 3)

def convert_to_vtk(plotdata):
    """Converte dados do myfempy para o formato VTK de forma vetorizada."""
    numnodes = len(plotdata["coord"])
    numelem = len(plotdata["inci"])
    
    # Pré-calcula tamanho da lista de células
    inci_subset = plotdata["inci"][:, 4:]
    cell_size = (plotdata["nodecon"] + 1) * numelem

    with open(plotdata["filename"] + ".vtk", "w") as file_object:
        # Cabeçalho
        header_lines = [
            "# vtk DataFile Version 4.0",
            "vtk output from myfempy solver",
            "ASCII",
            "DATASET UNSTRUCTURED_GRID",
            f"POINTS {numnodes} double"
        ]
        file_object.write("\n".join(header_lines) + "\n")
        
        # Coordenadas de forma vetorizada
        coord_str = "\n".join(
            " ".join(row) for row in plotdata["coord"][:, 1:].astype(str)
        )
        file_object.write(coord_str + "\n\n")
        
        # Células
        file_object.write(f"CELLS {numelem} {cell_size}\n")
        cells_lines = []
        for ii in range(numelem):
            listnodes = inci_subset[ii]
            valid_nodes = listnodes[listnodes.nonzero()]
            list2write = (valid_nodes - 1).astype(int).astype(str)
            cells_lines.append(f"{len(list2write)} " + " ".join(list2write))
        file_object.write("\n".join(cells_lines) + "\n\n")
        
        # Tipos de Células
        file_object.write(f"CELL_TYPES {numelem}\n")
        types_lines = [str(meshid2vtkid(plotdata["inci"][ii, 1])) for ii in range(numelem)]
        file_object.write("\n".join(types_lines) + "\n\n")
        
        # -----------------------------------------------------
        # POINT DATA
        file_object.write(f"POINT_DATA {numnodes}\n")

        if "displ_POINT_DATA_val" in plotdata:
            val = plotdata["displ_POINT_DATA_val"]
            file_object.write("FIELD FieldData 1\n")
            file_object.write(f"{plotdata['displ_POINT_DATA_title']} 3 {len(val)} float\n")
            disp_str = "\n".join(" ".join(row) for row in val.astype(str))
            file_object.write(disp_str + "\n")

        if "modes_POINT_DATA" in plotdata:
            for md_data in plotdata["modes_POINT_DATA"]:
                val = md_data["val"]
                file_object.write("\nFIELD FieldData 1\n")
                file_object.write(f"{md_data['title']} 3 {len(val)} float\n")
                modes_str = "\n".join(" ".join(row) for row in val.astype(str))
                file_object.write(modes_str + "\n")

        if "frf_POINT_DATA" in plotdata:
            for frf_data in plotdata["frf_POINT_DATA"]:
                val = frf_data["val"]
                file_object.write("\nFIELD FieldData 1\n")
                file_object.write(f"{frf_data['title']} 3 {len(val)} float\n")
                frf_str = "\n".join(" ".join(row) for row in val.astype(str))
                file_object.write(frf_str + "\n")

        if "temp_POINT_DATA_val" in plotdata:
            val = plotdata["temp_POINT_DATA_val"]
            file_object.write(f"SCALARS {plotdata['temp_POINT_DATA_title']} float 1\n")
            file_object.write("LOOKUP_TABLE default\n")
            temp_str = "\n".join(" ".join(row) for row in val.astype(str))
            file_object.write(temp_str + "\n")

        file_object.write("\n")
        
        # -----------------------------------------------------
        # CELL DATA
        file_object.write(f"CELL_DATA {numelem}\n")
        
        if "material_CELL_DATA_val" in plotdata:
            mat_val = plotdata["material_CELL_DATA_val"]
            for jj, title in enumerate(plotdata["material_CELL_DATA_title"]):
                file_object.write(f"SCALARS {title} float 1\n")
                file_object.write("LOOKUP_TABLE default\n")
                file_object.write("\n".join(mat_val[:, jj].astype(str)) + "\n")

        if "stress_CELL_DATA_val" in plotdata:
            stress_val = plotdata["stress_CELL_DATA_val"]
            for jj, title in enumerate(plotdata["stress_CELL_DATA_title"]):
                file_object.write(f"SCALARS {title} float 1\n")
                file_object.write("LOOKUP_TABLE default\n")
                file_object.write("\n".join(stress_val[:, jj].astype(str)) + "\n")


def convert_from_vtk(filename):
    """Lê o arquivo VTK otimizando a leitura das linhas com list comprehensions."""
    file_imp = filename + ".vtk"
    with open(file_imp, "r") as file_object:
        # Pula as 4 primeiras linhas do cabeçalho
        for _ in range(4):
            file_object.readline()
            
        # Lê número de nós
        lineaux = file_object.readline().split()
        nnod = int(lineaux[1])
        
        nodelist = []
        for ii in range(nnod):
            lineaux = file_object.readline().split()
            nodelist.append([ii + 1, float(lineaux[0]), float(lineaux[1]), float(lineaux[2])])
            
        file_object.readline() # Linha em branco
        
        # Lê número de elementos
        lineaux = file_object.readline().split()
        nelm = int(lineaux[1])
        
        conec_elm = [
            [float(x) for x in file_object.readline().split()] 
            for _ in range(nelm)
        ]
        
    return conec_elm, nodelist