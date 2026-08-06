#!/usr/bin/env python
__doc__ = """
Plotter Prev Process
"""

import numpy as np
import vedo

from myfempy.io.iovtk import convert_to_vtk
from myfempy.plots.physics import (view_beam_crossSection, view_bondcond_point,
                                   view_listforce, view_text_point)
from myfempy.utils.utils import get_version

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

# versao atualizada/ otimizada
def preview_plot(Model, previewset, path, Physic=None):
    """
    Prepara os dados do modelo e dispara a construção do preview em VTK.
    """
    render = previewset.setdefault("RENDER", {})
    
    # 1. Configuração básica do plotdata para exportação VTK
    plotdata = {
        "coord": Model.coord,
        "inci": Model.inci,
        "nodecon": Model.modelinfo["nodecon"],
        "elemid": Model.modelinfo["elemid"],
        "filename": f"{path}/{render.get('filename', 'model_preview')}",
        "title": ["UNDEFORM_MESH"],
        "solution": np.ones((len(Model.inci), 1)),
        "material_CELL_DATA_val": np.array([Model.inci[:, 2]]).T.astype(int),
        "material_CELL_DATA_title": ["Material_Set"]
    }

    convert_to_vtk(plotdata)

    # 2. Tratamento do fator de escala
    scale_val = render.get("scale", 100)
    max_coord = max(np.max(np.abs(Model.coord[:, 1:4])) , 1e-6)
    render["scale"] = (scale_val / 100) * max_coord

    # 3. Padrões de renderização (linhas, cs, tags)
    render.setdefault("lines", False)
    render.setdefault("cs", False)
    render.setdefault("numb", False)
    
    # Processamento de plottags e regions de forma limpa
    plottags = render.get("plottags", {})
    region_map = {"point": 0, "line": 1, "plane": 2}
    previewset["regions"] = [[], []]
    for tag, idx in region_map.items():
        if plottags.get(tag, False):
            previewset["regions"] = Model.regions[idx]
            break

    # 4. Seções transversais (tabcs)
    previewset["tabcs"] = {
        "typSection": [int(g["ID"]) for g in Model.tabgeo],
        "dimSection": [[g["B"], g["H"], g["T"], g["D"]] for g in Model.tabgeo]
    }

    # 5. Física (Forças e Restrições) com tratamento seguro
    if Physic is not None:
        if hasattr(Physic, "forces"):
            previewset["forces"] = Physic.forces
        if hasattr(Physic, "constrains"):
            previewset["constrains"] = Physic.constrains

    # 6. Repasse direto de atributos estruturais essenciais
    previewset.update({
        "coord": Model.coord,
        "inci": Model.inci,
        "nnode": Model.modelinfo["nnode"],
        "nodecon": Model.modelinfo["nodecon"],
        "dofs": Model.modelinfo["dofs"]
    })

    build_preview(previewset, path)


def build_preview(previewset: dict, path):
    """
    Carrega a malha gerada e monta a cena gráfica utilizando o Vedo.
    """
    render = previewset["RENDER"]
    file_name = f"{path}/{render['filename']}.vtk"
    
    mesh = vedo.load(file_name)
    if not mesh:
        raise FileNotFoundError(f"Não foi possível carregar o arquivo: {file_name}")
    
    mesh.cmap("jet")
    if render.get("lines", True):
        mesh.lw(1.0).linecolor("black")
    else:
        mesh.lw(0)

    actors = [mesh]
    scala_view = render["scale"]
    coord = previewset["coord"]

    # 1. Forças
    if "forces" in previewset:
        dof_f = previewset["dofs"]["f"]
        key_list_fc, val_list_fc = list(dof_f.keys()), list(dof_f.values())
        
        for frcApy_vet in previewset["forces"]:
            fc_copy = frcApy_vet.copy()
            fc_copy[1] = __setLoadDof(key_list_fc[val_list_fc.index(fc_copy[1])])
            cone1, cone2 = view_listforce(coord, fc_copy, scala_view)
            actors.extend([cone1, cone2])

    if "constrains" in previewset.keys():
        key_list_bc = list(previewset["dofs"]["d"].keys())
        val_list_bc = list(previewset["dofs"]["d"].values())

        dimbclist = previewset["constrains"].shape[0]
        for num_bc in range(dimbclist):
            bondCond_vet = previewset["constrains"][[num_bc]][0]
            if int(bondCond_vet[1]) == 0:
                pass
            
            elif int(bondCond_vet[1]) == 11 or int(bondCond_vet[1]) == 12:
                bondCond_vet[1] = 4

            elif int(bondCond_vet[1]) == 13 or int(bondCond_vet[1]) == 14:
                bondCond_vet[1] = 1

            elif int(bondCond_vet[1]) == 15 or int(bondCond_vet[1]) == 16:
                bondCond_vet[1] = 2

            elif int(bondCond_vet[1]) == 17 or int(bondCond_vet[1]) == 18 or int(bondCond_vet[1]) == 19 or int(bondCond_vet[1]) == 20:
                bondCond_vet[1] = 0

            else:
                bondCond_vet[1] = __setBCDof(
                    key_list_bc[val_list_bc.index(int(bondCond_vet[1]))]
                )
                
            pt_actor, tdof_actor = view_bondcond_point(previewset["coord"], bondCond_vet, scala_view)
            actors.extend([pt_actor, tdof_actor])

    # 3. Tags de Regiões / Textos
    if render.get("plottags"):
        coordMax = np.max(coord[:, 1:4], axis=0)
        reg_title, reg_items = previewset["regions"]
        
        for reg_item in reg_items:
            node_indices = reg_item[1] - 1
            mean_coord = np.mean(coord[node_indices, 1:4], axis=0)
            txt_actor = view_text_point(mean_coord, coordMax, scala_view, [reg_title, str(reg_item[0])])
            actors.append(txt_actor)

    # 4. Seções Transversais de Vigas (CS)
    if render.get("cs", False):
        inci = previewset["inci"]
        for num_bcs, (typ_sec, dim_sec) in enumerate(zip(previewset["tabcs"]["typSection"], previewset["tabcs"]["dimSection"])):
            inci_bcs = inci[inci[:, 3] == num_bcs + 1]
            for row in inci_bcs:
                noi, noj = int(row[4]) - 1, int(row[5]) - 1
                coord_bcs = np.concatenate([coord[noi, 1:4], coord[noj, 1:4]])
                beam_actor = view_beam_crossSection(dim_sec, int(typ_sec), coord_bcs, scala_view)
                actors.append(beam_actor)

    # 5. Numeração de Nós e Elementos (numb)
    if render.get("numb", True):
        node_points = coord[:, 1:4]
        actors.append(vedo.Points(node_points, c="yellow", r=5))

        for idx, pt in enumerate(node_points):
            node_id_str = str(int(coord[idx, 0])) if coord.shape[1] > 3 else str(idx + 1)
            actors.append(vedo.Text3D(node_id_str, pos=pt, s=scala_view * 0.6, c="yellow"))

        if "inci" in previewset and previewset["inci"] is not None:
            for idx, row in enumerate(previewset["inci"]):
                nodes_in_elem = [int(row[4]) - 1, int(row[5]) - 1] if len(row) > 5 else []
                if nodes_in_elem and all(n < len(coord) for n in nodes_in_elem):
                    centroid = np.mean(coord[nodes_in_elem, 1:4], axis=0)
                    actors.append(vedo.Text3D(str(idx + 1), pos=centroid, s=scala_view * 0.8, c="cyan"))

    # 6. Plotter e Renderização Final
    plt = vedo.Plotter(size=(640, 480), bg='white', bg2='lightblue')
    
    version_str = get_version() if 'get_version' in globals() else ""
    actors.append(vedo.Text2D(f"MYFEMPY " + get_version() + 
                            '\n> press "w" to wireframe/surface view \n> press "r" to reset view \n> press "q" to exit and continue'
                              , pos="top-left", s=1, c="black"))
    
    plt.show(actors, axes=4, interactive=False, viewup='y', zoom=1.4)

    if render.get("savepng", False):
        plt.screenshot(f"{path}/{render['filename']}.png")

    if render.get("show", True):
        plt.interactive().close()


def __setLoadDof(forcedof):
    fdoftype = {
        "fx": 1,
        "fy": 2,
        "fz": 3,
        "tx": 4,
        "ty": 5,
        "tz": 6,
        "masspoint": 15,
        "spring2ground": 16,
        "damper2ground": 17,
        "torsion2ground": 17,
        "cg": 18,
        "heatflux": 1,
        "convectionfluid": 2,
        "heat2fluid": 15,
    }
    return fdoftype[forcedof]


def __setBCDof(bcdof):
    bcdoftype = {
        "ux": 1,
        "uy": 2,
        "uz": 3,
        "rx": 4,
        "ry": 5,
        "rz": 6,
        "full": 0,
        "t": 1,
    }
    return bcdoftype[bcdof]
