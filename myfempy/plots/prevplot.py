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
    plt = vedo.Plotter(size=(1280, 720), bg='black', bg2='lightblue')
    
    version_str = get_version() if 'get_version' in globals() else ""
    actors.append(vedo.Text2D(f"MYFEMPY {version_str}", pos="top-left", s=1, c="black"))
    
    plt.show(actors, axes=4, interactive=False, viewup='y', zoom=1.4)

    if render.get("savepng", False):
        plt.screenshot(f"{path}/{render['filename']}.png")

    if render.get("show", True):
        plt.interactive().close()




# def preview_plot(Model, previewset, path, Physic = None):
#     # path = os.getcwd()
#     plotdata = dict()
#     plotdata["coord"] = Model.coord
#     plotdata["inci"] = Model.inci
#     plotdata["nodecon"] = Model.modelinfo["nodecon"]
#     plotdata["elemid"] = Model.modelinfo["elemid"]
#     plotdata["filename"] = path + "/" + previewset["RENDER"]["filename"]
#     plotdata["title"] = ["UNDEFORM_MESH"]
#     plotdata["solution"] = np.ones((len(Model.inci), 1))
#     plotdata["material_CELL_DATA_val"] = (
#         (np.array([(Model.inci[:, 2])])).T
#     ).astype(int)
#     plotdata["material_CELL_DATA_title"] = ["Material_Set"]

#     convert_to_vtk(plotdata)

#     if "scale" in previewset["RENDER"].keys():
#         previewset["RENDER"]["scale"] = (previewset["RENDER"]["scale"] / 100) * max(
#             [
#                 max(abs(Model.coord[:, 1])),
#                 max(abs(Model.coord[:, 2])),
#                 max(abs(Model.coord[:, 3])),
#             ]
#         )
#     else:
#         previewset["RENDER"]["scale"] = 1

#     if "lines" in previewset["RENDER"].keys():
#         pass
#     else:
#         previewset["RENDER"]["lines"] = True

#     if "plottags" in previewset["RENDER"].keys():
#         if (
#             "point" in previewset["RENDER"]["plottags"].keys()
#             and previewset["RENDER"]["plottags"]["point"] == True
#         ):
#             previewset["regions"] = Model.regions[0]
#         else:
#             pass

#         if (
#             "line" in previewset["RENDER"]["plottags"].keys()
#             and previewset["RENDER"]["plottags"]["line"] == True
#         ):
#             previewset["regions"] = Model.regions[1]
#         else:
#             pass

#         if (
#             "plane" in previewset["RENDER"]["plottags"].keys()
#             and previewset["RENDER"]["plottags"]["plane"] == True
#         ):
#             previewset["regions"] = Model.regions[2]
#         else:
#             pass
#     else:
#         previewset["regions"] = [[], []]

#     if "cs" in previewset["RENDER"].keys():
#         previewset["RENDER"]["cs"] = True
#     else:
#         previewset["RENDER"]["cs"] = False

#     previewset["tabcs"] = dict()
#     previewset["tabcs"]["typSection"] = []
#     previewset["tabcs"]["dimSection"] = []
#     for gg in range(len(Model.tabgeo)):
#         previewset["tabcs"]["typSection"].append(int(Model.tabgeo[gg]["ID"]))
#         # previewset["tabcs"]["dimSection"].append(modelinfo["tabgeo"][gg][5:9])
#         # dim = [modelinfo["tabgeo"][gg]['B'], modelinfo["tabgeo"][gg]['B'], modelinfo["tabgeo"][gg]['B'], modelinfo["tabgeo"][gg]['B']]
#         previewset["tabcs"]["dimSection"].append(
#             [
#                 Model.tabgeo[gg]["B"],
#                 Model.tabgeo[gg]["H"],
#                 Model.tabgeo[gg]["T"],
#                 Model.tabgeo[gg]["D"],
#             ]
#         )

#     try:
#         previewset["forces"] = Physic.forces
#     except:
#         pass

#     try:
#         previewset["constrains"] = Physic.constrains
#     except:
#         pass

#     previewset["coord"] = Model.coord
#     previewset["inci"] = Model.inci
#     previewset["nnode"] = Model.modelinfo["nnode"]
#     previewset["nodecon"] = Model.modelinfo["nodecon"]
#     previewset["dofs"] = Model.modelinfo["dofs"]

#     build_preview(previewset, path)


# def build_preview(previewset: dict, path):
#     file_name = f"{path}/{previewset['RENDER']['filename']}.vtk"
    
#     # 1. Carregar a malha e configurar o mapeamento de cores com vedo
#     mesh = vedo.load(file_name)
#     if not mesh:
#         raise FileNotFoundError(f"Não foi possível carregar o arquivo: {file_name}")
    
#     # Configurar LUT (Jet) e range de escalares
#     mesh.cmap("jet")
    
#     # Visibilidade de arestas / linhas
#     show_edges = previewset["RENDER"].get("lines", True)
#     if show_edges:
#         mesh.lw(1.0).linecolor("black")
#     else:
#         mesh.lw(0)

#     # Lista para acumular todos os atores que irão para a cena
#     actors = [mesh]

#     scala_view = previewset["RENDER"]["scale"]

#     # 2. Forças
#     if "forces" in previewset:
#         key_list_fc = list(previewset["dofs"]["f"].keys())
#         val_list_fc = list(previewset["dofs"]["f"].values())
        
#         for num_lf in range(previewset["forces"].shape[0]):
#             frcApy_vet = previewset["forces"][num_lf].copy()
#             frcApy_vet[1] = __setLoadDof(key_list_fc[val_list_fc.index(frcApy_vet[1])])
            
#             cone1, cone2 = view_listforce(previewset["coord"], frcApy_vet, scala_view)
#             actors.extend([cone1, cone2])

#     # 3. Condições de Contorno (Constrains)
#     if "constrains" in previewset.keys():
#         key_list_bc = list(previewset["dofs"]["d"].keys())
#         val_list_bc = list(previewset["dofs"]["d"].values())

#         dimbclist = previewset["constrains"].shape[0]
#         for num_bc in range(dimbclist):
#             bondCond_vet = previewset["constrains"][[num_bc]][0]
#             if int(bondCond_vet[1]) == 0:
#                 pass
            
#             elif int(bondCond_vet[1]) == 11 or int(bondCond_vet[1]) == 12:
#                 bondCond_vet[1] = 4

#             elif int(bondCond_vet[1]) == 13 or int(bondCond_vet[1]) == 14:
#                 bondCond_vet[1] = 1

#             elif int(bondCond_vet[1]) == 15 or int(bondCond_vet[1]) == 16:
#                 bondCond_vet[1] = 2

#             elif int(bondCond_vet[1]) == 17 or int(bondCond_vet[1]) == 18 or int(bondCond_vet[1]) == 19 or int(bondCond_vet[1]) == 20:
#                 bondCond_vet[1] = 0

#             else:
#                 bondCond_vet[1] = __setBCDof(
#                     key_list_bc[val_list_bc.index(int(bondCond_vet[1]))]
#                 )
                
#             pt_actor, tdof_actor = view_bondcond_point(previewset["coord"], bondCond_vet, scala_view)
#             actors.extend([pt_actor, tdof_actor])

#     # 4. Tags de Regiões / Textos
#     if previewset["RENDER"].get("plottags", False):
#         coordMax = [
#             max(previewset["coord"][:, 1]),
#             max(previewset["coord"][:, 2]),
#             max(previewset["coord"][:, 3]),
#         ]
        
#         for reg_item in previewset["regions"][1]:
#             node_indices = reg_item[1] - 1
#             text = [previewset["regions"][0], str(reg_item[0])]
#             coord = [
#                 np.mean(previewset["coord"][node_indices, 1]),
#                 np.mean(previewset["coord"][node_indices, 2]),
#                 np.mean(previewset["coord"][node_indices, 3]),
#             ]
#             txt_actor = view_text_point(coord, coordMax, scala_view, text)
#             actors.append(txt_actor)

#     # 5. Seções Transversais de Vigas (CS)
#     if previewset["RENDER"].get("cs", True):
#         for num_bcs, (typ_sec, dim_sec) in enumerate(zip(previewset["tabcs"]["typSection"], previewset["tabcs"]["dimSection"])):
#             inci_bcs = previewset["inci"][previewset["inci"][:, 3] == num_bcs + 1]
#             typSection = int(typ_sec)
            
#             for row in inci_bcs:
#                 noi, noj = int(row[4]), int(row[5])
#                 coord_bcs = np.array([
#                     previewset["coord"][noi - 1, 1], previewset["coord"][noi - 1, 2], previewset["coord"][noi - 1, 3],
#                     previewset["coord"][noj - 1, 1], previewset["coord"][noj - 1, 2], previewset["coord"][noj - 1, 3]
#                 ])
#                 beam_actor = view_beam_crossSection(dim_sec, typSection, coord_bcs, scala_view)
#                 actors.append(beam_actor)

#     # 6 Subrotina para Numeração de Nós e Elementos (numb)
#     if previewset["RENDER"].get("numb", True):
#             node_points = previewset["coord"][:, 1:4]
            
#             # Plota os pontos dos nós
#             pts_actor = vedo.Points(node_points, c="yellow", r=5)
#             actors.append(pts_actor)

#             # Adiciona a numeração individual de cada nó usando Text3D
#             for idx, pt in enumerate(node_points):
#                 node_id_str = str(int(previewset["coord"][idx, 0])) if previewset["coord"].shape[1] > 3 else str(idx + 1)
#                 node_label = vedo.Text3D(node_id_str, pos=pt, s=scala_view * 0.6, c="yellow")
#                 actors.append(node_label)

#             # Numeração de elementos nos centróides
#             if "inci" in previewset and previewset["inci"] is not None:
#                 for idx, row in enumerate(previewset["inci"]):
#                     nodes_in_elem = [int(row[4])-1, int(row[5])-1] if len(row) > 5 else []
#                     if nodes_in_elem and all(n < len(previewset["coord"]) for n in nodes_in_elem):
#                         centroid = np.mean(previewset["coord"][nodes_in_elem, 1:4], axis=0)
#                         elem_label = vedo.Text3D(str(idx + 1), pos=centroid, s=scala_view * 0.8, c="cyan")
#                         actors.append(elem_label)

#     # 7 Configuração do Plotter e Texto de Instruções
#     # plt = vedo.Plotter(size=(640, 480), bg="black")
#     plt = vedo.Plotter(size=(1280, 720), bg='black', bg2='lightblue')
    
#     version_str = get_version() if 'get_version' in globals() else ""
#     help_text = (
#         f"MYFEMPY {version_str}\n"
#         # '> press "w" to wireframe view \n'
#         # '> press "s" to surface view \n'
#         # '> press "r" to reset view \n'
#         # '> press "q" to exit and continue'
#     )
    
#     # Texto 2D fixo no canto superior esquerdo usando âncora padrão do vedo
#     text_actor = vedo.Text2D(help_text, pos="top-left", s=1, c="black")
#     actors.append(text_actor)

#     # Exibe a cena com os atores acumulados e eixos ativados (axes=7)
#     plt.show(actors, axes=4, interactive=False)

#     # 8 Salvar PNG se solicitado
#     if previewset["RENDER"].get("savepng", False):
#         png_path = f"{path}/{previewset['RENDER']['filename']}.png"
#         plt.screenshot(png_path)

#     # 9 Exibir janela interativa se solicitado
#     if previewset["RENDER"].get("show", True):
#         plt.interactive().close()


# versao antiga
# import vtk
# def build_preview(previewset: dict, path):
#     # path = os.getcwd()
#     file_name = str(path + "/" + previewset["RENDER"]["filename"] + ".vtk")
#     renderer = vtk.vtkRenderer()
#     renderer_window = vtk.vtkRenderWindow()
#     renderer_window.SetSize(640, 480)
#     renderer.SetBackground(0.0, 0.0, 0.0)
#     reader = vtk.vtkUnstructuredGridReader()
#     reader.SetFileName(file_name)
#     reader.Update()  # Needed because of GetScalarRange
#     output = reader.GetOutput()
#     output_port = reader.GetOutputPort()
#     scalar_range = output.GetScalarRange()
#     lut = vtk.vtkLookupTable()
#     lut.SetNumberOfColors(8)
#     lut.SetHueRange(0.6, 1)  # jet color
#     lut.Build()
#     mapper = vtk.vtkDataSetMapper()
#     mapper.SetLookupTable(lut)
#     mapper.SetInputConnection(output_port)
#     mapper.SetScalarRange(scalar_range)
#     scala_view = previewset["RENDER"]["scale"]
#     dimfrlist = 0
#     if "forces" in previewset.keys():
#         key_list_fc = list(previewset["dofs"]["f"].keys())
#         val_list_fc = list(previewset["dofs"]["f"].values())
#         dimfrlist = previewset["forces"].shape[0]
#         for num_lf in range(dimfrlist):
#             frcApy_vet = previewset["forces"][[num_lf]][0]
#             frcApy_vet[1] = __setLoadDof(key_list_fc[val_list_fc.index(frcApy_vet[1])])
#             exec(
#                 f'fr_point_actor_cone1_{num_lf},fr_point_actor_cone2_{num_lf} = view_listforce(previewset["coord"],frcApy_vet,scala_view)'
#             )
#     dimbclist = 0
#     if "constrains" in previewset.keys():
#         key_list_bc = list(previewset["dofs"]["d"].keys())
#         val_list_bc = list(previewset["dofs"]["d"].values())

#         dimbclist = previewset["constrains"].shape[0]
#         for num_bc in range(dimbclist):
#             bondCond_vet = previewset["constrains"][[num_bc]][0]
#             if int(bondCond_vet[1]) == 0:
#                 pass
            
#             elif int(bondCond_vet[1]) == 11 or int(bondCond_vet[1]) == 12:
#                 bondCond_vet[1] = 4

#             elif int(bondCond_vet[1]) == 13 or int(bondCond_vet[1]) == 14:
#                 bondCond_vet[1] = 1

#             elif int(bondCond_vet[1]) == 15 or int(bondCond_vet[1]) == 16:
#                 bondCond_vet[1] = 2

#             elif int(bondCond_vet[1]) == 17 or int(bondCond_vet[1]) == 18 or int(bondCond_vet[1]) == 19 or int(bondCond_vet[1]) == 20:
#                 bondCond_vet[1] = 0

#             else:
#                 bondCond_vet[1] = __setBCDof(
#                     key_list_bc[val_list_bc.index(int(bondCond_vet[1]))]
#                 )
#             exec(
#                 f'bc_point_actor_cone_{num_bc}, bc_point_actor_tdof_{num_bc} = view_bondcond_point(previewset["coord"],bondCond_vet,scala_view)'
#             )
#     objs = 0
#     coordMax = [
#         max(previewset["coord"][:, 1]),
#         max(previewset["coord"][:, 2]),
#         max(previewset["coord"][:, 3]),
#     ]
#     # for num_objs in range(len(previewset["regions"])):
#     for num_objs in range(len(previewset["regions"][1])):
#         text = [previewset["regions"][0], str(previewset["regions"][1][num_objs][0])]
#         coord = [
#             sum(previewset["coord"][previewset["regions"][1][num_objs][1] - 1, 1])
#             / len(previewset["coord"][previewset["regions"][1][num_objs][1] - 1, 1]),
#             sum(previewset["coord"][previewset["regions"][1][num_objs][1] - 1, 2])
#             / len(previewset["coord"][previewset["regions"][1][num_objs][1] - 1, 2]),
#             sum(previewset["coord"][previewset["regions"][1][num_objs][1] - 1, 3])
#             / len(previewset["coord"][previewset["regions"][1][num_objs][1] - 1, 2]),
#         ]
#         exec(
#             f"bc_text_actor_{objs} = view_text_point(coord, coordMax, scala_view, text)"
#         )
#         objs += 1

#     i_node = 0
#     if previewset["RENDER"]["cs"] == True:
#         dimcs = len(previewset["tabcs"]["typSection"])
#         for num_bcs in range(dimcs):
#             inci_bcs = previewset["inci"][
#                 np.where(previewset["inci"][:, 3] == num_bcs + 1), :
#             ][0]
#             typSection = int(previewset["tabcs"]["typSection"][num_bcs])
#             dimSection = previewset["tabcs"]["dimSection"][num_bcs]
#             for node_bcs in range(len(inci_bcs)):
#                 noi = int(inci_bcs[node_bcs, 4])
#                 noj = int(inci_bcs[node_bcs, 5])
#                 noix = previewset["coord"][noi - 1, 1]
#                 noiy = previewset["coord"][noi - 1, 2]
#                 noiz = previewset["coord"][noi - 1, 3]
#                 nojx = previewset["coord"][noj - 1, 1]
#                 nojy = previewset["coord"][noj - 1, 2]
#                 nojz = previewset["coord"][noj - 1, 3]
#                 coord_bcs = [noix, noiy, noiz, nojx, nojy, nojz]
#                 exec(
#                     f"beam_extrude_actor_{i_node} = view_beam_crossSection(dimSection, typSection, coord_bcs, scala_view)"
#                 )
#                 i_node += 1
    
#     actor = vtk.vtkActor()
#     actor.SetMapper(mapper)
#     actor.GetProperty().EdgeVisibilityOn()
#     if previewset["RENDER"]["lines"] == False:
#         actor.GetProperty().EdgeVisibilityOff()
#     else:
#         actor.GetProperty().SetLineWidth(1.0)
    
#     text_logo = vtk.vtkTextActor()
#     text_logo.SetInput(
#         "MYFEMPY " + get_version() + 
#         '\n> press "w" to wireframe view \n> press "s" to surface view \n> press "r" to reset view \n> press "q" to exit and continue'
#     )
#     txtprop = text_logo.GetTextProperty()
#     txtprop.SetFontFamilyToArial()
#     txtprop.SetFontSize(20)
#     txtprop.SetColor(1, 1, 1)
#     text_logo.SetDisplayPosition(20, 370)
#     colors = vtk.vtkNamedColors()
#     # backgroundColor = colors.GetColor3d("DarkSlateGray")
#     # actorColor = colors.GetColor3d("Tomato")
#     axis1Color = colors.GetColor3d("red")
#     axis2Color = colors.GetColor3d("red")
#     axis3Color = colors.GetColor3d("red")
#     cubeAxesActor = vtk.vtkCubeAxesActor()
#     cubeAxesActor.SetBounds(actor.GetBounds())
#     cubeAxesActor.SetCamera(renderer.GetActiveCamera())
#     cubeAxesActor.GetTitleTextProperty(0).SetFontSize(25)
#     cubeAxesActor.GetLabelTextProperty(0).SetColor(axis1Color)
#     cubeAxesActor.GetLabelTextProperty(1).SetColor(axis2Color)
#     cubeAxesActor.GetLabelTextProperty(2).SetColor(axis3Color)
#     renderer.AddActor(actor)
#     renderer.AddActor(text_logo)
#     renderer.AddActor(cubeAxesActor)
#     for ipt in range(dimbclist):
#         exec(f"renderer.AddActor(bc_point_actor_cone_{ipt})")
#         exec(f"renderer.AddActor(bc_point_actor_tdof_{ipt})")
#     for bcs in range(i_node):
#         exec(f"renderer.AddActor(beam_extrude_actor_{bcs})")
#     for ff in range(dimfrlist):
#         exec(f"renderer.AddActor(fr_point_actor_cone1_{ff})")
#         exec(f"renderer.AddActor(fr_point_actor_cone2_{ff})")
#     if "plottags" in previewset["RENDER"].keys():
#         for txt in range(objs):
#             exec(f"renderer.AddActor(bc_text_actor_{txt})")
#     else:
#         pass
#     renderer_window.AddRenderer(renderer)
#     renderer.ResetCamera()
#     if previewset["RENDER"]["savepng"] == False:
#         pass
#     else:
#         im = vtk.vtkWindowToImageFilter()
#         writer = (
#             vtk.vtkPNGWriter()
#         )  # vtkSTLWriter()#vtkVRMLExporter()#vtkPolyDataWriter()#
#         im.SetInput(renderer_window)
#         im.Update()
#         writer.SetInputConnection(im.GetOutputPort())
#         writer.SetFileName((path + "/" + previewset["RENDER"]["filename"] + ".png"))
#         writer.Write()
#     if previewset["RENDER"]["show"] == False:
#         pass
#     else:
#         interactor = vtk.vtkRenderWindowInteractor()
#         interactor.SetRenderWindow(renderer_window)
#         interactor.RemoveObservers("RightButtonPressEvent")
#         interactor.Initialize()
#         renderer_window.Render()
#         interactor.Start()
#         interactor.GetRenderWindow().Finalize()


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
